
__version__ = "0.0.5"



import argparse
import runpy
import inspect
import re
import subprocess
import time
import shutil
import shlex
from itertools import product
from pathlib import Path
from contextlib import suppress
import sys, os
import json



### PATH load saveguard (was required in an old version, may be removed in the future)
HERE = Path(__file__).resolve().parent



### Fail more gracefully
def _first_existing_parent(p: Path) -> Path | None:
    p = p.resolve()
    for anc in [p, *p.parents]:
        if anc.exists():
            return anc
    return None

def _hint_nearby(p: Path, max_items: int = 12) -> str:
    parent = _first_existing_parent(p)
    if parent is None:
        return "No part of this path exists."
    if parent.is_file():
        return f"Nearest existing path is a file: {parent}"
    try:
        items = sorted([x.name for x in parent.iterdir()])
    except Exception:
        return f"Nearest existing directory: {parent} (cannot list contents)"
    if not items:
        return f"Nearest existing directory is empty: {parent}"
    shown = items[:max_items]
    more = "" if len(items) <= max_items else f" (and {len(items) - max_items} more)"
    return f"Nearest existing directory: {parent}\nContents: {', '.join(shown)}{more}"

def _require_dir(arg_value: str, label: str) -> tuple[bool, Path | None]:
    if not arg_value:
        print(f"[Error!] Missing required argument: --{label}", file=sys.stderr)
        return False, None

    p = Path(arg_value).expanduser().resolve()

    if not p.exists():
        print(f"[Error!] {label} directory not found: {p}", file=sys.stderr)
        return False, None

    if not p.is_dir():
        print(f"[Error!] {label} is not a directory: {p}", file=sys.stderr)
        return False, None

    return True, p



### calculation imports
XTB_OPT_SH = shutil.which("xtb_opt.sh") or str(HERE / "calculation" / "xtb_opt.sh")
from calculation import obabel_opt as OBABEL_OPT_SH

### utility imports
from utility_functions import json_outline_creator as JSON_OUTLINE
from utility_functions import coordinate_adder as COORD_ADDER
from utility_functions import CBU_to_GBU as CBU2GBU

### assembler imports
from cage_assemblers import cuboctaeder as CUBO_ASM
from cage_assemblers import octaeder as OCTA_ASM
from cage_assemblers import cuboc_d3h as CUBO_D3H_ASM
from cage_assemblers import lantern as LANT_ASM

### interfaces imports
from interfaces.menu_mode import run_menu
from interfaces.cli_mode import run_cli



### shape definition
SHAPES = {
    "cuboctaeder": CUBO_ASM,
    "octaeder": OCTA_ASM,
    "cuboctaeder_d3h": CUBO_D3H_ASM,
    "lantern" : LANT_ASM
}



### Information read-in from filenames
Q_TOKEN = re.compile(r"Q([+-]?\d+)", re.IGNORECASE)
STOICH_M = re.compile(r"M(\d+)", re.IGNORECASE)
STOICH_L = re.compile(r"L(\d+)", re.IGNORECASE)



def _load_json_lenient(path: Path):
    txt = path.read_text(encoding="utf-8")
    try:
        return json.loads(txt)
    except json.JSONDecodeError:
        return json.loads(re.sub(r",(\s*[}\]])", r"\1", txt))



def _find_project_root(start: Path) -> Path | None:
    for anc in [start, *start.parents]:
        if (anc / "nodes_json").is_dir() or (anc / "linkers_json").is_dir():
            return anc
    return None



def _parse_ml_from_names(struct_path: Path) -> tuple[int | None, int | None]:
    ### Try stem first
    stem = struct_path.stem
    m = int(STOICH_M.search(stem).group(1)) if STOICH_M.search(stem) else None
    l = int(STOICH_L.search(stem).group(1)) if STOICH_L.search(stem) else None
    ### climb parents if missing
    if m is None:
        for anc in [struct_path.parent, *struct_path.parents]:
            mm = STOICH_M.search(anc.name)
            if mm:
                m = int(mm.group(1))
                break
    if l is None:
        for anc in [struct_path.parent, *struct_path.parents]:
            ll = STOICH_L.search(anc.name)
            if ll:
                l = int(ll.group(1))
                break
    ### fallback for M: reuse metal counter but neutral factor=1
    if m is None:
        m = infer_charge_from_path(struct_path, factor=1)
    return m, l



### charge detection section

def unit_charge_from_json(path: Path) -> int | None:
    try:
        data = _load_json_lenient(path)
    except Exception:
        return None
    ch = (data.get("composition") or {}).get("charge")
    try:
        return None if ch is None else int(round(float(ch)))
    except Exception:
        return None

def infer_charge_from_path(p: Path, factor: int = 2) -> int | None:
    ### Return 2x metal-count where metals is parsed from M# tokens in file/dir names (only for now, change to read from node json)
    names = [p.stem if hasattr(p, "stem") else p.name, p.name] + [anc.name for anc in p.parents]
    metals = []
    for name in names:
        for m in STOICH_M.finditer(name):
            try:
                metals.append(int(m.group(1)))
            except ValueError:
                pass
    return factor * max(metals) if metals else None

def guess_total_charge_from_metadata(struct_path: Path) -> int | None:
    ### Q token (charge) directly in file name (e.g. Q2 for +2)
    mQ = Q_TOKEN.search(struct_path.stem)
    if mQ:
        return int(mQ.group(1))

    ### Sidecar meta JSON: same stem, .meta.json suffix
    meta = struct_path.with_suffix(".meta.json")
    if meta.exists():
        try:
            data = _load_json_lenient(meta)
            q = (
                data.get("total_charge_suggested")
                or (data.get("charges") or {}).get("total")
            )
            if q is not None:
                return int(round(float(q)))

            ### fallback: compute from fields in meta
            M = int((data.get("stoichiometry") or {}).get("M"))
            L = int((data.get("stoichiometry") or {}).get("L"))
            ch = data.get("charges") or {}
            q_node  = int(round(float(ch.get("node_unit_charge", 0))))
            q_link0 = int(round(float(ch.get("linker_unit_charge_intrinsic", 0))))
            n_dep   = int(round(float(ch.get("linker_deprotonations_per_linker", 0))))
            return M * q_node + L * (q_link0 - n_dep)
        except Exception:
            pass

    ### walk up to project root and use node/linker JSONs
    root = _find_project_root(struct_path.parent)
    if root is None:
        return None

    ### Expected pattern: NODE__LINKER__...<ext>
    parts = struct_path.stem.split("__")
    node_name   = parts[0] if len(parts) > 0 else None
    linker_name = parts[1] if len(parts) > 1 else None

    node_q = linker_q = None
    if node_name:
        p = root / "nodes_json" / f"{node_name}.json"
        if p.exists():
            node_q = unit_charge_from_json(p)
    if linker_name:
        p = root / "linkers_json" / f"{linker_name}.json"
        if p.exists():
            linker_q = unit_charge_from_json(p)

    M, L = _parse_ml_from_names(struct_path)

    total = None
    if (node_q is not None) and (M is not None):
        total = M * node_q
    if (linker_q is not None) and (L is not None):
        total = (total or 0) + L * linker_q
    return total



def run(cmd, cwd=None):
    if not cmd:
        raise ValueError("empty command")
    target, *args = cmd
    args = [str(a) for a in args]

    saved_argv = sys.argv
    saved_cwd = Path.cwd()
    try:
        if cwd is not None:
            os.chdir(cwd)

        if hasattr(target, "main"):
            sys.argv = [getattr(target, "__file__", getattr(target, "__name__", "module")), *args]

            if hasattr(target, "parse_args"):
                ns = target.parse_args()
                try:
                    n_params = len(inspect.signature(target.main).parameters)
                except (ValueError, TypeError):
                    n_params = 0

                if n_params == 1:
                    return target.main(ns)
                else:
                    setattr(target, "args", ns)
                    return target.main()

            return target.main()

        script_path = Path(target)
        if script_path.suffix.lower() != ".py":
            raise ValueError(f"expected a .py script, got {script_path}")
        sys.argv = [str(script_path), *args]
        runpy.run_path(str(script_path), run_name="__main__")

    finally:
        sys.argv = saved_argv
        with suppress(FileNotFoundError, NotADirectoryError, PermissionError):
            os.chdir(saved_cwd)



### [future] here maybe split the next 4 functions to separate submodule to add more than xyz/json input of linkers
### [future] functions effected: xyz_to_json_pipeline and mol_to_json_pipeline, collect_linkers, collect_nodes
### [future] add compatibility for crystallographic data here

def xyz_to_json_pipeline(xyz_path: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{xyz_path.stem}.json"
    run([COORD_ADDER, xyz_path, json_path])
    run([CBU2GBU, json_path, "--write"])
    return json_path

def mol_to_json_pipeline(mol_path: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{mol_path.stem}.json"
    run([COORD_ADDER, mol_path, json_path])
    run([CBU2GBU, json_path, "--write"])
    return json_path

def collect_linkers(linkers_dir: Path, out_root: Path) -> list[Path]:
    ### Collect linker JSONs, converting from XYZ or mol if needed
    ### If more filetypes are added, respective for loops have to be appended here (similar to for xyz in ...)
    dest = out_root / "linkers_json"
    dest.mkdir(parents=True, exist_ok=True)

    src_jsons = sorted(p for p in linkers_dir.glob("*.json") if p.is_file())
    src_xyz   = sorted(p for p in linkers_dir.glob("*.xyz")  if p.is_file())
    src_mol  = sorted(p for p in linkers_dir.glob("*.mol") if p.is_file())

    if not src_jsons and not src_xyz and not src_mol:
        existing = sorted(p for p in dest.glob("*.json") if p.is_file())
        if existing:
            return existing
        raise FileNotFoundError(
            f"No linker .json, .xyz or .mol files in {linkers_dir} and none in {dest}"
        )

    out_jsons: list[Path] = []

    ### copy existing JSONs
    for src in src_jsons:
        dst = dest / src.name
        dst.write_bytes(src.read_bytes())
        out_jsons.append(dst)

    ### convert XYZ to JSON if not already present
    for xyz in src_xyz:
        name = f"{xyz.stem}.json"
        if (linkers_dir / name).is_file() or (dest / name).is_file():
            out_jsons.append(dest / name)
            continue
        out_jsons.append(xyz_to_json_pipeline(xyz, dest))

    ### convert mol to JSON if not already present
    for mol in src_mol:
        name = f"{mol.stem}.json"
        if (linkers_dir / name).is_file() or (dest / name).is_file():
            out_jsons.append(dest / name)
            continue
        out_jsons.append(mol_to_json_pipeline(mol, dest))
    
    print("\n")

    ### deduplicate
    out_jsons = sorted({p.resolve() for p in out_jsons if p.exists()})
    if not out_jsons:
        existing = sorted(p for p in dest.glob("*.json") if p.is_file())
        if existing:
            return existing
        raise FileNotFoundError(f"No linker JSONs available in {dest}")
    return out_jsons

def collect_nodes(nodes_dir: Path, out_root: Path) -> list[Path]:
    dest = out_root / "nodes_json"
    dest.mkdir(parents=True, exist_ok=True)

    out_jsons: list[Path] = []

    for src in sorted(p for p in nodes_dir.glob("*.json") if p.is_file()):
        dst = dest / src.name
        dst.write_bytes(src.read_bytes())
        out_jsons.append(dst)
    if not out_jsons:
        raise FileNotFoundError(f"No node .json files in {nodes_dir}")

    return sorted(out_jsons)



def assemble(node_json: Path, linker_json: Path, out_dir: Path, shape: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    asm = SHAPES.get(shape)
    if asm is None:
        raise ValueError(f"unknown shape: {shape}")
    run([asm, node_json, linker_json, out_dir])



def normalize_shapes(arg) -> list[str]:
    if not arg:
        return list(SHAPES.keys())
    if isinstance(arg, list):
        toks = []
        for v in arg:
            toks += [t.strip() for t in str(v).split(",") if t.strip()]
    else:
        toks = [t.strip() for t in str(arg).split(",") if t.strip()]
    invalid = [t for t in toks if t not in SHAPES]
    if invalid:
        allowed = ", ".join(SHAPES.keys())
        raise ValueError(f"invalid shape(s): {', '.join(invalid)}; allowed: {allowed}")
    return toks



### [future] here maybe split the next 2 functions to separate submodule to
### [future] functions effected: run_xtb_for_out_dir and run_obabel_for_out_dir
### [future] generally is regarded towards calculation submodules (not only opt)

def run_xtb_for_out_dir(
    out_dir: Path,
    *,
    extra_flags: list[str] | None = None,
    ff_preopt: bool = False,
    ff_flags: list[str] | None = None) -> None:

    if extra_flags is None:
        extra_flags = []
    if ff_flags is None:
        ff_flags = []

    print("\n")

    for sub in out_dir.iterdir():
        if not sub.is_dir():
            continue
        
        ### If this folder already has a successful xTB result, skip it entirely
        xtb_out = sub / "xtb_opt.mol"
        if xtb_out.exists() and xtb_out.stat().st_size > 0:
            print(f"\n[xTB] skip folder (xtb_opt.mol exists): \n   cwd={sub}")
            continue

        ### Charge fallback inferred once per subdir
        fallback_charge = infer_charge_from_path(sub)

        ### Iterate over "base" mol files only, skip any previous FF/xTB products
        base_mol_files = sorted(
            p
            for p in sub.glob("*.mol")
            if p.name not in {"ff_opt.mol", "xtb_opt.mol"}
        )
        if not base_mol_files:
            continue

        for base_mol in base_mol_files:
            xtb_input: Path = base_mol

            ### Optional FF preoptimization (OpenBabel) on base_mol (removes preexisting files from previous runs; change later)
            if ff_preopt:
                print(f"\n[FF] start (Open Babel): {base_mol.name}  cwd={sub}")
                print("")

                ff_mol = sub / "ff_opt.mol"
                ff_fail = sub / "ff_failed_last.mol"

                ### Reuseing existing result if present
                if ff_mol.exists() and ff_mol.stat().st_size > 0:
                    xtb_input = ff_mol
                    print(f"\n[FF] reuse: {base_mol.name} -> ff_opt.mol (already exists)")
                else:
                    cmd_ff = [
                        sys.executable, "-m", "calculation.obabel_opt",
                        "--path", str(base_mol),
                        "--out", str(sub),
                        *ff_flags,
                    ]
                    res_ff = subprocess.run(cmd_ff, cwd=sub, check=False)

                    if ff_mol.exists() and ff_mol.stat().st_size > 0:
                        xtb_input = ff_mol
                        print(f"\n[FF] done:   {base_mol.name}  -> ff_opt.mol")
                    else:
                        print(
                            f"\n[FF] Open Babel preopt failed in {sub} "
                            f"(no ff_opt.mol; exit={res_ff.returncode})"
                        )

            ### Charge determination (always from the original base_mol)
            auto = guess_total_charge_from_metadata(base_mol)
            if auto is not None:
                chrg = auto
                print("\n" + "=" * 64)
                print(
                    f"\n[xTB] charge from JSON/meta: {chrg}  "
                    f"({base_mol.name})"
                )
            else:
                chrg = (
                    fallback_charge
                    if fallback_charge is not None
                    else infer_charge_from_path(base_mol)
                )
                if chrg is None:
                    print(
                        f"\n[xTB] skip: cannot infer charge for {base_mol}"
                    )
                    continue
                print(
                    f"\n[xTB] fallback charge: {chrg}  ({base_mol.name})"
                )

            ### xTB optimization on xtb_input (base or FF-optimized)

            print("\n" + "=" * 64)
            print(
                f"\n[xTB] start: {xtb_input.name}  "
                f"chrg={chrg}  cwd={sub}"
            )

            cmd = [
                "/bin/bash",
                XTB_OPT_SH,
                "--path", str(xtb_input.resolve()),
                "--out", str(sub.resolve()),
                "--charge", str(chrg),
                "--mult", "1",
                "--save-trj",
                *extra_flags,
            ]
            
            subprocess.run(["/bin/chmod", "+x", XTB_OPT_SH], check=False)
            res = subprocess.run(cmd, cwd=sub, check=False)
            opt_xyz  = sub / "xtb_opt.xyz"
            xtb_fail = sub / "xtb_failed_last.xyz"
            xtb_opt = sub / "xtb_opt.mol"

            if res.returncode == 0 and opt_xyz.exists():
                print(f"\n[xTB] done:   {xtb_input.name}")

            elif xtb_fail.exists():
                print(
                    f"\n[xTB] failed: {xtb_input.name}  "
                    f"(wrote xtb_failed_last.xyz; see {sub/'xtb_opt.log'})"
                )
            else:
                print(
                    f"\n[xTB] failed: {xtb_input.name}  "
                    f"(no xtb_opt.mol or xtb_failed_last.xyz; exit={res.returncode})"
                )
            
            ### Remove xtb_opt.xyz
            if xtb_opt.exists() and xtb_opt.stat().st_size > 0:
                with suppress(FileNotFoundError):
                    opt_xyz.unlink()

def run_obabel_for_out_dir(out_dir: Path, *, extra_flags: list[str] | None = None) -> None:
    print("\n")
    if extra_flags is None:
        extra_flags = []

    for sub in out_dir.iterdir():
        if not sub.is_dir():
            continue

        ### FF-only mode: if this folder already has FF results, skip the folder entirely
        ff_out = sub / "ff_opt.mol"
        if ff_out.exists() and ff_out.stat().st_size > 0:
            print(f"\n[OBABEL] skip folder (ff_opt.mol exists): \n   cwd={sub}")
            continue

        ### Also skip folders that already have xTB results
        xtb_out = sub / "xtb_opt.mol"
        if xtb_out.exists() and xtb_out.stat().st_size > 0:
            print(f"\n[OBABEL] skip folder (xtb_opt.mol exists): \n   cwd={sub}")
            continue

        inputs: list[Path] = []
        inputs.extend(sorted(sub.glob("*.mol")))
        inputs.extend(sorted(sub.glob("*.xyz")))
        if not inputs:
            continue

        for inp in inputs:
            print("\n" + "="*64)
            print(f"\n[OBABEL] start: {inp.name} \n   cwd={sub}")
            print("")
            cmd = [
                sys.executable, "-m", "calculation.obabel_opt",
                "--path", str(inp),
                "--out", str(sub),
                *extra_flags,
            ]
            res = subprocess.run(cmd, cwd=sub, check=False)

            ff_opt  = sub / "ff_opt.mol"
            ff_fail = sub / "ff_failed_last.mol"

            if res.returncode == 0 and ff_opt.exists():
                print(f"\n[OBABEL] done: \n   {inp.name}")
            elif ff_fail.exists():
                print(f"\n[OBABEL] failed: \n   {inp.name}  (wrote ff_failed_last.mol; see {sub/'ff_opt.log'})")
            else:
                print(
                    f"\n[OBABEL] failed: \n   {inp.name}  "
                    f"(no ff_opt.mol or ff_failed_last.mol; exit={res.returncode})"
                )



### Entry Point
def parse_args(argv=None):
    fmt = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=32)
    ap = argparse.ArgumentParser(
        description="Cageinator: batch pipeline from Linkers and Nodes to assembled cages.",
        formatter_class=fmt,
    )

    ap.add_argument("--menu",
                    action="store_true",
                    help="Interactive menu mode. Default: CLI batch mode")
    ap.add_argument("--nodes",
                    help="Folder with node JSON files")
    ap.add_argument("--linkers",
                    help="Folder with linker XYZ/mol files")
    ap.add_argument("--out",
                    help="Destination folder for JSONs and assemblies")
    ap.add_argument("--shape",
                    metavar="SHAPE",
                    nargs="+",
                    help=f"Only build these shapes (comma or space separated). "
                    f"Allowed: {', '.join(SHAPES.keys())}. Omit to build all")
    ap.add_argument("--xtb-opt",
                    action="store_true",
                    help="xTB-Optimization of all built cages")
    ap.add_argument("--xtb-flags",
                    action="append",
                    default=[],
                    metavar='"FLAGS"',
                    help='Extra flags for xtb_opt.sh. Quote as a chunk. Can repeat. Example: --xtb-flags "--method gfnff --threads 8"')
    ap.add_argument("--xtb-help",
                    action="store_true",
                    help="Show optimizer (xtb_opt.sh) help and exit")
    ap.add_argument("--obabel-opt",
                    action="store_true",
                    help="Optimize all built cages using Open Babel (obminimize)")
    ap.add_argument("--obabel-flags",
                    action="append",
                    default=[],
                    metavar='"FLAGS"',
                    help='Extra flags for obabel_opt.sh. Quote as a chunk. Can repeat. '
                        'Example: --obabel-flags "--ff mmff94 --steps 500"')
    ap.add_argument("--version", "-v",
                    action="version",
                    version=f"The Cageinator, Version: {__version__}",
                    help="Print the Version Number of the Cageinator")

    return ap.parse_args(argv)



def main(argv=None) -> int:
    args = parse_args(argv)

    if args.menu:
        run_menu(
            collect_linkers=collect_linkers,
            collect_nodes=collect_nodes,
            assemble=assemble,
            run_xtb_for_out_dir=run_xtb_for_out_dir,
            SHAPES=SHAPES,
        )
        return 0

    if getattr(args, "xtb_help", False):
        subprocess.run(["/bin/chmod", "+x", XTB_OPT_SH], check=False)
        subprocess.run([XTB_OPT_SH, "--help"], check=False)
        return 0

    if not (args.nodes and args.linkers and args.out):
        print(
            "Usage: cageinator.py --nodes NODES_DIR --linkers LINKERS_DIR --out OUT_DIR",
            file=sys.stderr,
        )
        return 2

    ok_nodes, nodes_dir = _require_dir(args.nodes, "nodes")
    ok_linkers, linkers_dir = _require_dir(args.linkers, "linkers")

    out_dir = Path(args.out).expanduser().resolve()
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"[ERROR] Cannot create/access out directory:\n  {out_dir}\nReason: {e}", file=sys.stderr)
        return 2

    if not (ok_nodes and ok_linkers):
        return 2

    try:
        run_cli(
            args=args,
            collect_linkers=collect_linkers,
            collect_nodes=collect_nodes,
            assemble=assemble,
            run_xtb_for_out_dir=run_xtb_for_out_dir,
            normalize_shapes=normalize_shapes,
            run_obabel_for_out_dir=run_obabel_for_out_dir,
        )
        return 0

    except KeyboardInterrupt:
        return 130

    except (NotADirectoryError, FileNotFoundError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2



if __name__ == "__main__":
    raise SystemExit(main())




