
from __future__ import annotations
import shlex
from itertools import product
from pathlib import Path
import sys



def run_cli(
    *,
    args,
    collect_linkers,
    collect_nodes,
    assemble,
    run_xtb_for_out_dir,
    normalize_shapes,
    run_obabel_for_out_dir=None) -> None:
    
    nodes_dir   = Path(args.nodes).expanduser().resolve()
    linkers_dir = Path(args.linkers).expanduser().resolve()
    out_root    = Path(args.out).expanduser().resolve()
    shapes      = normalize_shapes(args.shape)

    if not nodes_dir.is_dir():
        raise NotADirectoryError(f"nodes not found: {nodes_dir}")
    if not linkers_dir.is_dir():
        raise NotADirectoryError(f"linkers not found: {linkers_dir}")

    linker_json = collect_linkers(linkers_dir, out_root)
    node_json   = collect_nodes(nodes_dir, out_root)

    for nj, lj in product(node_json, linker_json):
        for shape in shapes:
            out_dir = out_root / "assemblies"
            assemble(nj, lj, out_dir, shape)

            for f in out_dir.glob("[!.]*"):
                if not f.is_file():
                    continue
                dest = out_dir / f.stem
                dest.mkdir(parents=True, exist_ok=True)
                f.replace(dest / f.name)

    ### optimization stage(s)
    xtb_flags: list[str] = []
    for chunk in (getattr(args, "xtb_flags", []) or []):
        xtb_flags += shlex.split(chunk)

    obabel_flags: list[str] = []
    for chunk in (getattr(args, "obabel_flags", []) or []):
        obabel_flags += shlex.split(chunk)

    use_xtb = bool(getattr(args, "xtb_opt", False))
    use_ff  = bool(getattr(args, "obabel_opt", False))

    if use_xtb:
        ### If both flags are set, do FF preopt inside the xTB driver
        run_xtb_for_out_dir(
            out_dir,
            extra_flags=xtb_flags,
            ff_preopt=use_ff,
            ff_flags=obabel_flags if use_ff else None,
        )

    ### Pure FF-only mode: run Open Babel optimization only if xTB is NOT requested
    if use_ff and not use_xtb and run_obabel_for_out_dir is not None:
        run_obabel_for_out_dir(
            out_dir,
            extra_flags=obabel_flags,
        )
        