#!/usr/bin/env python3

"""
Uses the Open Babel Python API (openbabel-wheel) instead of system calling obabel/obminimize

Usage:

  obabel_opt.py --path INPUT.mol [--out OUTDIR]
                [--ff uff|mmff94|mmff94s|gaff|ghemical]
                [--algo cg|sd] [--steps N] [--crit THRESH]
                [--keep-work] [--quiet]

Defaults:
  ff=gaff  algo=cg  steps=5000  crit=1e-8

Outputs (in OUTDIR):
  ff_opt.mol   final optimized structure
  ff_opt.log   log (settings, energies, status)

On failure:
  ff_failed_last.mol is written to OUTDIR (copy of input), and exit code is 0
"""



import argparse
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path
try:
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob



def parse_args(argv=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--path", required=True, dest="input_path")
    parser.add_argument("--out", dest="outdir", default=None)
    parser.add_argument("--ff", dest="ff", default="gaff",
                        choices=["uff", "mmff94", "mmff94s", "gaff", "ghemical"])
    parser.add_argument("--algo", dest="algo", default="cg", choices=["cg", "sd"])
    parser.add_argument("--steps", dest="steps", type=int, default=5000)
    parser.add_argument("--crit", dest="crit", default="1e-8")
    parser.add_argument("--keep-work", dest="keep_work", action="store_true", default=False)
    parser.add_argument("--quiet", dest="quiet", action="store_true", default=False)

    parser.add_argument("--obmin-bin", dest="obmin_bin", default=None)

    parser.add_argument("-h", "--help", action="help", help="show this help message and exit")

    args = parser.parse_args(argv)

    if args.steps <= 0:
        parser.error("--steps must be a positive integer")

    try:
        args.crit = float(args.crit)
        if not (args.crit > 0.0 and math.isfinite(args.crit)):
            parser.error("--crit must be a positive finite float")
    except ValueError:
        parser.error("--crit must be a positive finite float")

    return args



def get_in_format(input_path: Path) -> str:
    ext = input_path.suffix.lower().lstrip(".")
    return ext if ext else "mol"



def get_forcefield(ff_name: str) -> ob.OBForceField:
    mapping = {
        "uff": "UFF",
        "mmff94": "MMFF94",
        "mmff94s": "MMFF94s",
        "gaff": "GAFF",
        "ghemical": "Ghemical",
    }
    ob_name = mapping.get(ff_name.lower())
    if ob_name is None:
        raise RuntimeError(f"Unsupported force field: {ff_name}")

    ff = ob.OBForceField.FindForceField(ob_name)
    if ff is None:
        raise RuntimeError(f"Open Babel could not find force field '{ob_name}'")
    return ff



def optimize_with_openbabel(
    input_path: Path,
    outdir: Path,
    ff_name: str,
    algo: str,
    steps: int,
    crit: float,
    keep_work: bool,
    quiet: bool) -> bool:

    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "ff_opt.log"
    opt_mol_path = outdir / "ff_opt.mol"

    ### Scratch dir
    workdir = Path(
        tempfile.mkdtemp(prefix="obmin_run.", dir=str(outdir))
    )

    try:
        ### Copy input into scratch
        work_input = workdir / "input.mol"
        shutil.copy2(input_path, work_input)

        ### Open log
        with log_path.open("w") as log:
            log.write("BACKEND=python-openbabel\n")
            log.write(f"INPUT={input_path}\n")
            log.write(f"WORKDIR={workdir}\n")
            log.write(f"FF={ff_name}\n")
            log.write(f"ALGO={algo}\n")
            log.write(f"STEPS={steps}\n")
            log.write(f"CRIT={crit}\n")

            ### Set up conversion
            in_fmt = get_in_format(input_path)
            conv = ob.OBConversion()
            if not conv.SetInFormat(in_fmt):
                raise RuntimeError(f"Could not set input format '{in_fmt}'")

            if not conv.SetOutFormat("mol"):
                raise RuntimeError("Could not set output format 'mol'")

            mol = ob.OBMol()
            if not conv.ReadFile(mol, str(work_input)):
                raise RuntimeError(f"Failed to read input file {work_input} as {in_fmt}")

            ### Force field
            ff = get_forcefield(ff_name)
            if not ff.Setup(mol):
                raise RuntimeError("Force field Setup(mol) failed")

            initial_energy = ff.Energy()
            log.write(f"ENERGY_INITIAL={initial_energy:.10f}\n")

            ### Optimization
            if algo == "cg":
                log.write("CMD: python-openbabel ConjugateGradients\n")
                n_steps = ff.ConjugateGradients(steps, crit)
            elif algo == "sd":
                log.write("CMD: python-openbabel SteepestDescent\n")
                n_steps = ff.SteepestDescent(steps, crit)
            else:
                raise RuntimeError(f"Unknown optimization algorithm: {algo}")

            ### Apply final coordinates back to mol
            ff.GetCoordinates(mol)
            final_energy = ff.Energy()

            log.write(f"ENERGY_FINAL={final_energy:.10f}\n")
            log.write(f"STEPS_TAKEN={n_steps}\n")

            ### Write output mol
            if not conv.WriteFile(mol, str(opt_mol_path)):
                raise RuntimeError(f"Failed to write output mol {opt_mol_path}")

            ### Quick sanity check
            if not opt_mol_path.is_file() or opt_mol_path.stat().st_size == 0:
                raise RuntimeError("ff_opt.mol was not written or is empty")

        if not quiet:
            print(f"\nOptimization finished. Output: {opt_mol_path}", file=sys.stderr)

        return True

    except Exception as e:
        ### Log error
        with log_path.open("a") as log:
            log.write(f"ERROR: {e}\n")

        if not quiet:
            print(f"Minimization failed: {e}", file=sys.stderr)

        return False

    finally:
        if not keep_work and workdir.is_dir():
            shutil.rmtree(workdir, ignore_errors=True)



def main(argv=None) -> int:
    args = parse_args(argv)

    input_path = Path(args.input_path).resolve()
    if not input_path.is_file():
        print(f"Input not found: {input_path}", file=sys.stderr)
        return 2

    if args.outdir is None:
        outdir = input_path.parent.resolve()
    else:
        outdir = Path(args.outdir).resolve()

    success = optimize_with_openbabel(
        input_path=input_path,
        outdir=outdir,
        ff_name=args.ff,
        algo=args.algo,
        steps=args.steps,
        crit=args.crit,
        keep_work=args.keep_work,
        quiet=args.quiet,
    )

    ### On failure, write ff_failed_last.mol and return 0
    if not success:
        failed_copy = outdir / "ff_failed_last.mol"
        try:
            shutil.copy2(input_path, failed_copy)
        except Exception:
            pass
        return 0

    return 0



if __name__ == "__main__":
    raise SystemExit(main())
