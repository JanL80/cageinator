#!/usr/bin/env python3

import os
import re
import sys
import shutil
import signal
import tempfile
import subprocess
from pathlib import Path



def usage(METHOD, OPTLEVEL, MAXITER, SOLVMODEL, CHARGE, MULT, VERBOSE, KEEP_WORK):
    sys.stdout.write(
        f"""Usage:
  xtb_opt.sh --path INPUT.mol [--out OUTDIR] [--charge Q] [--mult M]
             [--method gfn2|gfn1|gfnff] [--opt loose|normal|tight|vtight]
             [--maxiter N] [--solvent NAME] [--solvent-model alpb|gbsa]
             [--threads N] [--xtb-bin /path/to/xtb]
             [--namespace NAME] [--save-trj] [--keep-work] [--quiet]

Defaults:
  method={METHOD}  opt={OPTLEVEL}  maxiter={MAXITER}  solvent_model={SOLVMODEL}
  charge={CHARGE}  mult={MULT}  verbose={VERBOSE}  keep_work={KEEP_WORK}

Outputs (in OUTDIR):
  xtb_opt.mol          final optimized structure (MOL)
  xtb_opt.log          xTB stdout/stderr log
  xtb_opt.trj          trajectory if available (XYZ frames) (with --save-trj)
  xtb_failed_last.mol  input MOL copied if xtb fails

Examples:
  xtb_opt.sh --path a.mol --out ./run --charge 2 --mult 1 --method gfn2 --save-trj
  xtb_opt.sh --path b.mol --method gfnff --opt normal --threads 8 --quiet
"""
    )



def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]



    ### defaults
    INPUT = ""
    OUTDIR = ""
    CHARGE = 0
    MULT = 1
    METHOD = "gfnff"
    OPTLEVEL = "tight"
    MAXITER = 1000
    SOLVENT = ""
    SOLVMODEL = "gbsa"
    THREADS = ""
    XTB_BIN = os.environ.get("XTB_BIN", "")
    KEEP_WORK = 0
    SAVE_TRJ = 0
    VERBOSE = 1
    NAMESPACE = ""



    ### parse
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--path":
            INPUT = argv[i + 1]
            i += 2
        elif a == "--out":
            OUTDIR = argv[i + 1]
            i += 2
        elif a == "--charge":
            CHARGE = argv[i + 1]
            i += 2
        elif a == "--mult":
            MULT = argv[i + 1]
            i += 2
        elif a == "--method":
            METHOD = argv[i + 1]
            i += 2
        elif a == "--opt":
            OPTLEVEL = argv[i + 1]
            i += 2
        elif a == "--maxiter":
            MAXITER = argv[i + 1]
            i += 2
        elif a == "--solvent":
            SOLVENT = argv[i + 1]
            i += 2
        elif a == "--solvent-model":
            SOLVMODEL = argv[i + 1]
            i += 2
        elif a == "--threads":
            THREADS = argv[i + 1]
            i += 2
        elif a == "--xtb-bin":
            XTB_BIN = argv[i + 1]
            i += 2
        elif a == "--namespace":
            NAMESPACE = argv[i + 1]
            i += 2
        elif a == "--save-trj":
            SAVE_TRJ = 1
            i += 1
        elif a == "--keep-work":
            KEEP_WORK = 1
            i += 1
        elif a == "--quiet":
            VERBOSE = 0
            i += 1
        elif a in ("-h", "--help"):
            usage(METHOD, OPTLEVEL, MAXITER, SOLVMODEL, CHARGE, MULT, VERBOSE, KEEP_WORK)
            return 0
        else:
            sys.stderr.write(f"Unknown arg: {a}\n")
            usage(METHOD, OPTLEVEL, MAXITER, SOLVMODEL, CHARGE, MULT, VERBOSE, KEEP_WORK)
            return 2

    if not INPUT:
        sys.stdout.write("Missing --path\n")
        usage(METHOD, OPTLEVEL, MAXITER, SOLVMODEL, CHARGE, MULT, VERBOSE, KEEP_WORK)
        return 2

    in_path = Path(INPUT)
    if not in_path.is_file():
        sys.stderr.write(f"Input not found: {INPUT}\n")
        return 2



    ### Open Babel conversion via Python API (generic)
    def ob_convert_py(infile, outfile, infmt, outfmt):
        infile = str(infile)
        outfile = str(outfile)
        infmt = str(infmt)
        outfmt = str(outfmt)

        try:
            from openbabel import openbabel as ob
        except Exception:
            import openbabel as ob

        conv = ob.OBConversion()
        if not conv.SetInFormat(infmt):
            sys.stderr.write(f"OpenBabel: cannot set input format '{infmt}'\n")
            raise SystemExit(2)

        if outfmt.lower() in {"mdl", "sdf", "sd"}:
            outfmt = "mol"

        if not conv.SetOutFormat(outfmt):
            sys.stderr.write(f"OpenBabel: cannot set output format '{outfmt}'\n")
            raise SystemExit(2)

        ### Force V3000 CTAB
        conv.AddOption("3", ob.OBConversion.OUTOPTIONS)

        mol = ob.OBMol()
        if not conv.ReadFile(mol, infile):
            sys.stderr.write(f"OpenBabel: failed to read '{infile}' as '{infmt}'\n")
            raise SystemExit(2)

        if not conv.WriteFile(mol, outfile):
            sys.stderr.write(f"OpenBabel: failed to write '{outfile}' as '{outfmt}'\n")
            raise SystemExit(2)

        return 0



    ### xtb binary
    if not XTB_BIN:
        found = shutil.which("xtb")
        if not found:
            sys.stderr.write("xtb not found in PATH and --xtb-bin not given\n")
            return 127
        XTB_BIN = found

    if not (Path(XTB_BIN).is_file() and os.access(XTB_BIN, os.X_OK)):
        sys.stderr.write(f"xtb not executable: {XTB_BIN}\n")
        return 126



    ### out dir
    if not OUTDIR:
        OUTDIR = str(in_path.parent.resolve())
    out_dir = Path(OUTDIR)
    out_dir.mkdir(parents=True, exist_ok=True)



    ### threads
    if THREADS:
        os.environ["OMP_NUM_THREADS"] = str(THREADS)



    ### opt level check
    if OPTLEVEL not in ("loose", "normal", "tight", "vtight"):
        sys.stderr.write(f"Unknown --opt {OPTLEVEL}\n")
        return 2



    ### multiplicity UHF
    if not re.fullmatch(r"[0-9]+", str(MULT)):
        sys.stderr.write("--mult must be integer\n")
        return 2
    MULT_int = int(MULT)
    if MULT_int > 0:
        UHF = MULT_int - 1
    else:
        UHF = 0



    ### scratch
    WORKDIR = Path(tempfile.mkdtemp(prefix="xtb_run.", dir=str(out_dir)))
    _work_keep = (KEEP_WORK == 1)

    _exit_due_to_signal = {"code": None}

    def _signal_handler(signum, frame):

        if signum == signal.SIGINT:
            _exit_due_to_signal["code"] = 130
        elif signum == signal.SIGTERM:
            _exit_due_to_signal["code"] = 143
        raise SystemExit(_exit_due_to_signal["code"])

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)

    try:
        ### Copy input MOL verbatim; DO NOT normalize
        shutil.copyfile(str(in_path), str(WORKDIR / "input.mol"))

        os.chdir(str(WORKDIR))


        ### extract last XYZ frame without pre-known natoms
        _float_re = re.compile(r"^[-+]?([0-9]*\.[0-9]+|[0-9]+)([eE][-+]?[0-9]+)?$")

        def extract_last_xyz(infile, out_xyz):
            infile = Path(infile)
            out_xyz = Path(out_xyz)

            def isint_line(s):
                return bool(re.fullmatch(r"\s*[0-9]+\s*", s))

            def isfloat(tok):
                return bool(_float_re.fullmatch(tok))

            def isatomline(line):
                if "*" in line:
                    return False
                line = line.strip()
                if not line:
                    return False
                parts = re.split(r"\s+", line)
                if len(parts) < 4:
                    return False
                if not re.fullmatch(r"[A-Za-z]{1,3}", parts[0]):
                    return False
                return isfloat(parts[1]) and isfloat(parts[2]) and isfloat(parts[3])

            last_buf = None

            try:
                with infile.open("r", encoding="utf-8", errors="replace") as f:
                    lines = f.readlines()
            except Exception:
                return False

            idx = 0
            nlines = len(lines)
            while idx < nlines:
                line0 = lines[idx]
                if isint_line(line0):
                    n = int(line0.strip() or "0")
                    buf = [line0]
                    if idx + 1 >= nlines:
                        idx += 1
                        continue
                    line1 = lines[idx + 1]
                    if "*" in line1:
                        idx += 1
                        continue
                    buf.append(line1)

                    ok = True
                    start = idx + 2
                    end = start + n
                    if end > nlines:
                        ok = False
                    else:
                        for j in range(start, end):
                            if not isatomline(lines[j]):
                                ok = False
                                break
                            buf.append(lines[j])

                    if ok:
                        last_buf = "".join(buf)
                        idx = end
                        continue
                idx += 1

            if last_buf:
                try:
                    tmp = out_xyz.with_suffix(out_xyz.suffix + ".tmp")
                    tmp.write_text(last_buf, encoding="utf-8")
                    if tmp.stat().st_size > 0:
                        tmp.replace(out_xyz)
                        return True
                except Exception:
                    return False
            return False


        ### build argv
        cmd = [str(XTB_BIN), "input.mol", "--opt", str(OPTLEVEL), "--cycles", str(MAXITER)]
        cmd += ["--chrg", str(CHARGE)]

        m = str(METHOD)
        if m == "gfn2":
            cmd += ["--gfn", "2"]
        elif m == "gfn1":
            cmd += ["--gfn", "1"]
        elif m in ("gfnff", "GFNFF", "ff"):
            cmd += ["--gfnff"]
        else:
            sys.stderr.write(f"Unknown --method {METHOD}\n")
            return 2

        if UHF > 0 and m not in ("gfnff", "GFNFF", "ff"):
            cmd += ["--uhf", str(UHF)]

        if SOLVENT:
            if SOLVMODEL == "alpb":
                cmd += ["--alpb", str(SOLVENT)]
            elif SOLVMODEL == "gbsa":
                cmd += ["--gbsa", str(SOLVENT)]
            else:
                sys.stderr.write(f"Unknown --solvent-model {SOLVMODEL}\n")
                return 2

        if NAMESPACE:
            cmd += ["--namespace", str(NAMESPACE)]
        if VERBOSE == 1:
            cmd += ["--verbose"]


        ### run
        OPT_OK = 1
        log_path = out_dir / "xtb_opt.log"

        def _write_cmd_header():
            with log_path.open("w", encoding="utf-8") as lf:
                lf.write("CMD:")
                for a in cmd:
                    lf.write(f" {a}")
                lf.write("\n")

        try:
            if VERBOSE == 1:
                _write_cmd_header()

                with log_path.open("a", encoding="utf-8") as lf:
                    p = subprocess.Popen(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        bufsize=1,
                    )
                    assert p.stdout is not None
                    for line in p.stdout:
                        sys.stdout.write(line)
                        sys.stdout.flush()
                        lf.write(line)
                        lf.flush()
                    rc = p.wait()
                if rc != 0:
                    OPT_OK = 0
            else:
                with log_path.open("w", encoding="utf-8") as lf:
                    p = subprocess.Popen(
                        cmd,
                        stdout=lf,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                    rc = p.wait()
                OPT_OK = 1 if rc == 0 else 0

        except SystemExit:
            raise
        except Exception:
            OPT_OK = 0


        #### on failure: write xtb_failed_last.mol and continue (deprecated, will be removed in the future)
        #if [ "$OPT_OK" -ne 1 ]; then
        #  echo "xtb failed. Extracting last valid frame and converting to xtb_failed_last.mol." >&2
        #
        #  OUT_FAIL_MOL="$OUTDIR/xtb_failed_last.mol"
        #
        #  ### prefer native trajectory first
        #  SRC=""
        #  if [ -n "$NAMESPACE" ] && [ -s "${NAMESPACE}.trj" ]; then
        #    SRC="${NAMESPACE}.trj"
        #  else
        #    for cand in xtbopt.trj xtb.trj md.trj opt.trj *.trj; do
        #      [ -s "$cand" ] || continue
        #      SRC="$cand"
        #      break
        #    done
        #  fi
        #
        #  ### 1) try trajectory -> MOL (filtered)
        #  if [ -n "$SRC" ] && extract_last_mol "$SRC" "$OUT_FAIL_MOL"; then
        #    :
        #  else
        #    ### 2) fall back to log -> MOL (filtered)
        #    if [ -f "xtbopt.log" ]; then
        #      SRC_LOG="xtbopt.log"
        #    else
        #      SRC_LOG="$OUTDIR/xtb_opt.log"
        #    fi
        #
        #    if [ -s "$SRC_LOG" ] && extract_last_mol "$SRC_LOG" "$OUT_FAIL_MOL"; then
        #      :
        #    else
        #      echo "No valid last frame found in trajectory/log; falling back to input.mol" >&2
        #      cp -f "$WORKDIR/input.mol" "$OUT_FAIL_MOL" || true
        #    fi
        #  fi
        #fi


        ### harvest native trajectory of xTB
        if SAVE_TRJ == 1:
            TRJ_OK = 0
            if NAMESPACE and (WORKDIR / f"{NAMESPACE}.trj").is_file() and (WORKDIR / f"{NAMESPACE}.trj").stat().st_size > 0:
                shutil.copyfile(str(WORKDIR / f"{NAMESPACE}.trj"), str(out_dir / "xtb_opt.trj"))
                TRJ_OK = 1
            if TRJ_OK != 1:
                for cand in ["xtbopt.trj", "xtb.trj", "md.trj", "opt.trj"]:
                    pth = WORKDIR / cand
                    if pth.is_file() and pth.stat().st_size > 0:
                        shutil.copyfile(str(pth), str(out_dir / "xtb_opt.trj"))
                        TRJ_OK = 1
                        break
                if TRJ_OK != 1:
                    for pth in WORKDIR.glob("*.trj"):
                        if pth.is_file() and pth.stat().st_size > 0:
                            shutil.copyfile(str(pth), str(out_dir / "xtb_opt.trj"))
                            TRJ_OK = 1
                            break


        ### writing final structure as xtb_opt.mol
        if OPT_OK == 1:
            SRC_TRJ = ""

            ### Choose trajectory source (WORKDIR)
            if NAMESPACE and (WORKDIR / f"{NAMESPACE}.trj").is_file() and (WORKDIR / f"{NAMESPACE}.trj").stat().st_size > 0:
                SRC_TRJ = str(WORKDIR / f"{NAMESPACE}.trj")
            else:
                for cand in ["xtbopt.trj", "xtb.trj", "md.trj", "opt.trj"]:
                    pth = WORKDIR / cand
                    if pth.is_file() and pth.stat().st_size > 0:
                        SRC_TRJ = str(pth)
                        break
                if not SRC_TRJ:
                    for pth in WORKDIR.glob("*.trj"):
                        if pth.is_file() and pth.stat().st_size > 0:
                            SRC_TRJ = str(pth)
                            break

            if SRC_TRJ and extract_last_xyz(SRC_TRJ, out_dir / "xtb_opt.xyz"):
                pass
            else:
                SRC_LOG = ""
                if (WORKDIR / "xtbopt.log").is_file() and (WORKDIR / "xtbopt.log").stat().st_size > 0:
                    SRC_LOG = str(WORKDIR / "xtbopt.log")
                elif log_path.is_file() and log_path.stat().st_size > 0:
                    SRC_LOG = str(log_path)

                if SRC_LOG and extract_last_xyz(SRC_LOG, out_dir / "xtb_opt.xyz"):
                    pass
                else:
                    sys.stderr.write("[ERROR] xTB succeeded but no xtb_opt.xyz could be extracted from trajectory/log.\n")
                    return 1

        else:
            ### if xTB failed, write xtb_failed_last.xyz
            SRC_TRJ = ""

            if NAMESPACE and (WORKDIR / f"{NAMESPACE}.trj").is_file() and (WORKDIR / f"{NAMESPACE}.trj").stat().st_size > 0:
                SRC_TRJ = str(WORKDIR / f"{NAMESPACE}.trj")
            else:
                for cand in ["xtbopt.trj", "xtb.trj", "md.trj", "opt.trj"]:
                    pth = WORKDIR / cand
                    if pth.is_file() and pth.stat().st_size > 0:
                        SRC_TRJ = str(pth)
                        break
                if not SRC_TRJ:
                    for pth in WORKDIR.glob("*.trj"):
                        if pth.is_file() and pth.stat().st_size > 0:
                            SRC_TRJ = str(pth)
                            break

            if SRC_TRJ and extract_last_xyz(SRC_TRJ, out_dir / "xtb_failed_last.xyz"):
                pass
            else:
                SRC_LOG = ""
                if (WORKDIR / "xtbopt.log").is_file() and (WORKDIR / "xtbopt.log").stat().st_size > 0:
                    SRC_LOG = str(WORKDIR / "xtbopt.log")
                elif log_path.is_file() and log_path.stat().st_size > 0:
                    SRC_LOG = str(log_path)
                if SRC_LOG:
                    try:
                        extract_last_xyz(SRC_LOG, out_dir / "xtb_failed_last.xyz")
                    except Exception:
                        pass


        ### extract optimization trajectory from log if present (XYZ frames)
        if SAVE_TRJ == 1:
            trj_out = out_dir / "xtb_opt.trj"
            if not (trj_out.is_file() and trj_out.stat().st_size > 0):
                SRC_LOG = ""
                if (WORKDIR / "xtbopt.log").is_file():
                    SRC_LOG = str(WORKDIR / "xtbopt.log")
                if not SRC_LOG:
                    SRC_LOG = str(log_path)

                srcp = Path(SRC_LOG)
                if srcp.is_file() and srcp.stat().st_size > 0:
                    def isint_line(s):
                        return bool(re.fullmatch(r"\s*[0-9]+\s*", s))

                    try:
                        with srcp.open("r", encoding="utf-8", errors="replace") as f:
                            lines = f.readlines()
                    except Exception:
                        lines = []

                    out_chunks = []
                    idx = 0
                    nlines = len(lines)
                    while idx < nlines:
                        if isint_line(lines[idx]):
                            try:
                                n = int(lines[idx].strip())
                            except Exception:
                                idx += 1
                                continue
                            buf = [lines[idx]]
                            if idx + 1 >= nlines:
                                idx += 1
                                continue
                            buf.append(lines[idx + 1])
                            ok = True
                            start = idx + 2
                            end = start + n
                            if end > nlines:
                                ok = False
                            else:
                                buf.extend(lines[start:end])
                            if ok:
                                out_chunks.append("".join(buf))
                                idx = end
                                continue
                        idx += 1

                    if out_chunks:
                        try:
                            trj_out.write_text("".join(out_chunks), encoding="utf-8")
                        except Exception:
                            pass

                    if not (trj_out.is_file() and trj_out.stat().st_size > 0):
                        try:
                            trj_out.unlink()
                        except Exception:
                            pass

        return 0

    finally:
        if not _work_keep:
            try:
                shutil.rmtree(str(WORKDIR))
            except Exception:
                pass



if __name__ == "__main__":
    raise SystemExit(main())
    