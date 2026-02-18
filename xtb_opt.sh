#!/bin/sh

set -eu

usage() {
  cat <<USAGE
Usage:
  xtb_opt.sh --path INPUT.mol [--out OUTDIR] [--charge Q] [--mult M]
             [--method gfn2|gfn1|gfnff] [--opt loose|normal|tight|vtight]
             [--maxiter N] [--solvent NAME] [--solvent-model alpb|gbsa]
             [--threads N] [--xtb-bin /path/to/xtb]
             [--namespace NAME] [--save-trj] [--keep-work] [--quiet]

Defaults:
  method=$METHOD  opt=$OPTLEVEL  maxiter=$MAXITER  solvent_model=$SOLVMODEL
  charge=$CHARGE  mult=$MULT  verbose=$VERBOSE  keep_work=$KEEP_WORK

Outputs (in OUTDIR):
  xtb_opt.mol          final optimized structure (MOL)
  xtb_opt.log          xTB stdout/stderr log
  xtb_opt.trj          trajectory if available (XYZ frames) (with --save-trj)
  xtb_failed_last.mol  input MOL copied if xtb fails

Examples:
  xtb_opt.sh --path a.mol --out ./run --charge 2 --mult 1 --method gfn2 --save-trj
  xtb_opt.sh --path b.mol --method gfnff --opt normal --threads 8 --quiet
USAGE
}



### defaults
INPUT=""
OUTDIR=""
CHARGE=0
MULT=1
METHOD="gfnff"
OPTLEVEL="tight"
MAXITER=1000
SOLVENT=""
SOLVMODEL="gbsa"
THREADS=""
XTB_BIN="${XTB_BIN:-}"
KEEP_WORK=0
SAVE_TRJ=0
VERBOSE=1
NAMESPACE=""



### parse
while [ $# -gt 0 ]; do
  case "$1" in
    --path) INPUT=$2; shift 2;;
    --out) OUTDIR=$2; shift 2;;
    --charge) CHARGE=$2; shift 2;;
    --mult) MULT=$2; shift 2;;
    --method) METHOD=$2; shift 2;;
    --opt) OPTLEVEL=$2; shift 2;;
    --maxiter) MAXITER=$2; shift 2;;
    --solvent) SOLVENT=$2; shift 2;;
    --solvent-model) SOLVMODEL=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --xtb-bin) XTB_BIN=$2; shift 2;;
    --namespace) NAMESPACE=$2; shift 2;;
    --save-trj) SAVE_TRJ=1; shift 1;;
    --keep-work) KEEP_WORK=1; shift 1;;
    --quiet) VERBOSE=0; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[ -n "$INPUT" ] || { echo "Missing --path"; usage; exit 2; }
[ -f "$INPUT" ] || { echo "Input not found: $INPUT" >&2; exit 2; }



### Open Babel conversion via Python API (generic)
ob_convert_py() {
  infile="$1"
  outfile="$2"
  infmt="$3"
  outfmt="$4"

  python3 - "$infile" "$outfile" "$infmt" "$outfmt" <<'PY'
import sys

infile, outfile, infmt, outfmt = sys.argv[1:5]

try:
    from openbabel import openbabel as ob
except Exception:
    import openbabel as ob

conv = ob.OBConversion()
if not conv.SetInFormat(infmt):
    sys.stderr.write(f"OpenBabel: cannot set input format '{infmt}'\n")
    sys.exit(2)

# Important: use "mol" as output format and force V3000 writing.
# This avoids V2000 fixed-width coordinate overflow that manifests as "*****".
if outfmt.lower() in {"mdl", "sdf", "sd"}:
    outfmt = "mol"

if not conv.SetOutFormat(outfmt):
    sys.stderr.write(f"OpenBabel: cannot set output format '{outfmt}'\n")
    sys.exit(2)

# Force V3000 CTAB
conv.AddOption("3", ob.OBConversion.OUTOPTIONS)

mol = ob.OBMol()
if not conv.ReadFile(mol, infile):
    sys.stderr.write(f"OpenBabel: failed to read '{infile}' as '{infmt}'\n")
    sys.exit(2)

if not conv.WriteFile(mol, outfile):
    sys.stderr.write(f"OpenBabel: failed to write '{outfile}' as '{outfmt}'\n")
    sys.exit(2)

sys.exit(0)
PY
}



### xtb binary
if [ -z "$XTB_BIN" ]; then
  if ! XTB_BIN="$(command -v xtb 2>/dev/null)"; then
    echo "xtb not found in PATH and --xtb-bin not given" >&2; exit 127
  fi
fi
[ -x "$XTB_BIN" ] || { echo "xtb not executable: $XTB_BIN" >&2; exit 126; }



### out dir
if [ -z "$OUTDIR" ]; then OUTDIR=$(cd "$(dirname "$INPUT")" && pwd); fi
mkdir -p "$OUTDIR"



### threads
if [ -n "$THREADS" ]; then export OMP_NUM_THREADS="$THREADS"; fi



### opt level check
case "$OPTLEVEL" in
  loose|normal|tight|vtight) :;;
  *) echo "Unknown --opt $OPTLEVEL" >&2; exit 2;;
esac



### multiplicity -> UHF
case "$MULT" in ''|*[!0-9]*) echo "--mult must be integer" >&2; exit 2;; esac
if [ "$MULT" -gt 0 ]; then UHF=$((MULT - 1)); else UHF=0; fi



### scratch
WORKDIR=$(mktemp -d "${OUTDIR%/}/xtb_run.XXXXXX")
cleanup() {
  status=$?
  if [ "$KEEP_WORK" -ne 1 ]; then rm -rf "$WORKDIR"; fi
  exit $status
}
trap cleanup EXIT INT TERM



### Copy input MOL verbatim; DO NOT normalize
cp -f "$INPUT" "$WORKDIR/input.mol"

cd "$WORKDIR"



### extract last XYZ frame without pre-known natoms
extract_last_xyz() {
  in="$1"
  out_xyz="$2"
  tmp="$(mktemp)"

  awk '
    function isint(s){ return (s ~ /^[[:space:]]*[0-9]+[[:space:]]*$/) }
    function isfloat(s){
      return (s ~ /^[-+]?([0-9]*\.[0-9]+|[0-9]+)([eE][-+]?[0-9]+)?$/)
    }
    function isatomline(line,   a,n){
      if (line ~ /\*/) return 0
      gsub(/^[[:space:]]+|[[:space:]]+$/,"",line)
      n = split(line, a, /[[:space:]]+/)
      if (n < 4) return 0
      if (a[1] !~ /^[A-Za-z]{1,3}$/) return 0
      return isfloat(a[2]) && isfloat(a[3]) && isfloat(a[4])
    }
    {
      if (isint($0)) {
        n = $1 + 0
        buf = $0 ORS

        if ((getline line) <= 0) next
        if (line ~ /\*/) next
        buf = buf line ORS

        ok = 1
        for (i=0; i<n; i++) {
          if ((getline line) <= 0) { ok=0; break }
          if (!isatomline(line))   { ok=0; break }
          buf = buf line ORS
        }
        if (ok) last = buf
      }
    }
    END { if (last != "") printf "%s", last }
  ' "$in" > "$tmp"

  if [ -s "$tmp" ]; then
    mv -f "$tmp" "$out_xyz"
    return 0
  fi

  rm -f "$tmp"
  return 1
}



### build argv
set -- "$XTB_BIN" "input.mol" "--opt" "$OPTLEVEL" "--cycles" "$MAXITER"

set -- "$@" "--chrg" "$CHARGE"

case "$METHOD" in
  gfn2) set -- "$@" "--gfn" "2" ;;
  gfn1) set -- "$@" "--gfn" "1" ;;
  gfnff|GFNFF|ff) set -- "$@" "--gfnff" ;;
  *) echo "Unknown --method $METHOD" >&2; exit 2;;
esac

if [ "$UHF" -gt 0 ] && [ "$METHOD" != "gfnff" ] && [ "$METHOD" != "GFNFF" ] && [ "$METHOD" != "ff" ]; then
  set -- "$@" "--uhf" "$UHF"
fi

if [ -n "$SOLVENT" ]; then
  case "$SOLVMODEL" in
    alpb) set -- "$@" "--alpb" "$SOLVENT" ;;
    gbsa) set -- "$@" "--gbsa" "$SOLVENT" ;;
    *) echo "Unknown --solvent-model $SOLVMODEL" >&2; exit 2;;
  esac
fi

[ -n "$NAMESPACE" ] && set -- "$@" "--namespace" "$NAMESPACE"
[ "$VERBOSE" -eq 1 ] && set -- "$@" "--verbose"



### run
OPT_OK=1
if [ "$VERBOSE" -eq 1 ]; then
  printf "CMD:" > "$OUTDIR/xtb_opt.log"
  for a in "$@"; do printf " %s" "$a" >> "$OUTDIR/xtb_opt.log"; done
  printf "\n" >> "$OUTDIR/xtb_opt.log"

  PIPE="$WORKDIR/.xtb.pipe"
  rm -f "$PIPE"; mkfifo "$PIPE"
  tee -a "$OUTDIR/xtb_opt.log" < "$PIPE" &
  TEE_PID=$!

  if "$@" >"$PIPE" 2>&1; then
    CMD_STATUS=0
  else
    CMD_STATUS=$?
  fi

  wait "$TEE_PID" 2>/dev/null || true
  rm -f "$PIPE"

  [ "$CMD_STATUS" -eq 0 ] || OPT_OK=0
else
  if "$@" > "$OUTDIR/xtb_opt.log" 2>&1; then
    OPT_OK=1
  else
    OPT_OK=0
  fi
fi



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
if [ "$SAVE_TRJ" -eq 1 ]; then
  TRJ_OK=0
  if [ -n "$NAMESPACE" ] && [ -s "${NAMESPACE}.trj" ]; then
    cp -f "${NAMESPACE}.trj" "$OUTDIR/xtb_opt.trj"
    TRJ_OK=1
  fi
  if [ "$TRJ_OK" -ne 1 ]; then
    for cand in xtbopt.trj xtb.trj md.trj opt.trj *.trj; do
      [ -f "$cand" ] || continue
      [ -s "$cand" ] || continue
      cp -f "$cand" "$OUTDIR/xtb_opt.trj"
      TRJ_OK=1
      break
    done
  fi
fi



### writing final structure as xtb_opt.mol
if [ "$OPT_OK" -eq 1 ]; then
  SRC_TRJ=""

  # Choose trajectory source (WORKDIR)
  if [ -n "$NAMESPACE" ] && [ -s "${NAMESPACE}.trj" ]; then
    SRC_TRJ="${NAMESPACE}.trj"
  else
    for cand in xtbopt.trj xtb.trj md.trj opt.trj; do
      [ -s "$cand" ] || continue
      SRC_TRJ="$cand"
      break
    done
    if [ -z "$SRC_TRJ" ]; then
      for cand in ./*.trj; do
        [ -s "$cand" ] || continue
        SRC_TRJ="$cand"
        break
      done
    fi
  fi

  # 1) trajectory -> xtb_opt.xyz
  if [ -n "$SRC_TRJ" ] && extract_last_xyz "$SRC_TRJ" "$OUTDIR/xtb_opt.xyz"; then
    :
  else
    # 2) fallback: xtb log -> xtb_opt.xyz
    SRC_LOG=""
    if [ -s "xtbopt.log" ]; then
      SRC_LOG="xtbopt.log"
    elif [ -s "$OUTDIR/xtb_opt.log" ]; then
      SRC_LOG="$OUTDIR/xtb_opt.log"
    fi

    if [ -n "$SRC_LOG" ] && extract_last_xyz "$SRC_LOG" "$OUTDIR/xtb_opt.xyz"; then
      :
    else
      echo "[ERROR] xTB succeeded but no xtb_opt.xyz could be extracted from trajectory/log." >&2
      exit 1
    fi
  fi
else
  # xTB failed: write xtb_failed_last.xyz (best-effort)
  SRC_TRJ=""

  if [ -n "$NAMESPACE" ] && [ -s "${NAMESPACE}.trj" ]; then
    SRC_TRJ="${NAMESPACE}.trj"
  else
    for cand in xtbopt.trj xtb.trj md.trj opt.trj; do
      [ -s "$cand" ] || continue
      SRC_TRJ="$cand"
      break
    done
    if [ -z "$SRC_TRJ" ]; then
      for cand in ./*.trj; do
        [ -s "$cand" ] || continue
        SRC_TRJ="$cand"
        break
      done
    fi
  fi

  if [ -n "$SRC_TRJ" ] && extract_last_xyz "$SRC_TRJ" "$OUTDIR/xtb_failed_last.xyz"; then
    :
  else
    SRC_LOG=""
    if [ -s "xtbopt.log" ]; then
      SRC_LOG="xtbopt.log"
    elif [ -s "$OUTDIR/xtb_opt.log" ]; then
      SRC_LOG="$OUTDIR/xtb_opt.log"
    fi
    if [ -n "$SRC_LOG" ]; then
      extract_last_xyz "$SRC_LOG" "$OUTDIR/xtb_failed_last.xyz" || true
    fi
  fi
fi



### extract optimization trajectory from log if present (XYZ frames)
if [ "$SAVE_TRJ" -eq 1 ] && [ ! -s "$OUTDIR/xtb_opt.trj" ]; then
  SRC_LOG=""
  if [ -f "xtbopt.log" ]; then SRC_LOG="xtbopt.log"; fi
  if [ -z "$SRC_LOG" ]; then SRC_LOG="$OUTDIR/xtb_opt.log"; fi

  if [ -s "$SRC_LOG" ]; then
    awk '
      function isint(s){ return (s ~ /^[[:space:]]*[0-9]+[[:space:]]*$/) }
      {
        if (isint($0)) {
          n = $1
          buf = $0 ORS
          if ((getline line) <= 0) next
          buf = buf line ORS
          ok = 1
          for (i=0; i<n; i++) {
            if ((getline line) <= 0) { ok=0; break }
            buf = buf line ORS
          }
          if (ok) printf "%s", buf
        }
      }
    ' "$SRC_LOG" > "$OUTDIR/xtb_opt.trj" || true

    if [ ! -s "$OUTDIR/xtb_opt.trj" ]; then
      rm -f "$OUTDIR/xtb_opt.trj"
    fi
  fi
fi



exit 0