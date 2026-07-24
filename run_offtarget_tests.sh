#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_offtarget_tests.sh — one-shot test harness for the Unified CRISPR
# Off-Target Workflow (-entry OFFTARGET).
#
# Runs the cheap-to-medium tiers that DON'T need a real cohort / scheduler:
#   Tier 0  static parse    python ast on bin/*.py  +  nextflow -preview
#   Tier 1  glue unit tests  pytest tests/test_offtarget_glue.py (the ECS<->WGS join)
#   Tier 2  model contract   wgs_shape_model.pkl loads + sklearn version guard
#   Tier 3  stub end-to-end  nextflow -profile stub -stub-run (DAG wiring, no CRAMs)
#
# Tiers 1 & 2 run INSIDE the pipeline's own container (Apptainer), so they use
# the exact pinned deps the pipeline runs — python 3.11, scikit-learn 1.8.0
# (the version wgs_shape_model.pkl was pickled under), pandas, pysam. That makes
# Tier 2 a true mirror of production and sidesteps the repo-root vendored
# numpy/pandas dirs that otherwise shadow a bare interpreter.
#
# Best run on a normal compute node (module + apptainer available). It auto-loads
# the `apptainer` and `nextflow` Lmod modules; you can also pre-load them yourself:
#     module load apptainer nextflow
# Anything genuinely missing (no apptainer, no java) is reported SKIP, not FAIL.
# The real AAVS1 acceptance run needs a cohort and stays separate:
#     sbatch run_offtarget.sh --input offtarget_samplesheet_aavs1.csv --outdir results_offtarget_aavs1
#
# Usage:   bash run_offtarget_tests.sh
#   OFFTARGET_CONTAINER=ghcr.io/dhslab/docker-scge-offtarget:TAG  bash run_offtarget_tests.sh
#   OFFTARGET_TEST_CACHE=/path/to/cache  bash run_offtarget_tests.sh   # where the SIF is cached
#   INSTALL_DEPS=1 bash run_offtarget_tests.sh                         # allow the venv fallback
# Exit code is nonzero only if a tier FAILED (SKIP never fails the run).
# ---------------------------------------------------------------------------
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"
SS="assets/offtarget_samplesheet_template.csv"
SYS_PY="${PYTHON:-python3}"          # dependency-free; used only for the ast parse tier

# Pin the container straight from the pipeline modules so the harness can't drift
# from what actually runs. Override with OFFTARGET_CONTAINER.
CONTAINER="$(grep -rhoE 'ghcr\.io/dhslab/docker-scge-offtarget:[0-9A-Za-z._-]+' modules/local/*.nf 2>/dev/null | sort -u | head -1)"
CONTAINER="${OFFTARGET_CONTAINER:-${CONTAINER:-ghcr.io/dhslab/docker-scge-offtarget:260710}}"
CACHE="${OFFTARGET_TEST_CACHE:-$REPO/.offtarget_testcache}"
SIF="$CACHE/$(printf '%s' "$CONTAINER" | tr '/:' '__').sif"
PYX="$CACHE/pyextra"                 # pure-python pytest overlay (pytest isn't in the prod image)

pass=0; fail=0; skip=0; FAILED_NAMES=()
c_g=$'\033[32m'; c_r=$'\033[31m'; c_y=$'\033[33m'; c_b=$'\033[1m'; c_0=$'\033[0m'
hdr(){ printf '\n%s========== %s ==========%s\n' "$c_b" "$1" "$c_0"; }
ok(){   printf '  %sPASS%s %s\n' "$c_g" "$c_0" "$1"; pass=$((pass+1)); }
no(){   printf '  %sFAIL%s %s\n' "$c_r" "$c_0" "$1"; fail=$((fail+1)); FAILED_NAMES+=("$1"); }
sk(){   printf '  %sSKIP%s %s (%s)\n' "$c_y" "$c_0" "$1" "$2"; skip=$((skip+1)); }
# a nextflow run failing on the JVM (missing / <17) is an env issue, not a pipeline defect -> SKIP
nf_env_broke(){ grep -qiE 'Cannot find Java|wrong version|java: command not found|needs a Java|Java .* (or later|is required)' "$1"; }

# --- Lmod plumbing ---------------------------------------------------------
# `module` is a shell function that a non-login `bash script.sh` doesn't inherit;
# re-source the init if it's gone. MODULEPATH is an env var and *is* inherited, so
# the lab's own modulefiles stay visible either way.
init_modules(){
  command -v module >/dev/null 2>&1 && return 0
  local i
  for i in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh \
           "${LMOD_PKG:-/usr/share/lmod/lmod}/init/bash" /usr/share/lmod/lmod/init/bash; do
    [ -f "$i" ] && . "$i" 2>/dev/null && command -v module >/dev/null 2>&1 && return 0
  done
  command -v module >/dev/null 2>&1
}
newest_module(){ module -t avail 2>&1 | grep -iE "^$1/[0-9]" | sort -V | tail -1; }
# load_tool <exe> <module-basename>: make <exe> resolvable, loading the newest module if needed
load_tool(){
  command -v "$1" >/dev/null 2>&1 && return 0
  init_modules || return 1
  local m; m="$(newest_module "$2")" || return 1
  [ -n "$m" ] && module load "$m" >/dev/null 2>&1
  command -v "$1" >/dev/null 2>&1
}

# --- python backend: container (preferred) or venv (fallback) --------------
MODE_PY=none; CTR=""; VENVPY=""; SKLEARN_OK=0
# run python with the pipeline deps, always from a neutral cwd (/): the repo root
# has broken vendored numpy/pandas/... dirs that would otherwise shadow real packages
# (same reason the glue-test conftest invokes bin scripts with cwd=tmp).
run_py(){
  case "$MODE_PY" in
    container) ( cd / && $CTR env ${PYX:+PYTHONPATH="$PYX"} REPO="$REPO" python "$@" );;
    venv)      ( cd / && REPO="$REPO" "$VENVPY" "$@" );;
    *) return 127;;
  esac
}

setup_container(){
  local rt=""
  load_tool apptainer apptainer && rt=apptainer
  [ -z "$rt" ] && load_tool singularity singularity && rt=singularity
  [ -z "$rt" ] && return 1
  mkdir -p "$CACHE"
  if [ ! -s "$SIF" ]; then
    echo "  pulling $CONTAINER (one-time, ~300MB)..."
    APPTAINER_CACHEDIR="$CACHE/.pull" SINGULARITY_CACHEDIR="$CACHE/.pull" \
      "$rt" pull --force "$SIF" "docker://$CONTAINER" >"$CACHE/pull.log" 2>&1 \
      || { echo "  container pull failed (see $CACHE/pull.log)"; return 1; }
  fi
  CTR="$rt exec --bind $REPO:$REPO --bind $CACHE:$CACHE $SIF"
  MODE_PY=container
  SKLEARN_OK=1   # image carries scikit-learn at the pickling version -> Tier 2 is a true mirror
  # pytest is test-only and not in the prod image; drop a pure-python overlay beside the SIF
  if ( cd / && $CTR python -c "import pytest" ) >/dev/null 2>&1; then
    PYX=""
  elif [ ! -d "$PYX/_pytest" ]; then
    $CTR pip install --quiet --target="$PYX" pytest >"$CACHE/pytest_install.log" 2>&1 \
      || echo "  note: could not add pytest overlay (Tier 1 will SKIP; see $CACHE/pytest_install.log)"
  fi
  return 0
}

# Fallback for hosts with no container runtime. Opt-in (INSTALL_DEPS=1) because it
# builds a venv. A venv on shared storage hard-references a per-node interpreter, so
# we rebuild (--clear) whenever the recorded python no longer executes here.
setup_venv(){
  [ "${INSTALL_DEPS:-0}" = "1" ] || return 1
  mkdir -p "$CACHE"
  local py="${PYTHON:-python3}" pyver
  pyver="$("$py" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null || echo x)"
  local venv="$REPO/.offtarget_testenv-$pyver"
  VENVPY="$venv/bin/python"
  if [ ! -x "$VENVPY" ] || ! ( cd / && "$VENVPY" -c "import pandas,pytest,joblib,matplotlib" ) 2>/dev/null; then
    "$py" -m venv --clear "$venv" >"$CACHE/venv.log" 2>&1 || return 1
    "$VENVPY" -m pip install --quiet --upgrade pip >>"$CACHE/venv.log" 2>&1
    "$VENVPY" -m pip install --quiet pytest pandas joblib matplotlib >>"$CACHE/venv.log" 2>&1 || return 1
    "$VENVPY" -m pip install --quiet "scikit-learn==1.8.0" >>"$CACHE/venv.log" 2>&1 && SKLEARN_OK=1
  else
    ( cd / && "$VENVPY" -c "import sklearn" ) 2>/dev/null && SKLEARN_OK=1
  fi
  MODE_PY=venv
  return 0
}

NF=""
setup_nextflow(){
  if load_tool nextflow nextflow; then NF="$(command -v nextflow)"; return 0; fi
  if [ -x "$REPO/nextflow" ] && command -v java >/dev/null 2>&1; then NF="$REPO/nextflow"; return 0; fi
  return 1
}

printf '%sOFF-TARGET PIPELINE TESTS%s   repo=%s\n' "$c_b" "$c_0" "$REPO"
printf 'container=%s\n' "$CONTAINER"

# --- environment setup -----------------------------------------------------
hdr "environment setup"
init_modules >/dev/null 2>&1 || true
if   setup_container; then echo "  python deps: container ($CONTAINER)"
elif setup_venv;      then echo "  python deps: venv ($VENVPY)"
else echo "  python deps: unavailable (load the apptainer module, or re-run with INSTALL_DEPS=1)"; fi
if setup_nextflow; then echo "  nextflow: $NF  |  $(java -version 2>&1 | head -1)"
else echo "  nextflow/java: unavailable (nextflow tiers will SKIP)"; fi

# ===========================================================================
hdr "Tier 0 — static parse (fast, no container)"

# 0a. python syntax of every bin/ script (ast only, no third-party deps)
BAD=""
for s in "$REPO"/bin/*.py; do
  "$SYS_PY" -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$s" 2>/dev/null \
    || BAD="$BAD $(basename "$s")"
done
[ -z "$BAD" ] && ok "all bin/*.py parse" || no "python syntax errors in:$BAD"

# 0b. nextflow can compile the OFFTARGET entry + resolve the ris profile
if [ -z "$NF" ]; then
  sk "nextflow -preview" "no nextflow / java"
else
  if "$NF" run . -entry OFFTARGET -profile ris --input "$SS" -preview >/tmp/off_preview.log 2>&1; then
    ok "nextflow compiles OFFTARGET entry (-preview, ris profile)"
  elif nf_env_broke /tmp/off_preview.log; then
    sk "nextflow -preview" "Java 17+ not available in this env"
  else
    no "nextflow -preview (see /tmp/off_preview.log)"; tail -n 15 /tmp/off_preview.log | sed 's/^/      /'
  fi
fi

# ===========================================================================
hdr "Tier 1 — glue-logic unit tests (ECS<->WGS join + off-target combiner)"
if [ "$MODE_PY" = none ]; then
  sk "pytest tests/test_offtarget_glue.py tests/test_combine_offtarget.py" "no python backend (load apptainer, or INSTALL_DEPS=1)"
elif ! run_py -c "import pytest, pandas" >/dev/null 2>&1; then
  sk "pytest tests/test_offtarget_glue.py tests/test_combine_offtarget.py" "pytest/pandas unavailable"
else
  if run_py -m pytest -q "$REPO/tests/test_offtarget_glue.py" "$REPO/tests/test_combine_offtarget.py" >/tmp/off_pytest.log 2>&1; then
    ok "$(grep -Eo '[0-9]+ passed' /tmp/off_pytest.log | tail -1) — glue + combiner tests"
  else
    no "pytest glue/combiner tests (see /tmp/off_pytest.log)"; tail -n 20 /tmp/off_pytest.log | sed 's/^/      /'
  fi
fi

# ===========================================================================
hdr "Tier 2 — model contract (sklearn pin / version guard)"
if [ "$MODE_PY" = none ]; then
  sk "model load + version guard" "no python backend"
elif [ "$SKLEARN_OK" != 1 ] || ! run_py -c "import joblib, sklearn, pandas" >/dev/null 2>&1; then
  sk "model load + version guard" "scikit-learn 1.8.0 unavailable (use the container)"
else
  run_py - >/tmp/off_model.log 2>&1 <<'PY'
import os, sys
repo = os.environ["REPO"]
sys.path.insert(0, os.path.join(repo, "bin"))   # for features.py (absolute; cwd is neutral)
import joblib
from features import check_sklearn_version
b = joblib.load(os.path.join(repo, "assets/models/wgs_shape_model.pkl"))
model = b["model"] if isinstance(b, dict) else b
check_sklearn_version(model, name="wgs_shape_model.pkl")   # raises if runtime older than pickle
print("features:", b.get("features") if isinstance(b, dict) else None)
print("OK")
PY
  if grep -q '^OK' /tmp/off_model.log; then
    ok "wgs_shape_model.pkl loads + version guard passes"
    grep -i 'features:' /tmp/off_model.log | sed 's/^/      /'
  else
    no "model load / version guard (see /tmp/off_model.log)"; tail -n 15 /tmp/off_model.log | sed 's/^/      /'
  fi
fi

# ===========================================================================
hdr "Tier 3 — stub end-to-end (DAG wiring, no CRAMs)"
if [ -z "$NF" ]; then
  sk "nextflow -stub-run" "no nextflow / java"
else
  rm -rf "$REPO/results_stub"
  if "$NF" run . -entry OFFTARGET -profile stub -stub-run \
        --input "$SS" --outdir "$REPO/results_stub" >/tmp/off_stub.log 2>&1; then
    ok "stub run completed (both arms + paired join scheduled)"
  elif nf_env_broke /tmp/off_stub.log; then
    sk "nextflow -stub-run" "Java 17+ not available in this env"
  else
    no "stub run (see /tmp/off_stub.log)"; tail -n 20 /tmp/off_stub.log | sed 's/^/      /'
  fi

  # HOTSPOTS entry (gRNA -> targets) — pure-shell stubs, tiny committed fixture FASTA.
  rm -rf "$REPO/results_hotspots_stub"
  if "$NF" run . -entry HOTSPOTS -profile stub -stub-run \
        --input assets/grna_samplesheet_template.csv --fasta assets/stub/tiny.fa \
        --outdir "$REPO/results_hotspots_stub" >/tmp/hotspots_stub.log 2>&1; then
    ok "HOTSPOTS stub run completed (Cas-OFFinder -> combine -> targets.vcf)"
  elif nf_env_broke /tmp/hotspots_stub.log; then
    sk "HOTSPOTS -stub-run" "Java 17+ not available in this env"
  else
    no "HOTSPOTS stub run (see /tmp/hotspots_stub.log)"; tail -n 20 /tmp/hotspots_stub.log | sed 's/^/      /'
  fi
fi

# ===========================================================================
hdr "SUMMARY"
printf '  %sPASS %d%s   %sFAIL %d%s   %sSKIP %d%s\n' \
  "$c_g" "$pass" "$c_0" "$c_r" "$fail" "$c_0" "$c_y" "$skip" "$c_0"
if [ "$fail" -gt 0 ]; then printf '  failed: %s\n' "${FAILED_NAMES[*]}"; fi
cat <<EOF

  Not covered here (needs a real cohort — run from a compute node):
    Tier 4  real AAVS1 acceptance run:  sbatch run_offtarget.sh --input offtarget_samplesheet_aavs1.csv --outdir results_offtarget_aavs1
            then check results_offtarget_aavs1/offtarget/ for:
              training.tsv         non-empty, labels {0,1}   (empty => join broke)
              offtarget_report.csv AAVS1 on-target: is_hotspot=1 & ecs_confirmed=1
              recall_vs_vaf.{csv,png}  ~1.0 in high-VAF bins
  Known blind spots (no data yet): sub-5% VAF floor; de-novo off-target positive control.
EOF
[ "$fail" -eq 0 ]
