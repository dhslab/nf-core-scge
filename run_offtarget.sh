#!/bin/bash
#SBATCH --job-name=offtarget
# Lab defaults (RIS Compute2) so `sbatch run_offtarget.sh ...` needs no extra flags.
# An #SBATCH line cannot read a variable, so these are literal — override them the
# normal SLURM way, which takes precedence: `sbatch --partition=X --account=Y ...`
# (or export SBATCH_PARTITION / SBATCH_ACCOUNT).
#SBATCH --partition=condo-dspencer
#SBATCH --account=compute2-dspencer
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=offtarget.%j.log
#SBATCH --error=offtarget.%j.err
# ---------------------------------------------------------------------------
# Unified CRISPR Off-Target Workflow (-entry OFFTARGET), SLURM + a container
# engine. ONE wrapper for every cohort: pass the samplesheet and an output dir;
# everything else has a sane default and can be overridden.
#
#   sbatch run_offtarget.sh --input <samplesheet.csv> --outdir <dir> [options]
#
# Defaults are the lab's working RIS Compute2 setup, so lab members need no flags:
#     sbatch run_offtarget.sh --input s.csv --outdir out
# Nothing is locked to that cluster, though — every site setting resolves as:
#     command-line flag  >  environment variable  >  inherited from SLURM  >  lab default
# and passing an EMPTY value ("") means "pass nothing, use the cluster's own default".
# On a non-RIS cluster you will typically want --partition/--account/--profile
# (see the last example).
#
# Options (all optional except --input):
#   --input     FILE   samplesheet (required)
#   --outdir    DIR    results dir                    (default: ./results_offtarget)
#   --partition NAME   SLURM partition for the PIPELINE'S OWN jobs   [$OFFTARGET_PARTITION]
#   --account   NAME   SLURM account for those jobs                  [$OFFTARGET_ACCOUNT]
#   --profile   NAME   nextflow profile(s)            (default: ris2,apptainer) [$OFFTARGET_PROFILE]
#   --work-dir  DIR    nextflow work dir              (default: <rundir>/work)  [$OFFTARGET_WORKDIR]
#   --cache     DIR    container image cache          [$OFFTARGET_CACHE]
#   --bind      PATHS  comma-separated container bind paths          [$OFFTARGET_BIND]
#   --modules   LIST   modules to `module load` ("" = none)          [$OFFTARGET_MODULES]
#   --snapshots        render IGV-style tumor|normal pileup PNGs for LIKELY EDITs
#   --no-resume        start fresh instead of -resume
#   --                 everything after is passed straight through to `nextflow run`
#
# PARTITION — read this once, it is the usual stumbling block. There are TWO
# levels of job, set in different places:
#   1. this wrapper job    -> the #SBATCH line above, or `sbatch --partition=NAME ...`
#   2. the pipeline's jobs -> `--partition NAME` here, forwarded to nextflow as
#      --slurm_partition, which the SLURM profile maps to process.queue
# You rarely set (2): it is INHERITED from (1) via $SLURM_JOB_PARTITION. So moving a
# run to another partition takes ONE flag and both levels follow:
#      sbatch --partition=mypart run_offtarget.sh --input s.csv --outdir out
# Account behaves identically (inherits $SLURM_JOB_ACCOUNT).
#
# Examples:
#   # lab default — RIS Compute2, no site flags needed
#   sbatch run_offtarget.sh --input offtarget_samplesheet_aavs1.csv \
#          --outdir results_offtarget_aavs1
#
#   # another partition; the pipeline's jobs follow automatically
#   sbatch --partition=general-cpu run_offtarget.sh --input s.csv --outdir out
#
#   # a completely different SLURM cluster: singularity, its own cache, no RIS binds,
#   # no Lmod, and the site's default account
#   sbatch --partition=batch run_offtarget.sh --input s.csv --outdir out \
#          --profile singularity --account "" \
#          --cache "$HOME/.singularity_cache" --bind "" --modules ""
#
# Launch from a node that can sbatch (login/compute) — NOT an interactive exec
# node, which cannot dispatch. If nested submission is blocked on your partition,
# run the body under nohup on a compute node instead.
#
# Per-process resources are right-sized in conf/modules.config; the ECS
# unevaluable-reads debug log is off by default (offtarget_ecs_unevaluable_log),
# so the whole run stays a few GB and lives comfortably in the repo-local ./work.
#
# For the RIS Compute1 (LSF) path instead, use: nextflow run . -entry OFFTARGET
# -profile ris  (see docs/OFFTARGET_WORKFLOW.md).
# ---------------------------------------------------------------------------
set -euo pipefail

INPUT=""
OUTDIR="./results_offtarget"
SNAP=""
RESUME="-resume"
EXTRA=()

# Defaults are the lab's working RIS Compute2 values, so a lab member runs this with
# no extra flags. Resolution order is
#     --flag  >  $OFFTARGET_*  >  inherited from this job  >  lab default
# Inheriting from $SLURM_JOB_* before the lab default is what makes
# `sbatch --partition=X run_offtarget.sh ...` place the pipeline's own jobs in X too,
# from that one flag. Use ${VAR-default} (no colon) throughout so that setting a
# variable to EMPTY means "pass nothing, use the cluster default" — the non-RIS path.
PARTITION="${OFFTARGET_PARTITION-${SLURM_JOB_PARTITION-condo-dspencer}}"
ACCOUNT="${OFFTARGET_ACCOUNT-${SLURM_JOB_ACCOUNT-compute2-dspencer}}"
PROFILE="${OFFTARGET_PROFILE-ris2,apptainer}"
WORKDIR="${OFFTARGET_WORKDIR-}"
CACHE="${OFFTARGET_CACHE-/scratch2/fs1/dspencer/apptainer_cache}"
BIND="${OFFTARGET_BIND-/storage2,/scratch2}"
MODULES="${OFFTARGET_MODULES-nextflow/25.10.4 apptainer/1.4.5}"
NXF_MEM="${OFFTARGET_NXF_OPTS--Xmx8g}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT="${2:?--input needs a value}"; shift 2 ;;
        --outdir)    OUTDIR="${2:?--outdir needs a value}"; shift 2 ;;
        --partition) PARTITION="${2-}"; shift 2 ;;
        --account)   ACCOUNT="${2-}"; shift 2 ;;
        --profile)   PROFILE="${2:?--profile needs a value}"; shift 2 ;;
        --work-dir)  WORKDIR="${2:?--work-dir needs a value}"; shift 2 ;;
        --cache)     CACHE="${2-}"; shift 2 ;;
        --bind)      BIND="${2-}"; shift 2 ;;
        --modules)   MODULES="${2-}"; shift 2 ;;
        --snapshots) SNAP="--offtarget_snapshots true"; shift ;;
        --no-resume) RESUME=""; shift ;;
        --)          shift; EXTRA=("$@"); break ;;
        # print the banner block between the two '# ---' rules (no brittle line numbers)
        -h|--help)   awk '/^# -{10,}/{n++; next} n==1' "$0"; exit 0 ;;
        *)           echo "run_offtarget.sh: unknown option '$1' (see --help)" >&2; exit 2 ;;
    esac
done
[[ -n "$INPUT" ]] || { echo "ERROR: --input <samplesheet.csv> is required (see --help)" >&2; exit 2; }

# Environment modules are a site convention, not a given: load them only if this
# cluster has Lmod AND the caller has not opted out with --modules "".
if [[ -n "$MODULES" ]]; then
    source /etc/profile.d/lmod.sh 2>/dev/null || true
    if command -v module >/dev/null 2>&1; then
        # shellcheck disable=SC2086
        module load $MODULES || echo "WARN: 'module load $MODULES' failed; " \
            "continuing with whatever is already on PATH" >&2
    else
        echo "WARN: no 'module' command here; skipping --modules '$MODULES'" >&2
    fi
fi
command -v nextflow >/dev/null 2>&1 || {
    echo "ERROR: nextflow is not on PATH. Load it yourself, or pass --modules '<your modules>'." >&2
    exit 127; }

# Container cache + bind paths. Both are site-specific, so an empty value means
# "leave the engine's own default alone" rather than exporting an empty setting.
if [[ -n "$CACHE" ]]; then
    export NXF_APPTAINER_CACHEDIR="$CACHE" APPTAINER_CACHEDIR="$CACHE"
    export NXF_SINGULARITY_CACHEDIR="$CACHE" SINGULARITY_CACHEDIR="$CACHE"
    mkdir -p "$CACHE"
fi
[[ -n "$BIND" ]] && export APPTAINER_BINDPATH="$BIND" SINGULARITY_BIND="$BIND"
[[ -n "$NXF_MEM" ]] && export NXF_OPTS="$NXF_MEM"

# Under sbatch, $0 is a spool copy — use the submission dir; fall back to the script dir (nohup).
RUNDIR="${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"
cd "$RUNDIR"
export NXF_WORK="${WORKDIR:-${RUNDIR}/work}"

# Forward the SLURM placement to the pipeline's own jobs. Only pass what we
# actually know: an unset partition/account must stay unset so the site default
# (or the profile's) applies instead of an empty -p/-A being submitted.
NF_SITE=()
[[ -n "$PARTITION" ]] && NF_SITE+=(--slurm_partition "$PARTITION")
[[ -n "$ACCOUNT"   ]] && NF_SITE+=(--slurm_account "$ACCOUNT")

echo "run_offtarget.sh: profile=${PROFILE} partition=${PARTITION:-<site default>} " \
     "account=${ACCOUNT:-<site default>} work=${NXF_WORK}"

nextflow run . -entry OFFTARGET -profile "$PROFILE" \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    -work-dir "$NXF_WORK" \
    ${NF_SITE[@]+"${NF_SITE[@]}"} \
    $SNAP $RESUME "${EXTRA[@]+"${EXTRA[@]}"}"
