#!/bin/bash
#SBATCH --job-name=offtarget
#SBATCH --account=compute2-dspencer
#SBATCH --partition=general-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=offtarget.%j.log
#SBATCH --error=offtarget.%j.err
# ---------------------------------------------------------------------------
# Unified CRISPR Off-Target Workflow (-entry OFFTARGET) on RIS Compute2
# (SLURM + Apptainer). ONE wrapper for every cohort: pass the samplesheet and
# an output dir; everything else has a sane default. Replaces the old per-cohort
# wrappers (run_offtarget_{aavs1,aavs1_slurm,cart_slurm}.sh).
#
#   sbatch run_offtarget.sh --input <samplesheet.csv> --outdir <dir> [options]
#
# Options:
#   --input   FILE   samplesheet (required)
#   --outdir  DIR    results dir            (default: ./results_offtarget)
#   --snapshots      render IGV-style tumor|normal pileup PNGs for LIKELY EDITs
#   --no-resume      start fresh instead of -resume
#   --           everything after is passed straight through to `nextflow run`
#
# Launch from a node that can sbatch (login/compute) — NOT the interactive exec
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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT="${2:?--input needs a value}"; shift 2 ;;
        --outdir)    OUTDIR="${2:?--outdir needs a value}"; shift 2 ;;
        --snapshots) SNAP="--offtarget_snapshots true"; shift ;;
        --no-resume) RESUME=""; shift ;;
        --)          shift; EXTRA=("$@"); break ;;
        -h|--help)   sed -n '10,35p' "$0"; exit 0 ;;
        *)           echo "run_offtarget.sh: unknown option '$1' (see --help)" >&2; exit 2 ;;
    esac
done
[[ -n "$INPUT" ]] || { echo "ERROR: --input <samplesheet.csv> is required (see --help)" >&2; exit 2; }

# Lmod: the nextflow module brings its own java; apptainer runs the task containers.
source /etc/profile.d/lmod.sh 2>/dev/null || true
module load nextflow/25.10.4 apptainer/1.4.5

# Cache + binds must match the ris2 apptainer{} block in nextflow.config so the
# image is converted once and the shared filesystems are visible inside tasks.
export NXF_APPTAINER_CACHEDIR="/scratch2/fs1/dspencer/apptainer_cache"
export APPTAINER_CACHEDIR="$NXF_APPTAINER_CACHEDIR"
export APPTAINER_BINDPATH="/storage2,/scratch2"
export NXF_OPTS='-Xmx8g'
mkdir -p "$NXF_APPTAINER_CACHEDIR"

# Under sbatch, $0 is a spool copy — use the submission dir; fall back to the script dir (nohup).
RUNDIR="${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"
cd "$RUNDIR"
export NXF_WORK="${RUNDIR}/work"

nextflow run . -entry OFFTARGET -profile ris2,apptainer \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    -work-dir "$NXF_WORK" \
    $SNAP $RESUME "${EXTRA[@]+"${EXTRA[@]}"}"
