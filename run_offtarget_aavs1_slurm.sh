#!/bin/bash
#SBATCH --job-name=offtarget_aavs1
#SBATCH --account=compute2-dspencer
#SBATCH --partition=general-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=offtarget_aavs1.%j.log
#SBATCH --error=offtarget_aavs1.%j.err
# ---------------------------------------------------------------------------
# First real end-to-end run of the Unified CRISPR Off-Target Workflow on
# RIS Compute2 (SLURM + Apptainer). This is the SLURM analogue of the LSF
# wrapper run_offtarget_aavs1.sh; it targets -profile ris2,apptainer.
#
# Paired AAVS1 subset (2 WGS x 6 ECS replicates, guides AAVS1_site5 / _site14).
# The nextflow HEAD process runs here and submits the task jobs to SLURM, so
# this must be launched from a node that can sbatch:
#     sbatch run_offtarget_aavs1_slurm.sh
# NOT from the interactive JupyterLab exec node (it cannot dispatch). If nested
# submission is blocked on your partition, run the body under nohup on a compute
# node instead of via sbatch.
# ---------------------------------------------------------------------------
set -euo pipefail

# Lmod: nextflow module brings its own java; apptainer runs the task containers.
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
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"

nextflow run . -entry OFFTARGET -profile ris2,apptainer \
    --input offtarget_samplesheet_aavs1.csv \
    --outdir ./results_offtarget_aavs1 \
    -resume
