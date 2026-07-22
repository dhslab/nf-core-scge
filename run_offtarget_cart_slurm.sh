#!/bin/bash
#SBATCH --job-name=offtarget_cart
#SBATCH --account=compute2-dspencer
#SBATCH --partition=general-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=offtarget_cart.%j.log
#SBATCH --error=offtarget_cart.%j.err
# ---------------------------------------------------------------------------
# CART cohort through the Unified CRISPR Off-Target Workflow on RIS Compute2
# (SLURM + Apptainer). Paired ECS+WGS, 25 guides / 65 samples, incl. the two
# confirmed off-targets (PLCB2 chr12:32,679,410; CNNM3 chr10:102,919,984).
#
# --offtarget_snapshots true renders IGV-style tumor|normal pileup PNGs for every
# LIKELY EDIT (published to results_offtarget_cart/offtarget/), so the confirmed
# off-targets can be eyeballed for verification.
#
# Launch from a node that can sbatch (login/compute), NOT the interactive exec node:
#     sbatch run_offtarget_cart_slurm.sh
# ---------------------------------------------------------------------------
set -euo pipefail

source /etc/profile.d/lmod.sh 2>/dev/null || true
module load nextflow/25.10.4 apptainer/1.4.5

export NXF_APPTAINER_CACHEDIR="/scratch2/fs1/dspencer/apptainer_cache"
export APPTAINER_CACHEDIR="$NXF_APPTAINER_CACHEDIR"
export APPTAINER_BINDPATH="/storage2,/scratch2"
export NXF_OPTS='-Xmx8g'
mkdir -p "$NXF_APPTAINER_CACHEDIR"

# Work dir on /storage2 (repo-local ./work). With offtarget_ecs_unevaluable_log=false
# (the default), ECS_INDELS writes only its small .tsv/.vcf (~1 MB/sample) instead of the
# ~0.5-1 TB/sample unevaluable-reads debug log that previously filled storage2 (ENOSPC)
# AND the shared /scratch2 group quota (EDQUOT). The whole run is now a few GB, so it fits
# comfortably here and avoids the scratch2 group-quota dependency entirely.
export NXF_WORK="${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}/work"

cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"

nextflow run . -entry OFFTARGET -profile ris2,apptainer \
    --input /storage2/fs1/dspencer/Active/clinseq/projects/scge/cart_seq/offtarget_samplesheet_cart.csv \
    --outdir ./results_offtarget_cart \
    --offtarget_snapshots true \
    -work-dir "$NXF_WORK" \
    -resume
