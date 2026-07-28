#!/usr/bin/env bash
# build_crisprme_index.sh — build the one-time, offline CRISPRme/CRISPRitz genome index that
# the CRISPRME pipeline module consumes read-only (params.crisprme_index_dir).
#
# CRISPRme's per-guide `complete-search` needs a TST index under `genome_library/<PAM>_<bMax>_Genome/`,
# which is specific to the PAM and the *total* bulge budget (bMax = bDNA + bRNA). Building it is the
# slow, heavy step (tens of GB, minutes-to-hours for a full genome), so it is done ONCE here and then
# reused by every guide. This script is NOT part of the Nextflow pipeline — run it by hand on a
# compute node, then point --crisprme_index_dir at its output.
#
# It produces this layout (this whole directory is params.crisprme_index_dir):
#   <outdir>/
#     Genome/                       per-chromosome unzipped FASTAs (CRISPRme requires one file/contig)
#     genome_library/<PAM>_<bMax>_Genome/   the prebuilt TST index (.bin) reused by complete-search
#     <len>bp-<PAM>-<nuclease>.txt   the PAM file (its NAME encodes the nuclease — CRISPRme parses it)
#
# Usage (inside the CRISPRme container):
#   apptainer exec --writable-tmpfs -B /storage2,/scratch2 \
#       docker://pinellolab/crisprme:2.1.10 \
#       bash build_crisprme_index.sh <ref.fa> <outdir> [PAM] [SPACER_LEN] [bMax] [THREADS] [NUCLEASE]
#
# Example (SpCas9, 20 bp spacer, NGG, DNA+RNA bulge budget 2):
#   bash build_crisprme_index.sh hg38_PLVM_CD19_CARv4_cd34.fa /storage2/.../refdata/crisprme NGG 20 2 8 SpCas9
set -euo pipefail

FASTA=${1:?usage: build_crisprme_index.sh <ref.fa> <outdir> [PAM] [SPACER_LEN] [bMax] [THREADS] [NUCLEASE]}
OUTDIR=${2:?missing <outdir>}
PAM=${3:-NGG}
SPLEN=${4:-20}
BMAX=${5:-2}
THREADS=${6:-8}
NUCLEASE=${7:-SpCas9}

# CRISPRme ships in a conda env that its entrypoint activates; Nextflow/apptainer exec bypasses that,
# so make the tools importable here too.
export PATH=/opt/conda/bin:${PATH}

command -v crispritz.py >/dev/null || { echo "ERROR: crispritz.py not on PATH (run inside the CRISPRme container)"; exit 1; }

mkdir -p "${OUTDIR}/Genome"
cd "${OUTDIR}"

# --- 1. PAM file: '<SPLEN Ns><PAM> <pam_len>', named '<len>bp-<PAM>-<nuclease>.txt' -----------------
# The filename is load-bearing: CRISPRme derives the nuclease from basename.split('.')[0].split('-')[2].
PAMFILE="${SPLEN}bp-${PAM}-${NUCLEASE}.txt"
printf '%s%s %d\n' "$(printf 'N%.0s' $(seq 1 "${SPLEN}"))" "${PAM}" "${#PAM}" > "${PAMFILE}"
echo "PAM file: ${PAMFILE}  ->  $(cat "${PAMFILE}")"

# --- 2. split the reference into one unzipped FASTA per contig under Genome/ ------------------------
# CRISPRme requires per-chromosome files. Single awk pass over the (possibly multi-thousand-contig)
# reference; the contig name is the first whitespace token of the header.
if [ -z "$(ls -A Genome 2>/dev/null)" ]; then
    echo "Splitting ${FASTA} into per-contig FASTAs under Genome/ ..."
    awk '/^>/ { name=substr($1,2); close(f); f="Genome/" name ".fa" }
         { print > f }' "${FASTA}"
    echo "  wrote $(ls Genome | wc -l) contig files"
else
    echo "Genome/ already populated ($(ls Genome | wc -l) files) — skipping split"
fi

# --- 3. build the TST index (the slow part) --------------------------------------------------------
# Produces genome_library/<PAM>_<bMax>_Genome/. The genome-dir basename ('Genome') becomes the index
# suffix, so the pipeline module must also pass --genome <dir>/Genome for the reuse lookup to match.
INDEX="genome_library/${PAM}_${BMAX}_Genome"
if [ -d "${INDEX}" ] && [ -n "$(ls -A "${INDEX}" 2>/dev/null)" ]; then
    echo "Index ${INDEX} already exists — skipping build"
else
    echo "Building CRISPRitz index ${INDEX}  (bMax=${BMAX}, ${THREADS} threads) — this is the slow step ..."
    crispritz.py index-genome Genome Genome/ "${PAMFILE}" -bMax "${BMAX}" -th "${THREADS}"
fi

echo
echo "DONE. crisprme_index_dir = ${OUTDIR}"
echo "  Genome/          $(ls Genome | wc -l) contigs"
echo "  ${INDEX}/  $(ls "${INDEX}" 2>/dev/null | wc -l) index files"
echo "  PAM file:        ${PAMFILE}"
