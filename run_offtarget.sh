#!/bin/bash
# Unified CRISPR Off-Target Workflow on RIS Compute1 (LSF).
# Leaves the default SCGE pipeline untouched (selected via -entry OFFTARGET).
export LSF_DOCKER_VOLUMES="/storage2/fs1/dspencer/Active:/storage2/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME"
bsub -q dspencer -G compute-dspencer -g /dspencer/nextflow -o offtarget.%J.log -e offtarget.%J.err \
  -a "docker(ghcr.io/dhslab/docker-baseimage:latest)" -app docker1 \
  "export NXF_OPTS='-Xmx8g' && unset SSL_CERT_FILE && cd $PWD && \
   nextflow run . -entry OFFTARGET -profile ris \
     --input offtarget_samplesheet.csv --outdir ./results_offtarget -resume"
