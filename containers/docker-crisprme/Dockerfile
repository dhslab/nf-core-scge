# ghcr.io/dhslab/docker-crisprme
#
# Thin wrapper over the upstream CRISPRme image that adds `procps` (the `ps` command). Nextflow's
# per-task wrapper runs `ps` to collect resource metrics under `-euo pipefail`; the upstream
# micromamba-based pinellolab/crisprme image has no `ps`, so tasks abort before the payload runs
# ("Command 'ps' required by nextflow to collect task metrics cannot be found"). Everything else
# (crisprme.py, crispritz, samtools, bedtools in the base conda env) is inherited unchanged.
#
# To publish: copy this directory into the lab image repo (dhslab/dhslab-docker-images) as
# `docker-crisprme/` and push to main; push_docker.yaml auto-builds it to
# ghcr.io/dhslab/docker-crisprme with :latest + :YYMMDD tags. Then modules/local/crisprme.nf's
# `container` (ghcr.io/dhslab/docker-crisprme:latest) resolves.
FROM pinellolab/crisprme:latest

USER root
# Install procps via apt so `ps` lands in /usr/bin (on the container's default PATH). A
# micromamba/conda install would put it under /opt/conda/bin, which the image does NOT have on
# its default PATH, so Nextflow's wrapper (which bypasses conda activation) would not find it.
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
