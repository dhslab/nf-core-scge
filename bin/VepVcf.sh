#!/bin/bash

# bash script to submit a job that annotates a vcf file with VEP

reference="/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/sequence/hg38_mgi_patch.fa"
cytobands="/storage1/fs1/dspencer/Active/shared/refdata/hg38/hg38.cytoBandIdeo.bed.gz"
vepcache="/storage1/fs1/dspencer/Active/shared/refdata/hg38/VEP_cache"

# input arguments
vcf=$1
out=$2

jobgroup="/dspencer/adhoc"
computegroup="compute-dspencer"
queue="general"

bsub -g ${jobgroup} -G ${computegroup} -q ${queue} -oo %J.vep.log -eo %J.vep.err \
     -M 8000000 -R"select[mem>8000] rusage[mem=8000]" -q general -a "docker(ghcr.io/dhslab/docker-vep_release113:250810)" \
     /usr/bin/perl \
        -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \
        --format vcf \
        --vcf --fasta ${reference} \
        --hgvs \
        --symbol \
        --term SO \
        --flag_pick \
        --custom ${cytobands},cytobands,bed \
        -o ${out} \
        -i ${vcf} \
        --offline \
        --cache \
        --max_af --dir ${vepcache}
