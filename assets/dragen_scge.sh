#!/bin/bash

#bsub -g /dspencer/adhoc -G compute-dspencer -q dragen-4 -M 350G -n 30 -R "span[hosts=1] select[mem>350G] rusage[mem=350G]" -a 'docker(seqfu/oracle8-dragen-4.0.3:latest)'

/opt/edico/bin/dragen -r /storage1/fs1/dspencer/Active/clinseq/projects/scge/cart_seq/refdata/singh \
      -1 /storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/research/CLE-5280/demux_fastq_2/Donor_16_Unedited-DNA_S1_L006_R1_001.fastq.gz \
      -2 /storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/research/CLE-5280/demux_fastq_2/Donor_16_Unedited-DNA_S1_L006_R2_001.fastq.gz \
      --tumor-fastq1 /storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/research/CLE-5280/demux_fastq_1/Donor_16_Edited-DNA_S2_R1_001.fastq.gz \
      --tumor-fastq2 /storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/research/CLE-5280/demux_fastq_1/Donor_16_Edited-DNA_S2_R2_001.fastq.gz \
      --RGID Donor_16_Unedited-DNA.ATAGTCTAGC.ATGTCGTATT --RGSM Donor_16_Unedited-DNA --RGLB Donor_16_Unedited-DNA-lib1 \
      --RGID-tumor Donor_16_Edited-DNA.CCGCATATTC.TCAGTCTCGT --RGSM-tumor Donor_16_Edited-DNA --RGLB-tumor Donor_16_Edited-DNA-lib1 \
      --read-trimmers adapter \
      --trim-adapter-read1 /storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/t2t-chm13_adapter1.fa \
      --trim-adapter-read2 /storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/t2t-chm13_adapter2.fa \
      --enable-map-align true \
      --enable-map-align-output true \
      --enable-bam-indexing true \
      --enable-duplicate-marking true \
      --qc-coverage-ignore-overlaps true \
      --gc-metrics-enable true \
      --enable-variant-caller true \
      --vc-combine-phased-variants-distance 3 \
      --dbsnp /storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/dbsnp.vcf.gz \
      --enable-sv true \
      --sv-output-contigs true \
      --sv-hyper-sensitivity true \
      --sv-use-overlap-pair-evidence true \
      --enable-cnv true \
      --cnv-use-somatic-vc-baf true \
      --cnv-somatic-enable-het-calling true \
      --cnv-enable-ref-calls false \
      --output-format CRAM \
      --intermediate-results-dir /staging/ \
      --output-directory /storage1/fs1/dspencer/Active/clinseq/projects/scge/cart_seq/Donor_16_CART \
      --output-file-prefix Donor_16_CART
