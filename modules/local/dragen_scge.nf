process DRAGEN_SCGE {
    label 'dragen'
    container "${params.dragen_container}"

    input:
    tuple val(normal_meta), path(normal_fastqlist), val(tumor_meta), path(tumor_fastqlist)    

    output:
    tuple val(tumor_meta), path("${tumor_meta.RGSM}*"), emit: dragen_output
    path "versions.yml",    emit: versions

    script:
    """
    /opt/edico/bin/dragen -r /storage1/fs1/dspencer/Active/clinseq/projects/scge/cart_seq/refdata/singh \
      --tumor-fastq-list ${tumor_fastqlist} --fastq-list ${normal_fastqlist} \
      --RGID ${normal_meta.RGID} --RGSM ${normal_meta.RGSM} --RGLB ${normal_meta.RGLB} \
      --RGID-tumor ${tumor_meta.RGID} --RGSM-tumor ${tumor_meta.RGSM} --RGLB-tumor ${tumor_meta.RGLB} \
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
      --output-file-prefix ${tumor_meta.RGID}

      cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    """
    cp -r /storage1/fs1/dspencer/Active/clinseq/projects/scge/workflow/sample_data/Donor_16_CART/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat ${projectDir}/assets/stub/versions/dragen_version.txt)
    END_VERSIONS
    """

}