process DRAGEN_SCGE {
    label 'dragen'
    label 'dragenalign'
    container "${ext.dragen_aws_image}" ?: "${params.dragen_container}"

    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), val(type), path("*")
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    output:
    tuple val(meta), path("dragen/*"), emit: dragen_output
    path "versions.yml",    emit: versions

    script:
    def intermediate_dir = task.ext.intermediate_dir ? "--intermediate-results-dir ${task.ext.intermediate_dir}" : ""
    def args_license = task.ext.dragen_license_args ?: ''
    def hotspotvcf = dragen_inputs.hotspot_vcf != null ? "--vc-somatic-hotspots inputs/${dragen_inputs.hotspot_vcf}" : ""
    def input = ""
    if (type == 'fastq') {
        input = "--tumor-fastq-list fastq_list.csv --tumor-fastq-list-sample-id ${meta.tumor} --fastq-list fastq_list.csv --fastq-list-sample-id ${meta.normal}"
    } else if (type == 'cram') {
        input = "--tumor-cram-input ${meta.tumor} --cram-input ${meta.normal}"
    }
    if (type == 'bam') {
        input = "--tumor-bam-input ${meta.tumor} --bam-input ${meta.normal}"
    }
    """
    mkdir dragen && \\
    /opt/edico/bin/dragen -r inputs/${dragen_inputs.dragen_hash} ${intermediate_dir} ${input} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps=true \\
                --gc-metrics-enable true \\
                --enable-duplicate-marking ${params.mark_duplicates} \\
                --read-trimmers adapter \\
                --trim-adapter-read1 inputs/${dragen_inputs.dragen_adapter1} \\
                --trim-adapter-read2 inputs/${dragen_inputs.dragen_adapter2} \\
                --enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} \\
                --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} \\
                --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 ${hotspotvcf}\\
                --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 \\
                --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} \\
                --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false \\
                --output-format ${params.alignment_file_format} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id}
                
      cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def intermediate_dir = task.ext.intermediate_dir ? "--intermediate-results-dir ${task.ext.intermediate_dir}" : ""
    def args_license = task.ext.dragen_license_args ?: ''
    def hotspotvcf = dragen_inputs.hotspot_vcf != null ? "--vc-somatic-hotspots inputs/${dragen_inputs.hotspot_vcf}" : ""
    def input = ""
    if (type == 'fastq') {
        input = "--tumor-fastq-list fastq_list.csv --tumor-fastq-list-sample-id ${meta.tumor} --fastq-list fastq_list.csv --fastq-list-sample-id ${meta.normal}"
    } else if (type == 'cram') {
        input = "--tumor-cram-input ${meta.tumor} --cram-input ${meta.normal}"
    }
    if (type == 'bam') {
        input = "--tumor-bam-input ${meta.tumor} --bam-input ${meta.normal}"
    }
    """
    mkdir dragen && \\
    echo /opt/edico/bin/dragen -r inputs/${dragen_inputs.dragen_hash} ${intermediate_dir} ${input} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps=true \\
                --gc-metrics-enable true \\
                --enable-duplicate-marking ${params.mark_duplicates} \\
                --read-trimmers adapter \\
                --trim-adapter-read1 inputs/${dragen_inputs.dragen_adapter1} \\
                --trim-adapter-read2 inputs/${dragen_inputs.dragen_adapter2} \\
                --enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} \\
                --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} \\
                --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 ${hotspotvcf}\\
                --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 \\
                --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} \\
                --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false \\
                --output-format ${params.alignment_file_format} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} > ./dragen/${meta.id}.command.txt
    
    for i in /storage1/fs1/dspencer/Active/clinseq/projects/scge/workflow/sample_data/Donor_16_CART/*;
    do 
        touch ./dragen/\$(basename \$i);
    done    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat ${projectDir}/assets/stub/versions/dragen_version.txt)
    END_VERSIONS
    """

}