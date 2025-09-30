process DRAGEN_SCGE {
    label 'dragen'
    label 'dragenalign'
    container "${ ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? params.aws_dragen_container : params.dragen_container }"    
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), val(type), val(crams)
    path(hotspot_file)
    tuple val(intermediate_directory_value), path(intermediate_directory)
    path(reference_dir)
    path(adapter1)
    path(adapter2)
    path(cram_reference)
    path(sv_noise_file)
    path(snv_noise_file)
    path(cnv_population_vcf)

    output:
    tuple val(meta), path("dragen/*"),   emit: dragen_output
    path("dragen/${meta.id}_usage.txt"), emit: usage, optional: true
    path "versions.yml",                 emit: versions

    script:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "/opt/edico" : "/opt/dragen/4.3.6"

    def alignment_params = [
        task.ext.dragen_license_args                  ?: "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                         : "",
        adapter1                                      ? "--trim-adapter-read1 ${adapter1}"                                   : "",
        adapter2                                      ? "--trim-adapter-read2 ${adapter2}"                                   : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"    : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"   : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                             : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                            : "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"               : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"         : "",
        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                : ""
    ].join(' ').trim()

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
    mkdir -p dragen

    ${exe_path}/bin/dragen \\
        ${input} \\
        ${alignment_params} \\
        --force \\
        --enable-sv true \\
        --enable-cnv true \\
        --enable-sort true \\
        --output-format CRAM \\
        --read-trimmers adapter \\
        --enable-map-align true \\
        --gc-metrics-enable true \\
        --sv-output-contigs true \\
        --enable-bam-indexing true \\
        --cnv-enable-ref-calls false \\
        --enable-variant-caller true \\
        --cnv-use-somatic-vc-baf true \\
        --enable-map-align-output true \\
        --enable-duplicate-marking true \\
        --qc-coverage-ignore-overlaps true \\
        --vc-enable-triallelic-filter false \\
        --sv-use-overlap-pair-evidence true \\
        --cnv-somatic-enable-het-calling true \\
        --vc-combine-phased-variants-distance 3 \\
        --output-directory ./dragen \\
        --output-file-prefix ${meta.id}
                
      cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "/opt/edico" : "/opt/dragen/4.3.6"

    def alignment_params = [
        task.ext.dragen_license_args                  ?: "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                         : "",
        adapter1                                      ? "--trim-adapter-read1 ${adapter1}"                                   : "",
        adapter2                                      ? "--trim-adapter-read2 ${adapter2}"                                   : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"    : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"   : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                             : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                            : "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"               : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"         : "",
        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                : ""
    ].join(' ').trim()

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
    mkdir -p dragen

    echo ${exe_path}/bin/dragen \\
            ${input} \\
            ${alignment_params} \\
            --force \\
            --enable-sv true \\
            --enable-cnv true \\
            --enable-sort true \\
            --output-format CRAM \\
            --read-trimmers adapter \\
            --enable-map-align true \\
            --gc-metrics-enable true \\
            --sv-output-contigs true \\
            --enable-bam-indexing true \\
            --cnv-enable-ref-calls false \\
            --enable-variant-caller true \\
            --cnv-use-somatic-vc-baf true \\
            --enable-map-align-output true \\
            --enable-duplicate-marking true \\
            --qc-coverage-ignore-overlaps true \\
            --vc-enable-triallelic-filter false \\
            --sv-use-overlap-pair-evidence true \\
            --cnv-somatic-enable-het-calling true \\
            --vc-combine-phased-variants-distance 3 \\
            --output-directory ./dragen \\
            --output-file-prefix ${meta.id}
                
     > ./dragen/${meta.id}.command.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat ${projectDir}/assets/stub/versions/dragen_version.txt)
    END_VERSIONS
    """

}