process DRAGEN_SCGE {
    tag "${meta.id}"
    label 'dragen'
    container "${task.ext.dragen_container}"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(reads, stageAs: "fastq_files/*"), path(fastq_list), path(hotspot_vcf)
    tuple val(intermediate_directory_value), path(intermediate_directory)
    path(reference_dir)
    path(adapter1)
    path(adapter2)
    path(sv_noise_file)
    path(snv_noise_file)
    path(nirvana_path)

    output:
    tuple val(meta), path("dragen/*"),                      emit: dragen_output
    path("dragen/${meta.id}_usage.txt"),                    emit: usage, optional: true
    path "versions.yml",                                    emit: versions

    script:
    def exe_path = "${task.ext.dragen_path}"

    def input = [
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.tumor_id} --fastq-list ${fastq_list} --fastq-list-sample-id ${meta.normal_id}" :
            error("Input file is not a CSV file.")
        ].join(' ').trim()

    def alignment_params = [
        "--enable-variant-caller true",
        "--vc-enable-triallelic-filter false",
        "--vc-combine-phased-variants-distance 3",
        "--enable-sv true",
        "--sv-exome true",
        "--sv-output-contigs true",
        "--enable-cnv true",
        "--cnv-use-somatic-vc-baf true",
        "--cnv-enable-self-normalization true",
        "--cnv-enable-ref-calls false",
        "--cnv-somatic-enable-het-calling true",
        "--enable-ploidy-estimator true",
        "--enable-duplicate-marking ${params.mark_duplicates}",
        task.ext.dragen_license_args                  ?: "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        params.alignment_file_format                  ? "--output-format ${params.alignment_file_format}"                     : "",
        adapter1 && adapter2                          ? "--read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2}" : "",
        hotspot_vcf                                   ? "--vc-somatic-hotspots ${hotspot_vcf.min{ it.toString().length() }}"  : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
    ].join(' ').trim()

    """
    mkdir dragen && \\
    ${exe_path}/bin/dragen \\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                ${alignment_params} \\
                ${input} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id}

    # Copy and rename DRAGEN usage
    find dragen/ \\
        -type f \\
        -name "*_usage.txt" \\
        -exec mv "{}" "dragen/${meta.id}_usage.txt" \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${task.ext.dragen_exe_path}/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = "${task.ext.dragen_path}"

    def input = [
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.tumor_id} --fastq-list ${fastq_list} --fastq-list-sample-id ${meta.normal_id}" :
            error("Input file is not a CSV file.")
        ].join(' ').trim()

    def alignment_params = [
        "--enable-variant-caller true",
        "--vc-enable-triallelic-filter false",
        "--vc-combine-phased-variants-distance 3",
        "--enable-sv true",
        "--sv-exome true",
        "--sv-output-contigs true",
        "--enable-cnv true",
        "--cnv-use-somatic-vc-baf true",
        "--cnv-enable-self-normalization true",
        "--cnv-enable-ref-calls false",
        "--cnv-somatic-enable-het-calling true",
        "--enable-ploidy-estimator true",
        "--enable-duplicate-marking ${params.mark_duplicates}",
        task.ext.dragen_license_args                  ?: "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        params.alignment_file_format                  ? "--output-format ${params.alignment_file_format}"                     : "",
        adapter1 && adapter2                          ? "--read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2}" : "",
        hotspot_vcf                                   ? "--vc-somatic-hotspots ${hotspot_vcf.min{ it.toString().length() }}"  : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
    ].join(' ').trim()

    """
    mkdir dragen && \\
    echo ${exe_path}/bin/dragen \\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                ${alignment_params} \\
                ${input} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} > dragen/dragen_command.txt

    touch "dragen/${meta.id}_usage.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echo dragen_version
    END_VERSIONS
    """

}