process DRAGEN_SCGE {
    label 'dragen'
    label 'dragenalign'
    container "${ ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? params.aws_dragen_container : params.dragen_container }"    
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(tumor_reads, stageAs: "tumor_fastq_files/*"), path(tumor_fastq_list), path(normal_reads, stageAs: "normal_fastq_files/*"), path(normal_fastq_list)
    tuple val(intermediate_directory_value), path(intermediate_directory)
    path(reference_dir)
    path(adapter1)
    path(adapter2)
    path(cram_reference)
    path(sv_noise_file)
    path(snv_noise_file)
    path(cnv_population_vcf)
    path(hotspots)

    output:
    tuple val(meta), path("dragen/*"),   emit: dragen_output
    path("dragen/${meta.id}_usage.txt"), emit: usage, optional: true
    path "versions.yml",                 emit: versions

    script:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "/opt/edico" : "/opt/dragen/4.3.6"

    def input_args = []

    if (tumor_fastq_list?.toString()?.endsWith('.csv') && normal_fastq_list?.toString()?.endsWith('.csv')) {
        input_args << "--tumor-fastq-list ${tumor_fastq_list}"
        input_args << "--tumor-fastq-list-sample-id ${meta.edited_id}"
        input_args << "--fastq-list ${normal_fastq_list}"
        input_args << "--fastq-list-sample-id ${meta.control_id}"
    }
    if (input_args.isEmpty()) {
        error("No valid input provided. Expected a fastq_list.csv file, or one or more BAM/CRAM files.")
    }

    def input = input_args.join(' ').trim()

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
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"         : ""
//        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                : ""
    ].join(' ').trim()

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

    # Copy and rename DRAGEN usage
    find dragen/ \\
        -type f \\
        -name "*_usage.txt" \\
        -exec mv "{}" "dragen/${meta.id}_usage.txt" \\;

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        dragen: \${DRAGEN_VERSION}
    END_VERSIONS
    """

    stub:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "/opt/edico" : "/opt/dragen/4.3.6"

    def input = [
        edited_alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--tumor-bam-input ${it}"  }                                        ?:
        control_alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--bam-input ${it}"       }                                        ?:
        edited_alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--tumor-cram-input ${it}" }                                        ?:
        control_alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--cram-input ${it}"      }                                        ?:
        fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list fastq_list.csv --tumor-fastq-list-sample-id ${meta.edited_id} --fastq-list fastq_list.csv --fastq-list-sample-id ${meta.control_id}" :
        error("Input file is not a BAM, CRAM, or CSV file.")
    ].join(' ').trim()

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
    ${task.process}:
        dragen: \${DRAGEN_VERSION}
    END_VERSIONS
    """

}