process PARSE_INPUT_SAMPLESHEET {
//    tag "${task.ext.prefix.id}"
    label 'process_low'

    container 'docker.io/gregorysprenger/pandas-excel:v2.2.2'

    input:
    path(samplesheet)

    output:
    path("alignment_samples.csv")       , optional: true, emit: samples_to_align
    path("analysis_samples.csv")        , optional: true, emit: samples_to_analyze
    path("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_input_samplesheet.py \\
        --input_file ${samplesheet} \\
        --output_dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    parse_input_samplesheet.py \\
        --input_file ${samplesheet} \\
        --output_dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
