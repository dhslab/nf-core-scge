process PARSE_INPUT_SAMPLESHEET {
    //tag "${task.ext.prefix.id}"
    label 'process_low'

    container 'docker.io/gregorysprenger/pandas-excel:v2.2.2'

    input:
    path(samplesheet)

    output:
    path("demux_samples.csv")          , optional: true, emit: samples_to_demux
    path("analysis_samples.csv")       , optional: true, emit: samples_to_analyze
    path("alignment_samples.csv")      , optional: true, emit: samples_to_align
    path("versions.yml")               , emit: versions

    script:
    def args = [ 
        samplesheet ? "--input_file ${samplesheet}" : ""
    ].join(' ').trim()

    """
    parse_input_samplesheet.py ${args} --output_dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = [ 
        samplesheet ? "--input_file ${samplesheet}" : ""
    ].join(' ').trim()

    """
    parse_input_samplesheet.py ${args} --output_dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
