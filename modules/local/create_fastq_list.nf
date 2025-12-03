process CREATE_FASTQ_LIST {
    tag "${meta.id}"

    container "ghcr.io/dhslab/docker-python3:240604"

    input:
    tuple val(meta), path(read1file), path(read2file), path(runinfo)

    output:
    tuple val(meta), path("*fastq_list.csv"), emit: fastq_list
    path("versions.yml")   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastq_list_args = [
        read1file  ? "--read1 ${read1file}"      : "",
        read2file  ? "--read2 ${read2file}"      : "",
        runinfo    ? "--runinfo ${runinfo}"      : ""
    ].join(' ').trim()
    """
    create_fastq_list.py \\
        -i ${meta.id} \\
        ${fastq_list_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        \$(create_fastq_list.py -v)
    END_VERSIONS
    """

    stub:
    def fastq_list_args = [
        read1file  ? "--read1 ${read1file}"      : "",
        read2file  ? "--read2 ${read2file}"      : "",
        runinfo    ? "--runinfo ${runinfo}"      : "",
    ].join(' ').trim()
    """
    create_fastq_list.py \\
        -i ${meta.id} \\
        ${fastq_list_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        \$(create_fastq_list.py -v)
    END_VERSIONS
    """

}
