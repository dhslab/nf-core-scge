process MAKE_FASTQLIST {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(meta), val(fastqlist)

    output:
    tuple val(meta), path('fastq_list.csv'), emit: fastq_list

    script:
    """
    echo -e "${fastqlist}" > fastq_list.csv

    """

    stub:
    """
    echo -e "${fastqlist}" > fastq_list.csv

    """
}
