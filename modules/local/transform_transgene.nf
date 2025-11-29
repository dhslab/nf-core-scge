process TRANSFORM_TRANSGENE {
    label 'process_single'
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(meta), path(transgene_file)

    output:
    tuple val(meta), path("${meta.id}.circos_input.tsv"), emit: circos_input
    path "versions.yml"                               , emit: versions

    script:
    """
    transform_transgene.py --input $transgene_file --output ${meta.id}.circos_input.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
} 