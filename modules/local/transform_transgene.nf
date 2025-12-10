process TRANSFORM_TRANSGENE {
    tag "${meta.id}"
    label 'process_single'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(meta), path(transgene_file)

    output:
    tuple val(meta), path("${meta.id}.circos_input.tsv"), emit: circos_input
    path "versions.yml"                                 , emit: versions

    script:
    """
    export PATH=\$PATH:/usr/local/bin
    transform_transgene.py --input $transgene_file --output ${meta.id}.circos_input.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
} 