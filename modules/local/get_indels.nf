process GET_INDELS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(meta), path(files)
    path(targetfile)

    output:
    tuple val(meta), path("${meta.id}.indels.txt")
    path "versions.yml",    emit: versions

    script:
    """
    getIndelsFromBam.py ${targetfile} ${meta.id}_tumor.cram ${meta.id}.cram > ${meta.id}.indels.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}