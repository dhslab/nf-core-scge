process GET_INDELS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(files), path(hotspot_file)

    output:
    tuple val(meta), path("${meta.id}.indels.txt")
    path "versions.yml",    emit: versions

    script:
    """
    extract_variant_reads.py ${hotspot_file} ${meta.id}_tumor.cram ${meta.id}.cram > ${meta.id}.indels.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}