process GET_INDELS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.results.txt")
    path "versions.yml",    emit: versions

    script:
    """
    getIndelsFromBam.py "$projectDir/assets/accessory_files/trac_trbc_merged.idt_offtarget.txt" ${meta.id}_tumor.cram ${meta.id}.cram > ${meta.id}.results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}