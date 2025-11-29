process BND_FROM_INDELS_TO_VCF {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(indels_txt)

    output:
    tuple val(meta), path("${meta.id}.indels.bnd.vcf"), emit: bnd_vcf
    path "versions.yml",    emit: versions

    script:
    """
    bnd_from_indels_to_vcf.py \\
        --meta_id ${meta.id} \\
        --indels_path ${indels_txt}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}


