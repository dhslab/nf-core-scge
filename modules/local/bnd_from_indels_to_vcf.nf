process BND_FROM_INDELS_TO_VCF {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(indelfile)

    output:
    tuple val(meta), path("${meta.id}.offtarget_svs.vcf"), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    bnd_from_indels_to_vcf.py \\
        --meta_id ${meta.id} \\
        --indels_path ${indelfile} \\
        --outfile ${meta.id}.offtarget_svs.vcf

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.offtarget_svs.vcf

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}


