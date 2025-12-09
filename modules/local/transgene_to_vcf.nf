process TRANSGENE_TO_VCF {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(transgene_file)

    output:
    tuple val(meta), path("${meta.id}.transgene.vcf"), emit: transgene_vcf
    path "versions.yml", emit: versions

    script:
    """
    transgene2vcf.py $transgene_file ${meta.id}.transgene.vcf

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        transgene2vcf: \$(transgene2vcf.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.transgene.vcf

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        transgene2vcf: \$(transgene2vcf.py -v)
    END_VERSIONS
    """
}
