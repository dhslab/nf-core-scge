process ANNOTATE_TRANSGENE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), path(files), path(transgene_vcf)

    output:
    tuple val(meta), path("${meta.id}.transgene.annotated.tsv"), emit: annotated_transgene_variants
    path "versions.yml",    emit: versions

    script:
    """
    /opt/vep/src/ensembl-vep/vep \\
        --offline \\
        --cache \\
        --dir ${params.vep_cache} \\
        --fasta ${params.fasta} \\
        --symbol \\
        --term SO \\
        --flag_pick \\
        --everything \\
        --tab \\
        --fields Location,Consequence,SYMBOL,BIOTYPE,EXON,INTRON,STRAND,Canonical,Pick,Feature \\
        --plugin StructuralVariantOverlap \\
        -i ${transgene_vcf} \\
        -o ${meta.id}.transgene.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\\s*//g')
    END_VERSIONS
    """

}