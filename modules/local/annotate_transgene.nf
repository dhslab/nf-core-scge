process ANNOTATE_TRANSGENE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), path(transgene_vcf)
    path(fasta)
    path(vepcache)
    path(cytobands)


    output:
    tuple val(meta), path("${meta.id}.transgene.annotated.tsv"), emit: annotated_transgene_variants
    path "versions.yml",    emit: versions

    script:
    """
    if [ \$(grep -vc '^#' ${transgene_vcf}) -gt 0 ]; then
        /opt/vep/src/ensembl-vep/vep \\
            --offline \\
            --cache \\
            --dir ${vepcache} \\
            --fasta ${fasta} \\
            --symbol \\
            --term SO \\
            --flag_pick \\
            --everything \\
            --tab \\
            --fields Location,Consequence,SYMBOL,BIOTYPE,EXON,INTRON,STRAND,Canonical,Pick,Feature \\
            --plugin StructuralVariantOverlap,file=${cytobands} \\
            -i ${transgene_vcf} \\
            -o ${meta.id}.transgene.annotated.tsv
    else
        touch ${meta.id}.transgene.annotated.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\\s*//g')
    END_VERSIONS
    """
}