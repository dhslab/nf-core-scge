process ANNOTATE_TRANSGENE_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), path(transgene_vcf)
    path(reference)
    path(vepcache)
    path(cytobands)

    output:
    tuple val(meta), path("${meta.id}.transgene_junctions.annotated.tsv"), emit: annotated_junctions
    path "versions.yml",    emit: versions

    script:
    def annotate_args = [
        vepcache                             ? "--dir ${vepcache}"   : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta ${it}" } ?: "",
        cytobands.find{ it ==~ /.*\.bed\.gz$/ }?.with{ "--plugin StructuralVariantOverlap,file=${it}" } ?: "",
        transgene_vcf ? "-i ${transgene_vcf}" : ""
    ].join(' ').trim()

    """
    /opt/vep/src/ensembl-vep/vep \\
        --offline \\
        --cache \\
        --symbol \\
        --term SO \\
        --flag_pick \\
        --format vcf \\
        --tab \\
        --fields Location,Consequence,SYMBOL,BIOTYPE,EXON,INTRON,STRAND,Canonical,Pick,Feature \\
        ${annotate_args} \\
        -o ${meta.id}.transgene_junctions.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\\s*//g')
    END_VERSIONS
    """
}