process ANNOTATE_OFFTARGETS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep_release113:250810" // Use the same VEP container

    input:
    tuple val(meta), path(targetfile)
    path(vep_cache)
    path(fasta)

    output:
    tuple val(meta), path("${meta.id}.indels.annotated.tsv"), emit: targetfile
    path "versions.yml", emit: versions

    script:
    def vep_args = [
        targetfile                                          ? "-i ${targetfile}"    : "",
        vep_cache                                           ? "--dir ${vep_cache}"   : "",
        fasta.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: ""
    ].join(' ').trim()

    """
    head -n 1 ${targetfile} | tr -d '\\r' | sed 's/\$/,vep/' > "${meta.id}.indels.annotated.tsv"

    /opt/vep/src/ensembl-vep/vep \\
            --offline \\
            --cache \\
            ${vep_args} \\
            --force --symbol --term SO --per_gene --fields "Location,SYMBOL,DISTANCE,INTRON,EXON" --numbers -o stdout \\
            | add_vep2targetfile.pl >> "${meta.id}.indels.annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """
}


