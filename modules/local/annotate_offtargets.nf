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
    /opt/vep/src/ensembl-vep/vep \\
            --offline \\
            --cache \\
            ${vep_args} \\
            --force --symbol --term SO --per_gene --fields "Location,SYMBOL,DISTANCE,INTRON,EXON" --numbers -o stdout \\
            | perl -F'\t' -ane 'next if /^#/; @id=split(/,/,\$F[0]); \$xtra=\$F[13]; \$sym = (\$xtra =~ /SYMBOL=([^;]+)/)[0] // "."; \$dist = (\$xtra =~ /DISTANCE=([^;]+)/)[0] // "."; \\ 
                                \$int = (\$xtra =~ /INTRON=([^;]+)/)[0] // "."; \$ex = (\$xtra =~ /EXON=([^;]+)/)[0] // "."; \\
                                print join("\t", \$id[3], \$id[4]-1, \$id[4], "INS", \$id[5], \$F[0]) . ",SYMBOL=\$sym;GeneID=\$F[3];TranscriptID=\$F[4],Distance=\$dist;Intron=\$int;Exon=\$ex;Consequence=\$F[6];\n"' \\
            > "${meta.id}.indels.annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: echo "foo"
    END_VERSIONS
    """
}


