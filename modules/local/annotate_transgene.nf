process ANNOTATE_TRANSGENE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(files), path("${meta.id}.transgene_out.tsv")

    output:
    tuple val(meta), path("${meta.id}.transgene.annotated.tsv"), emit: annotated_transgenes
    path "versions.yml",    emit: versions

    script:
    """
    awk 'BEGIN {FS=OFS="\\t"} NR > 1 {sub(/^chr/, "", \$1); for (i = 1; i < NF; i++) {printf "%s%s", \$i, (i == NF-1 ? "\\n" : OFS)}}' ${meta.id}.transgene_out.tsv > ${meta.id}.edited.transgene.tsv && \\
    /opt/vep/src/ensembl-vep/vep --cache --dir /storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/VEP_cache --symbol --per_gene -i ${meta.id}.edited.transgene.tsv -o ${meta.id}.transgene.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    """

}