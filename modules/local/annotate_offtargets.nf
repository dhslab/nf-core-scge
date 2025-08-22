process ANNOTATE_OFFTARGETS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(hotspots_csv)

    output:
    tuple val(meta), path("${meta.id}.offtargets.annotated.tsv"), emit: annotated_offtargets
    path "versions.yml",    emit: versions

    script:
    """
    set -euo pipefail

    # Convert hotspots CSV to region list acceptable by VEP (chrom start end)
    # Expecting columns: Chromosome,Start,End,...
    awk -F, 'NR>1 {gsub(/^chr/, "", \$1); print \$1"\t"\$2"\t"\$3 }' ${hotspots_csv} > ${meta.id}.offtargets.regions.tsv

    /opt/vep/src/ensembl-vep/vep \
        --offline \
        --cache --dir ${params.vep_cache} \
        --fasta ${params.fasta} \
        --symbol --term SO --flag_pick \
        --tab \
        --fields Location,Consequence,SYMBOL,BIOTYPE,EXON,INTRON,STRAND,Canonical,Pick \
        --format region \
        -i ${meta.id}.offtargets.regions.tsv \
        -o ${meta.id}.offtargets.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    """
}


