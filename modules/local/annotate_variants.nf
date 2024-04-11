process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
    --format vcf --vcf --fasta ${params.fasta} --hgvs --symbol --term SO --flag_pick -o ${meta.id}.annotated.vcf \\
    -i ${meta.id}.sv.vcf.gz --custom ${params.cytobands},cytobands,bed --custom ${params.custom_annotation_vcf},${params.custom_annotation_parameters} --offline --cache --max_af --dir ${params.vepcache} && \\
    bgzip -c ${meta.id}.annotated.vcf > ${meta.id}.annotated.vcf.gz && \\
    tabix -p vcf ${meta.id}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    /opt/vep/src/ensembl-vep/vep --dir ${params.vepcache} --show_cache_info | awk '{ print "    "\$1": "\$2; }' >> versions.yml
    """

}