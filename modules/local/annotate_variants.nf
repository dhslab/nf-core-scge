process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(files)
    tuple val(assay_inputs), path("*")
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.hard-filtered.annotated.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
    --format vcf --vcf --fasta ${reference} --hgvs --symbol --term SO --flag_pick -o ${meta.id}.hard-filtered.annotated.vcf \\
    -i ${meta.id}.hard-filtered.vcf.gz --offline --cache --max_af --dir ${assay_inputs.vepcache} && \\
    bgzip -c ${meta.id}.hard-filtered.annotated.vcf > ${meta.id}.hard-filtered.annotated.vcf.gz && \\
    tabix -p vcf ${meta.id}.hard-filtered.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    /opt/vep/src/ensembl-vep/vep --dir ${assay_inputs.vepcache} --show_cache_info | awk '{ print "    "\$1": "\$2; }' >> versions.yml
    """

}