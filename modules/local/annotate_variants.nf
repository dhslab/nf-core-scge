process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(vcf_file)
    path(reference)
    path(vep_cache)

    output:
    tuple val(meta), path("${meta.id}.hard-filtered.annotated.vcf.gz"), optional: true, emit: vcf
    tuple val(meta), path("${meta.id}.hard-filtered.annotated.vcf.gz.tbi"), optional: true, emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
        --format vcf \\
        --vcf \\
        --fasta ${reference} \\
        --hgvs \\
        --symbol \\
        --term SO \\
        --flag_pick \\
        --force_overwrite \\
        -i ${vcf_file} \\
        --offline \\
        --cache \\
        --max_af \\
        --dir ${vep_cache} \\
        -o "${meta.id}.hard-filtered.annotated.vcf"

    bgzip -c ${meta.id}.hard-filtered.annotated.vcf > ${meta.id}.hard-filtered.annotated.vcf.gz
    tabix -p vcf ${meta.id}.hard-filtered.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.hard-filtered.annotated.vcf.gz
    touch ${meta.id}.hard-filtered.annotated.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: "105.0"
    END_VERSIONS
    """
}