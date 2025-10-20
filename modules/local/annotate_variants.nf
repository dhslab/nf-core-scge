process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), path(dragen_dir, stageAs: 'dragen/*'), path(reference), path(vep_cache)

    output:
    tuple val(meta), path("*.annotated.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    VCF_FILE=\$(find dragen -name "*hard-filtered.vcf.gz" | head -n 1)

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --hgvs \\
        --cache \\
        --cache_version 113 \\
        --species homo_sapiens \\
        --max_af \\
        --symbol \\
        --term SO \\
        --offline \\
        --flag_pick \\
        --format vcf \\
        --force_overwrite \\
        --dir "${vep_cache}" \\
        --fasta ${reference} \\
        -o "${meta.id}.hard-filtered.annotated.vcf" \\
        \$VCF_FILE

    bgzip -c ${meta.id}.hard-filtered.annotated.vcf > ${meta.id}.hard-filtered.annotated.vcf.gz
    tabix -p vcf ${meta.id}.hard-filtered.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """

    stub:
    def vcf_file = dragen_dir.find{ it ==~ /.*\\.hard-filtered.vcf.gz$/ }
    def annotate_args = [
        vep_cache                                 ? "--dir vep_cache"   : "",
        reference.find{ it ==~ /.*\\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
        vcf_file ? "-i dragen/${vcf_file.name}" : ""
    ].join(' ').trim()
    """
    touch ${meta.id}.hard-filtered.annotated.vcf.gz
    touch ${meta.id}.hard-filtered.annotated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """
}