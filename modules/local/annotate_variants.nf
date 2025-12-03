process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep_release113:250810"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*")
    path(reference)
    path(vep_cache)

    output:
    tuple val(meta), path("*.annotated.vcf.gz*"), emit: vcf
    path("versions.yml")                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def annotate_args = [
        vep_cache                                 ? "--dir ${vep_cache}"   : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
        dragen_files.find{ it ==~ /.*\.hard-filtered.vcf.gz$/ }?.with{ "-i $it" } ?: ""
    ].join(' ').trim()
    
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --hgvs \\
        --cache \\
        --max_af \\
        --symbol \\
        --term SO \\
        --offline \\
        --flag_pick \\
        --format vcf \\
        --force_overwrite \\
        ${annotate_args} \\
        -o "${meta.id}.hard-filtered.annotated.vcf"

    bgzip -c ${meta.id}.hard-filtered.annotated.vcf > ${meta.id}.hard-filtered.annotated.vcf.gz
    tabix -p vcf ${meta.id}.hard-filtered.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """

    stub:
    def annotate_args = [
        vep_cache                                 ? "--dir ${vep_cache}"   : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
        dragen_files.find{ it ==~ /.*\.hard-filtered.vcf.gz$/ }?.with{ "-i $it" } ?: ""
    ].join(' ').trim()
    
    """
    touch ${meta.id}.hard-filtered.annotated.vcf.gz
    touch ${meta.id}.hard-filtered.annotated.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: "105.0"
    END_VERSIONS
    """
}