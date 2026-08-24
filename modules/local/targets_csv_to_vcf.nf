// TARGETS_CSV_TO_VCF — expand a combined off-target sites CSV into the per-base hotspot
// VCF the OFFTARGET arm consumes as a target_file (windowed, merged, reference-ordered).
process TARGETS_CSV_TO_VCF {
    tag "${meta.id}"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    tuple val(meta), path(sites_csv)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.targets.vcf"), emit: vcf
    path "versions.yml",                             emit: versions

    script:
    """
    python ${projectDir}/bin/targets_csv_to_vcf.py \\
        --csv ${sites_csv} \\
        --fasta ${fasta} \\
        --window ${params.hotspot_window_size} \\
        -o ${meta.id}.targets.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.targets.vcf versions.yml"
}
