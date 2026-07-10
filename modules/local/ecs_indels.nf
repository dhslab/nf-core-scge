// ECS_INDELS — deep, error-corrected edit calls at the hotspot panel (the TRUTH arm).
// Wraps bin/find_edited_reads.py directly (explicit --edited-bam/--control-bam, no length
// heuristic), producing <id>.offtarget_analysis.tsv whose indel_fraction is the ECS VAF.
process ECS_INDELS {
    tag "$meta.id"
    label 'process_highmem'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    tuple val(meta), val(edited_cram), val(control_cram), val(target_file)
    val reference

    output:
    tuple val(meta), path("${meta.id}.offtarget_analysis.tsv"), emit: indels_file
    tuple val(meta), path("${meta.id}.offtarget_edits.vcf"),    emit: indels_vcf
    path "versions.yml", emit: versions

    script:
    """
    find_edited_reads.py \\
        --fasta ${reference} \\
        --edited-bam ${edited_cram} \\
        --control-bam ${control_cram} \\
        --target-file ${target_file} \\
        -u ${meta.id}.unevaluable_reads.txt \\
        --vcf-out ${meta.id}.offtarget_edits.vcf \\
        -o ${meta.id}.offtarget_analysis.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.offtarget_analysis.tsv ${meta.id}.offtarget_edits.vcf versions.yml"
}
