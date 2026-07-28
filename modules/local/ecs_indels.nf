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
    tuple val(meta), path("${meta.id}.tagged.bam*"), optional: true, emit: tagged_bam
    path "versions.yml", emit: versions

    script:
    // The per-read "unevaluable reads" log is a debug artifact that is NOT an emitted
    // output and is not consumed downstream, yet it grows to ~0.1-1 TB per sample and
    // was the sole cause of multi-TB work-dir bloat / ENOSPC. Off unless explicitly asked.
    def unevaluable = params.offtarget_ecs_unevaluable_log ? "-u ${meta.id}.unevaluable_reads.txt" : ""
    // Review aid, not a pipeline input: a BAM of the target windows in which every read
    // carries an XC tag naming the per-read call, for colouring the pileup in IGV. Window-
    // restricted, but it still scales with target count -- off unless explicitly asked.
    def tagged_bam = params.offtarget_tagged_bam ? "--tagged-bam-out ${meta.id}.tagged.bam" : ""
    """
    python ${projectDir}/bin/find_edited_reads.py \\
        --fasta ${reference} \\
        --edited-bam ${edited_cram} \\
        --control-bam ${control_cram} \\
        --target-file ${target_file} \\
        ${unevaluable} \\
        ${tagged_bam} \\
        --vcf-out ${meta.id}.offtarget_edits.vcf \\
        -o ${meta.id}.offtarget_analysis.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def tagged_bam = params.offtarget_tagged_bam ? "${meta.id}.tagged.bam ${meta.id}.tagged.bam.bai" : ""
    "touch ${meta.id}.offtarget_analysis.tsv ${meta.id}.offtarget_edits.vcf ${tagged_bam} versions.yml"
}
