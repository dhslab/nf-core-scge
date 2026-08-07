// REVIEW_FILTER — turn the raw off-target call table into the short list a human reviews.
//
// Wraps bin/review_filter.py. Four rules: matched control clean, indel within 10 bp of a PAM
// position, >=3 distinct indel lengths, and not a known-bad site. On the CAR-T WGS cohort this
// takes the review queue from 238 sites to 63 while retaining 61/61 real edits.
//
// The panel of normals is optional but strongly preferred: it is what makes rule 4 work for a
// SINGLE-GUIDE submission. Without it the process falls back to cross-guide recurrence, which
// needs several differently-guided samples in the same invocation and cannot fire at all for
// one guide. All samples are staged into a single call so that fallback has a chance to work.
process REVIEW_FILTER {
    tag "cohort"
    label 'process_single'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path analysis_tsvs
    path pon
    path repeat_beds

    output:
    path "review_queue.tsv",     emit: queue
    path "review_queue_all.tsv", emit: audit
    path "versions.yml",         emit: versions

    script:
    def pon_arg = pon.name != 'NO_FILE' ? "--pon ${pon}" : ''
    def rep_arg = repeat_beds ? "--repeats ${repeat_beds.join(' ')}" : ''
    """
    python ${projectDir}/bin/review_filter.py ${analysis_tsvs} \\
        ${pon_arg} ${rep_arg} \\
        --min-reads ${params.review_min_reads} \\
        --min-vaf ${params.review_min_vaf} \\
        --max-cut-dist ${params.review_max_cut_dist} \\
        --min-distinct-len ${params.review_min_distinct_len} \\
        --max-control-vaf ${params.review_max_control_vaf} \\
        -o review_queue.tsv

    # Same thresholds, nothing filtered: every gated row with a why_dropped column, so a
    # reviewer can audit what was removed and why without re-running anything.
    python ${projectDir}/bin/review_filter.py ${analysis_tsvs} \\
        ${pon_arg} ${rep_arg} \\
        --min-reads ${params.review_min_reads} \\
        --min-vaf ${params.review_min_vaf} \\
        --max-cut-dist ${params.review_max_cut_dist} \\
        --min-distinct-len ${params.review_min_distinct_len} \\
        --max-control-vaf ${params.review_max_control_vaf} \\
        --keep-all \\
        -o review_queue_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // versions.yml must carry real content even in a stub: nf-core's processVersionsFromYAML
    // does yaml.load(f).collectEntries{...}, and an empty file loads as null -> NPE.
    """
    touch review_queue.tsv review_queue_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
