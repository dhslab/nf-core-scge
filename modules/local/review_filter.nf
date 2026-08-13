// REVIEW_FILTER — turn the raw off-target call table into the short list a human reviews.
//
// Wraps bin/review_filter.py. Six rules, first match wins: matched control clean, indel within
// 10 bp of a PAM position, >=3 distinct indel lengths, not a known-bad site (panel of normals),
// not inside a repeat, and not on a recurrent indel-noise locus in the external DRAGEN panel.
// On the 32-sample CAR-T WGS cohort this takes 1,498 gated rows to an 88-site queue (87 with
// rule 6 enabled) while retaining all 81 on-target edits.
//
// The panel of normals is optional but strongly preferred: it is what makes rule 4 work for a
// SINGLE-GUIDE submission. Without it the process falls back to cross-guide recurrence, which
// needs several differently-guided samples in the same invocation and cannot fire at all for
// one guide. All samples are staged into a single call so that fallback has a chance to work.
//
// Rule 6 (--snv-noise) is off unless params.review_snv_noise is set. It streams a ~1 GB unindexed
// BED, and this process invokes the script twice, so enabling it costs roughly 6 minutes.
process REVIEW_FILTER {
    tag "cohort"
    label 'process_single'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path analysis_tsvs
    path pon
    path repeat_beds
    path snv_noise

    output:
    path "review_queue.tsv",     emit: queue
    path "review_queue_all.tsv", emit: audit
    path "versions.yml",         emit: versions

    script:
    def pon_arg = (pon.name != 'NO_FILE' && params.review_noise_model == 'off') ? "--pon ${pon}" : ''
    def rep_arg = repeat_beds ? "--repeats ${repeat_beds.join(' ')}" : ''
    def snv_arg = snv_noise.name != 'NO_FILE'
        ? "--snv-noise ${snv_noise} --snv-noise-min-donors ${params.review_snv_noise_min_donors}"
        : ''
    // The noise model REPLACES rule 4, so the PoN is not passed alongside it -- handing the script
    // both would silently pick one and make the run's provenance unreadable.
    def nm_arg = params.review_noise_model != 'off'
        ? "--noise-model ${params.review_noise_model} --aq-min ${params.review_aq_min}" +
          (params.review_depth_floor ? '' : ' --no-depth-floor')
        : ''
    def strict_arg = params.review_strict_fallback ? '--strict-fallback' : ''
    """
    python ${projectDir}/bin/review_filter.py ${analysis_tsvs} \\
        ${pon_arg} ${rep_arg} ${snv_arg} ${nm_arg} ${strict_arg} \\
        --min-reads ${params.review_min_reads} \\
        --min-vaf ${params.review_min_vaf} \\
        --max-cut-dist ${params.review_max_cut_dist} \\
        --min-distinct-len ${params.review_min_distinct_len} \\
        --max-control-vaf ${params.review_max_control_vaf} \\
        -o review_queue.tsv

    # Same thresholds, nothing filtered: every gated row with a why_dropped column, so a
    # reviewer can audit what was removed and why without re-running anything.
    python ${projectDir}/bin/review_filter.py ${analysis_tsvs} \\
        ${pon_arg} ${rep_arg} ${snv_arg} ${nm_arg} ${strict_arg} \\
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
