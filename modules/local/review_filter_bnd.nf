// REVIEW_FILTER_BND — the breakend counterpart of REVIEW_FILTER.
//
// Wraps bin/review_filter_bnd.py. Reads the same *.offtarget_analysis.tsv tables as REVIEW_FILTER
// and triages the BND calls rather than the indels. On the 32-sample CAR-T WGS cohort (the
// 2026-08-17 run): 1,022 breakends -> 25 pass the evidence gate -> 25 reach the queue, i.e. rules
// 2-4 drop nothing and the >= 3 read gate is the only filter that acts. Those 25 rows are 8
// junctions (23 multi-cut deletion rows across 7 samples spanning 414 bp to 126 kb, 2 deletions at
// a cut site, 0 off-target junctions). See docs/OFFTARGET.md for the full accounting.
//
// The earlier "2,404 -> 36 -> 25" in this header described the 2026-08-10 run, before the caller
// changed; do not quote it against current output.
//
// The rules differ from the indel side because a junction has two ends. Length diversity is
// meaningless here, so it is replaced by breakpoint promiscuity: a bin joining many distinct
// partners is a mapping hub rather than biology, with the real cut sites exempt.
//
// The SV noise panel (--sv-noise) is off unless params.review_sv_noise is set, and only the WGS
// v3.1.0 panel is safe -- the IDPF and FF_Heme panels flag 25/25 and 24/25 of the real junctions
// respectively. See the warning in nextflow.config next to review_sv_noise.
process REVIEW_FILTER_BND {
    tag "cohort"
    label 'process_single'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path analysis_tsvs
    path sv_noise

    output:
    path "bnd_review_queue.tsv",     emit: queue
    path "bnd_review_queue_all.tsv", emit: audit
    path "versions.yml",             emit: versions

    script:
    def sv_arg = sv_noise.name != 'NO_FILE'
        ? "--sv-noise ${sv_noise} --sv-noise-slop ${params.review_sv_noise_slop}"
        : ''
    """
    python ${projectDir}/bin/review_filter_bnd.py ${analysis_tsvs} \\
        ${sv_arg} \\
        --min-reads ${params.review_bnd_min_reads} \\
        --max-cut-dist ${params.review_bnd_max_cut_dist} \\
        --max-partners ${params.review_bnd_max_partners} \\
        --bin-size ${params.review_bnd_bin_size} \\
        -o bnd_review_queue.tsv

    # Same thresholds, nothing filtered: every gated row with a why_dropped column, so a
    # reviewer can audit what was removed and why without re-running anything.
    python ${projectDir}/bin/review_filter_bnd.py ${analysis_tsvs} \\
        ${sv_arg} \\
        --min-reads ${params.review_bnd_min_reads} \\
        --max-cut-dist ${params.review_bnd_max_cut_dist} \\
        --max-partners ${params.review_bnd_max_partners} \\
        --bin-size ${params.review_bnd_bin_size} \\
        --keep-all \\
        -o bnd_review_queue_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // versions.yml must carry real content even in a stub: nf-core's processVersionsFromYAML
    // does yaml.load(f).collectEntries{...}, and an empty file loads as null -> NPE.
    """
    touch bnd_review_queue.tsv bnd_review_queue_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
