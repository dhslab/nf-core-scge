// REVIEW_SNAPSHOTS — turn the review queue into a packet of IGV-style pileup images.
//
// One PNG per surviving site, edited sample on top and its matched unedited control below.
// That pairing is what makes the packet self-adjudicating: germline and shared alignment
// artifacts appear in both panels, a real edit appears in only one.
//
// CRAMs are read from an absolute-path map rather than staged, matching the convention
// SCORE_HOTSPOTS already uses for its cram list: the queue names arbitrary samples, so staging
// every cohort CRAM into this one task would copy terabytes to render a few dozen pictures.
process REVIEW_SNAPSHOTS {
    tag "cohort"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path review_queue
    path cram_map
    val  reference

    output:
    path "snapshots/*.png", optional: true, emit: snapshots
    path "versions.yml",                   emit: versions

    script:
    """
    python ${projectDir}/bin/review_snapshots.py \\
        --queue ${review_queue} \\
        --cram-map ${cram_map} \\
        --fasta ${reference} \\
        --window ${params.review_snapshot_window} \\
        --outdir snapshots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // versions.yml must carry real content even in a stub: nf-core's processVersionsFromYAML
    // does yaml.load(f).collectEntries{...}, and an empty file loads as null -> NPE.
    """
    mkdir -p snapshots && touch snapshots/stub_review_site.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
