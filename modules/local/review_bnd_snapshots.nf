// REVIEW_BND_SNAPSHOTS — the breakend half of the review packet.
//
// One PNG per JUNCTION, not per queue row: bnd_review_queue.tsv reports each event from
// both ends and at a few bp of jitter, so the rows outnumber the junctions several-fold
// (25 rows / 8 junctions on the CAR-T cohort). Each figure is a to-scale schematic of
// the excision over a 2x2 read grid — left and right breakpoint, each with the edited
// sample above its matched unedited control — with reads whose SA tag lands at the
// partner locus drawn green.
//
// Same CRAM-map-by-absolute-path convention as REVIEW_SNAPSHOTS, and the same reason:
// the queue names arbitrary samples, and staging every cohort CRAM to draw a handful of
// pictures would copy terabytes.
process REVIEW_BND_SNAPSHOTS {
    tag "cohort"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path bnd_queue
    path cram_map
    val  reference

    output:
    path "bnd_snapshots/*.png", optional: true, emit: snapshots
    path "versions.yml",                        emit: versions

    script:
    """
    python ${projectDir}/bin/bnd_snapshots.py \\
        --queue ${bnd_queue} \\
        --cram-map ${cram_map} \\
        --fasta ${reference} \\
        --window ${params.review_bnd_snapshot_window} \\
        --max-junctions ${params.review_bnd_max_junctions} \\
        --outdir bnd_snapshots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // versions.yml must carry real content even in a stub: nf-core's processVersionsFromYAML
    // does yaml.load(f).collectEntries{...}, and an empty file loads as null -> NPE.
    """
    mkdir -p bnd_snapshots && touch bnd_snapshots/stub_bnd_junction.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
