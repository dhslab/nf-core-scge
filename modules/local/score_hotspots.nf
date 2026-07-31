// SCORE_HOTSPOTS — score the WGS CRAMs at the ECS hotspot panel with the shape model.
// Wraps bin/score.py with the homology gate OFF and the stage-1 if/control gates opened
// (the table's indel_fraction is a gate-passer; real signal is recomputed from the WGS
// CRAM). Output is the per-hotspot WGS feature+score table.
process SCORE_HOTSPOTS {
    tag "cohort"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path table
    path cram_map
    path model
    val  reference

    output:
    path "wgs_hotspot_scores.csv", emit: scores
    // Per-read XC-tagged BAMs for IGV, one per WGS sample, covering the called sites only.
    // Optional because --offtarget_wgs_tagged_bam is off by default.
    path "wgs_tagged/*.bam*", optional: true, emit: tagged_bam
    path "versions.yml",           emit: versions

    script:
    """
    python ${projectDir}/bin/score.py \\
        --table ${table} \\
        --cram-list ${cram_map} \\
        --ref ${reference} \\
        --model ${model} \\
        --no-homology-gate \\
        --include-ontarget \\
        --min-ifrac -1 \\
        --max-control 2 \\
        --min-span ${params.offtarget_min_span} \\
        ${params.offtarget_max_cut_dist != null ? "--max-cut-dist ${params.offtarget_max_cut_dist}" : ''} \\
        ${params.offtarget_wgs_tagged_bam ? "--tagged-bam-dir wgs_tagged" : ''} \\
        ${params.offtarget_rescue ? "--rescue-min-ifrac ${params.offtarget_rescue_min_ifrac} --rescue-min-conc ${params.offtarget_rescue_min_conc} --rescue-min-span ${params.offtarget_rescue_min_span}" : '--no-rescue'} \\
        --out wgs_hotspot_scores.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // the tagged BAMs are an `optional:` output, but the stub must still create them when the
    // param is on, or a stub run silently exercises a different DAG than the real one
    def tagged = params.offtarget_wgs_tagged_bam
        ? "mkdir -p wgs_tagged && touch wgs_tagged/stub.wgs_tagged.bam wgs_tagged/stub.wgs_tagged.bam.bai"
        : "true"
    """
    touch wgs_hotspot_scores.csv versions.yml
    ${tagged}
    """
}
