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
    path "versions.yml",           emit: versions

    script:
    """
    score.py \\
        --table ${table} \\
        --cram-list ${cram_map} \\
        --ref ${reference} \\
        --model ${model} \\
        --no-homology-gate \\
        --include-ontarget \\
        --min-ifrac -1 \\
        --max-control 2 \\
        --min-span ${params.offtarget_min_span} \\
        --out wgs_hotspot_scores.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch wgs_hotspot_scores.csv versions.yml"
}
