// HOTSPOT_TO_TABLE — turn the ECS hotspot panel (+ ECS truth VAF) into a score.py input
// table keyed to the WGS samples of the same guide. Emits the scoring table and the ECS
// truth table for the downstream training-table / recall joins.
process HOTSPOT_TO_TABLE {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path ecs_tables      // all ECS *.offtarget_analysis.tsv
    path samplesheet

    output:
    path "wgs_hotspot_input_table.csv", emit: table
    path "ecs_hotspot_truth.csv",       emit: truth
    path "versions.yml",                emit: versions

    script:
    """
    hotspot_to_table.py \\
        --ecs-tables ${ecs_tables} \\
        --samplesheet ${samplesheet} \\
        --edit-threshold ${params.offtarget_ecs_edit_threshold} \\
        --out-table wgs_hotspot_input_table.csv \\
        --out-truth ecs_hotspot_truth.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch wgs_hotspot_input_table.csv ecs_hotspot_truth.csv versions.yml"
}
