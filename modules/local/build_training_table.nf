// BUILD_TRAINING_TABLE — join WGS hotspot features ⋈ ECS truth on (guide, chrom, start).
// Emits training.tsv (WGS features + ECS VAF + label) for the OFFLINE model trainer.
// Training is deliberately NOT in this DAG; the deployed model stays a fixed asset.
process BUILD_TRAINING_TABLE {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    path wgs_scores
    path truth
    path samplesheet

    output:
    path "training.tsv", emit: training
    path "versions.yml", emit: versions

    script:
    """
    join_training_table.py \\
        --wgs-scores ${wgs_scores} \\
        --truth ${truth} \\
        --samplesheet ${samplesheet} \\
        --out training.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch training.tsv versions.yml"
}
