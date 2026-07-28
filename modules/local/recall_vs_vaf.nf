// RECALL_VS_VAF — the honest limit, measured. From the training table, compute WGS
// recall of ECS-confirmed edits as a function of ECS VAF, and the VAF floor above which
// WGS-only detection is trustworthy. Curve + metrics CSV.
process RECALL_VS_VAF {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path training

    output:
    path "recall_vs_vaf.csv",             emit: metrics
    path "recall_vs_vaf.png", optional: true, emit: curve
    path "versions.yml",                  emit: versions

    script:
    """
    python ${projectDir}/bin/recall_vs_vaf.py \\
        --training ${training} \\
        --hi ${params.offtarget_hi_score} \\
        --target-recall ${params.offtarget_target_recall} \\
        --min-ecs-vaf ${params.offtarget_min_ecs_vaf} \\
        --min-ecs-reads ${params.offtarget_min_ecs_reads} \\
        --out-metrics recall_vs_vaf.csv \\
        --out-curve recall_vs_vaf.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch recall_vs_vaf.csv recall_vs_vaf.png versions.yml"
}
