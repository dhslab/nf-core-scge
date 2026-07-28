// OFFTARGET_METRICS — the metrics a recall-first diagnostic is judged on: PR-AUC and
// recall-weighted F-beta (F2, F5), plus precision/recall/F1 at the reported operating
// point. Computed against the ECS label in training.tsv, which has a genuine two-class
// label — NEVER against the manual review, which has no confirmed negatives (its NaNs
// mean "unreviewed", not "rejected"). Recall vs manual review stays with
// bin/validate_recall.py and stays recall-only. Both denominators are stamped into the
// outputs so the two can never be confused.
process OFFTARGET_METRICS {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path training

    output:
    path "offtarget_metrics.json",                     emit: metrics
    path "offtarget_metrics.txt",                      emit: report
    path "offtarget_pr_curve.png", optional: true,     emit: curve
    path "versions.yml",                               emit: versions

    script:
    """
    python ${projectDir}/bin/offtarget_metrics.py \\
        --training ${training} \\
        --hi ${params.offtarget_hi_score} \\
        --min-ecs-vaf ${params.offtarget_min_ecs_vaf} \\
        --min-ecs-reads ${params.offtarget_min_ecs_reads} \\
        --betas ${params.offtarget_metrics_betas} \\
        --negatives ${params.offtarget_metrics_negatives} \\
        --out-json offtarget_metrics.json \\
        --out-txt offtarget_metrics.txt \\
        --out-curve offtarget_pr_curve.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch offtarget_metrics.json offtarget_metrics.txt offtarget_pr_curve.png versions.yml"
}
