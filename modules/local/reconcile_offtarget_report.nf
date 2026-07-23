// RECONCILE_OFFTARGET_REPORT — annotate the genome-wide PoN worklist with is_hotspot /
// ecs_confirmed / ecs_if, so a reviewer sees whether a homology-free WGS hit is a known
// predicted (ECS-backed) site or a novel candidate. Truth is optional (wgs_only mode).
process RECONCILE_OFFTARGET_REPORT {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path worklist
    path truth      // ecs_hotspot_truth.csv, or the NO_FILE placeholder

    output:
    path "offtarget_report.csv", emit: report
    path "versions.yml",         emit: versions

    script:
    def truth_arg = truth.name != 'NO_FILE' ? "--truth ${truth}" : ""
    """
    python ${projectDir}/bin/reconcile_offtarget_report.py \\
        --worklist ${worklist} \\
        ${truth_arg} \\
        --pad ${params.offtarget_hotspot_pad} \\
        --out offtarget_report.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch offtarget_report.csv versions.yml"
}
