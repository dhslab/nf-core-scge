// PON_OFFTARGET_FILTER — Panel-of-Normals post-filter (the shallow-normal fix).
// Wraps bin/pon_filter.py. Pools ALL unedited normals in the cram map and demotes any
// LIKELY-EDIT candidate that shows the same indel in ANY donor's normal (germline/mosaic
// leak the single matched normal was too shallow to catch). Genuine gather step.
process PON_OFFTARGET_FILTER {
    tag "cohort"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    path worklist
    path cram_map
    val  reference

    output:
    path "wgs_offtarget_worklist_pon.csv", emit: worklist
    path "versions.yml",                   emit: versions

    script:
    """
    pon_filter.py \\
        --worklist ${worklist} \\
        --cram-list ${cram_map} \\
        --ref ${reference} \\
        --top ${params.offtarget_top} \\
        --out wgs_offtarget_worklist_pon.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch wgs_offtarget_worklist_pon.csv versions.yml"
}
