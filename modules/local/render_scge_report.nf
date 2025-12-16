process RENDER_SCGE_REPORT {
    tag "${meta.id}"
    label "process_high"
    label 'final_output'
    container "ghcr.io/dhslab/docker-quarto-chromoseq:latest"

    input:
    tuple val(meta), path(report_json), path(cna_plot), path(baf_plot)
    path(scge_report_qmd)

    output:
    path("${meta.id}.scge_report.html") , emit: scge_report
    path("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Initialize micromamba environment (container uses micromamba, not conda)
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate base

    # Add quarto to PATH
    export PATH="/opt/quarto/bin:\$PATH"

    # Debug: verify tools are available
    which quarto && which Rscript || { echo "ERROR: Required tools not found"; exit 1; }

    quarto render ${scge_report_qmd} -P report_json:"${report_json}" -P off_target_threshold:${params.off_target_threshold} --output "${meta.id}.scge_report.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
