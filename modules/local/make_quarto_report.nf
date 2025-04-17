process MAKE_QUARTO_REPORT {
    tag "${meta.id}"
    label "process_low"
    label 'final_output'
    container "ghcr.io/dhslab/docker-quarto-chromoseq:latest"

    input:
    tuple val(meta), path(circos_plot), path(indels)
    path(scge_report_qmd)

    output:
    path("${meta.id}.scge_report.html") , emit: scge_report
    path("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export PATH=\$PATH:/opt/conda/bin/

    quarto render make_scge_report.qmd -P offtargets:"$indels" -P transgene:"${params.transgene}" --output "${meta.id}.scge_report.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}