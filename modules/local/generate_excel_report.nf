process GENERATE_EXCEL_REPORT {
    tag "${meta.id}"
    label "process_low"
    label "final_output"
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(report_json), path(cna_plot), path(baf_plot), path(circos_plot), path(indel_plot), path(off_target_plot)

    output:
    path("${meta.id}.scge_report.xlsx"), emit: excel_report
    path("versions.yml")               , emit: versions

    script:
    """
    # Install xlsxwriter (pandas and openpyxl are in base image)
    # Use --user to avoid permission issues in /opt/conda
    pip install --user pandas xlsxwriter openpyxl
    
    # Ensure python can find the user installed packages
    export PYTHONPATH=\$(python3 -m site --user-site):\${PYTHONPATH:-}

    # Defensive copy/dereference logic for Circos
    CIRCOS_FILENAME=""
    if [ -s "${circos_plot}" ]; then
        CIRCOS_FILENAME=\$(basename "${circos_plot}")
        if [ -L "\${CIRCOS_FILENAME}" ]; then
            cp -L "\${CIRCOS_FILENAME}" "dereferenced_circos.png"
            mv "dereferenced_circos.png" "\${CIRCOS_FILENAME}"
        fi
    fi

    make_scge_excel.py \\
        --report_json "${report_json}" \\
        --circos_plot "\${CIRCOS_FILENAME}" \\
        --cna_plot "${cna_plot}" \\
        --baf_plot "${baf_plot}" \\
        --indel_freq_plot "${indel_plot}" \\
        --off_targets_plot "${off_target_plot}" \\
        --output "${meta.id}.scge_report.xlsx"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.scge_report.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}