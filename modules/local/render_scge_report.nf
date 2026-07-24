process RENDER_SCGE_REPORT {
    tag "${meta.id}"
    label "process_high"
    label 'final_output'
    container "ghcr.io/dhslab/docker-quarto-chromoseq:latest"

    input:
    tuple val(meta), path(report_json), path(cna_plot), path(baf_plot), path(circos_plot)
    path(scge_report_qmd)

    output:
    tuple val(meta), path("${meta.id}.scge_report.html") , emit: scge_report
    tuple val(meta), path("indel_freq.png"), emit: indel_plot
    tuple val(meta), path("off_targets.png"), emit: off_targets_plot
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

    # Create params file for Quarto to avoid shell quoting issues
    # Ensure circos_plot is just the filename if it exists
    CIRCOS_FILENAME=""
    if [ -s "${circos_plot}" ]; then
        CIRCOS_FILENAME=\$(basename "${circos_plot}")
        # Resolve symlink to physical file to prevent container mount issues
        if [ -L "\${CIRCOS_FILENAME}" ]; then
            cp -L "\${CIRCOS_FILENAME}" "dereferenced_circos.png"
            mv "dereferenced_circos.png" "\${CIRCOS_FILENAME}"
        fi
    fi

    # Copy QMD to current directory to ensure relative paths work correctly
    # (Quarto can sometimes resolve symlinked QMDs to their original dir)
    cp "${scge_report_qmd}" "report.qmd"

    # Copy all input images to current directory if not already there
    # This ensures they are in the same directory as the .qmd file being rendered
    if [ -s "${cna_plot}" ]; then
        CNA_BASENAME=\$(basename "${cna_plot}")
        if [ ! -f "\${CNA_BASENAME}" ]; then
            cp -L "${cna_plot}" "\${CNA_BASENAME}"
        fi
    fi
    if [ -s "${baf_plot}" ]; then
        BAF_BASENAME=\$(basename "${baf_plot}")
        if [ ! -f "\${BAF_BASENAME}" ]; then
            cp -L "${baf_plot}" "\${BAF_BASENAME}"
        fi
    fi
    # Also copy off-target plots if they exist
    if [ -f "indel_freq.png" ]; then
        cp -L "indel_freq.png" "indel_freq_copy.png"
        mv "indel_freq_copy.png" "indel_freq.png"
    fi
    if [ -f "off_targets.png" ]; then
        cp -L "off_targets.png" "off_targets_copy.png"
        mv "off_targets_copy.png" "off_targets.png"
    fi

    cat <<EOF > params.yml
    report_json: "${report_json}"
    off_target_threshold: ${params.off_target_threshold}
    cna_plot_file: "${cna_plot}"
    baf_plot_file: "${baf_plot}"
    circos_plot_file: "\${CIRCOS_FILENAME}"
    transgene_name: "${params.transgene_name ?: ''}"
    EOF

    # Explicitly set execution dir to current dir
    quarto render "report.qmd" --execute-dir . --execute-params params.yml --output "${meta.id}.scge_report.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.scge_report.html
    touch indel_freq.png
    touch off_targets.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
