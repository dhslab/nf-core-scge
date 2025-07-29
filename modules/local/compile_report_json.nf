process COMPILE_REPORT_JSON {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-clerbase:250719'

    input:
    tuple val(meta),
          path(cna_plot),
          path(baf_plot),
          path(circos_plot),
          path(on_target_sv_transgene),
          path(off_target_indels)

    output:
    tuple val(meta), path("report_input.json"), emit: json
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def transgene_str = meta.transgene ?: "N/A"
    def circos_arg = circos_plot.name == "NO_FILE.png" ? "" : "--circos_plot ${circos_plot}"
    """
    compile_report_data.py \\
        --sample_id ${meta.id} \\
        --transgene "${transgene_str}" \\
        --cna_plot ${cna_plot} \\
        --baf_plot ${baf_plot} \\
        \${circos_arg} \\
        --on_target_sv_transgene ${on_target_sv_transgene} \\
        --off_target_indels ${off_target_indels} \\
        --output report_input.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
} 