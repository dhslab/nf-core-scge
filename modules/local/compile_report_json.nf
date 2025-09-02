process COMPILE_REPORT_JSON {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-baseimage:latest'

    input:
    tuple val(meta),
          path(cna_plot),
          path(baf_plot),
          val(circos_plot),
          path(on_target_sv_transgene),
          path(vcf_tsv),
          path(off_target_indels),
          path(tumor_cov),
          path(normal_cov),
          val(timestamp)

    output:
    tuple val(meta), path("report_input.json"), emit: json
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def transgene_str = params.transgene ?: (meta.transgene ?: "N/A")
    def control_sample = params.control_sample ?: (meta.normal ?: "N/A")
    def grnas_str = params.grnas ?: meta.id
    def hotspot_file_arg = meta.hotspot_file ? "--hotspot_file ${meta.hotspot_file}" : ""
    def circos_arg = (circos_plot && circos_plot.toString().endsWith(".png") && file(circos_plot).exists()) ? "--circos_plot ${circos_plot}" : ""
    def tumor_coverage_arg = tumor_cov ? "--tumor_coverage ${tumor_cov}" : ""
    def normal_coverage_arg = normal_cov ? "--normal_coverage ${normal_cov}" : ""
    """
    python3 ${projectDir}/bin/compile_report_data.py \\
        --sample_id ${meta.id} \\
        --transgene "${transgene_str}" \\
        --control_sample "${control_sample}" \\
        --grnas "${grnas_str}" \\
        ${hotspot_file_arg} \\
        --cna_plot ${cna_plot} \\
        --baf_plot ${baf_plot} \\
        ${circos_arg} \\
        --on_target_sv_transgene ${on_target_sv_transgene} \\
        --vcf_tsv ${vcf_tsv} \\
        --off_target_indels ${off_target_indels} \\
        ${tumor_coverage_arg} \\
        ${normal_coverage_arg} \\
        --output report_input.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
} 