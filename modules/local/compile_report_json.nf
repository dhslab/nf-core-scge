process COMPILE_REPORT_JSON {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-baseimage:latest'

    input:
    tuple val(meta), path(files)
    val(timestamp)

    output:
    tuple val(meta), path("${meta.id}.scge_report.json"), emit: json
    path "versions.yml"                                 , emit: versions


    script:
    def args = task.ext.args ?: ''
    def control_sample = params.control_sample ?: (meta.normal ?: "N/A")
    def input = [
        files.find{ it ==~ /.*\.wgs_overall_mean_cov_tumor\.csv$/ }?.with{ "--tumor_coverage $it" } ?: "",
        files.find{ it ==~ /.*\.wgs_overall_mean_cov_normal\.csv$/ }?.with{ "--normal_coverage $it" } ?: "",
        files.find{ it ==~ /.*\.cna_plot\.png$/ }?.with{ "--cna_plot $it" } ?: "",
        files.find{ it ==~ /.*\.baf_plot\.png$/ }?.with{ "--baf_plot $it" } ?: "",
        files.find{ it ==~ /.*\.annotated_transgene_insertions\.tsv$/ }?.with{ "--transgene_insertions $it" } ?: "",
        files.find{ it ==~ /.*\.transgene_insertions_circos\.png$/ }?.with{ "--circos_plot $it" } ?: "",
        files.find{ it ==~ /.*\.hard-filtered\.annotated\.tsv$/ }?.with{ "--somatic_variants $it" } ?: "",
        files.find{ it ==~ /.*\.offtarget_analysis\.tsv$/ }?.with{ "--offtarget_indels $it" } ?: "",
        files.find{ it ==~ /.*\.offtarget_svs\.tsv$/ }?.with{ "--offtarget_svs $it" } ?: "",
    ].join(' ').trim()

    """
    compile_report_data.py \\
        --sample_id ${meta.id} \\
        --control_sample "${control_sample}" \\
        ${input} \\
        --output ${meta.id}.scge_report.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
} 