process GENERATE_CNA_BAF_PLOTS {
    tag "${meta.id}"
    label 'process_low'
    container 'ghcr.io/dhslab/docker-rbase4.4.0:251223'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*")

    output:
    tuple val(meta), path("${meta.id}.cna_plot.png"), path("${meta.id}.baf_plot.png"), emit: plots
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input = [
        "--id ${meta.id}",
        dragen_files.find{ it ==~ /.*\.baf\.bedgraph\.gz$/ }?.with{ "--baf $it" } ?: "",
        dragen_files.find{ it ==~ /.*\.tn\.tsv\.gz$/ }?.with{ "--cn $it" } ?: ""
    ].join(' ').trim()
    """
    generate_cna_baf_plots.R ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_cna_baf_plots: \$(generate_cna_baf_plots.R --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cna_plot.png
    touch ${meta.id}.baf_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_cna_baf_plots: \$(generate_cna_baf_plots.R --version)
    END_VERSIONS
    """
}