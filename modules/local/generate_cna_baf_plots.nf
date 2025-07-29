process GENERATE_CNA_BAF_PLOTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-clerbase:250719'

    input:
    tuple val(meta), path(dragen_baf), path(dragen_cnv)

    output:
    tuple val(meta), path("cna_plot.png"), emit: cna_plot
    tuple val(meta), path("baf_plot.png"), emit: baf_plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    generate_cna_baf_plots.R \\
        ${meta.id} \\
        ${dragen_baf} \\
        ${dragen_cnv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed 's/.*version \\([0-9.]*\\).*/\\1/')
        ggplot2: \$(R -e "library(ggplot2); sessionInfo()" | grep ggplot2 | awk '{print \$2}')
        dplyr: \$(R -e "library(dplyr); sessionInfo()" | grep dplyr | awk '{print \$2}')
        cowplot: \$(R -e "library(cowplot); sessionInfo()" | grep cowplot | awk '{print \$2}')
        genomicranges: \$(R -e "library(GenomicRanges); sessionInfo()" | grep GenomicRanges | awk '{print \$2}')
    END_VERSIONS
    """
} 