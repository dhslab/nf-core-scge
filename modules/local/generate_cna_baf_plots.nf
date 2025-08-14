process GENERATE_CNA_BAF_PLOTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-clerbase:250719'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*")

    output:
    tuple val(meta), path("cna_plot.png"), emit: cna_plot
    tuple val(meta), path("baf_plot.png"), emit: baf_plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    baf_path="ls dragen_files/*.baf.bedgraph.gz"
    cnv_path="ls dragen_files/*.tn.tsv.gz"

    if [ -n "\$baf_path" ] && [ -n "\$cnv_path" ] && [ -e "\$baf_path" ] && [ -e "\$cnv_path" ]; then
        generate_cna_baf_plots.R \
            ${meta.id} \
            \$baf_path \
            \$cnv_path
    else
        echo "[GENERATE_CNA_BAF_PLOTS] Missing BAF/CNV inputs for ${meta.id}; creating placeholder plots" >&2
        R --vanilla <<'RSCRIPT'
        png("cna_plot.png", width=1800, height=900, res=150)
        par(mar=c(0,0,0,0))
        plot.new(); text(0.5, 0.5, "CNA plot unavailable", cex=2)
        dev.off()
        png("baf_plot.png", width=1800, height=900, res=150)
        par(mar=c(0,0,0,0))
        plot.new(); text(0.5, 0.5, "BAF plot unavailable", cex=2)
        dev.off()
RSCRIPT
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed 's/.*version \\([0-9.]*\\).*/\\1/')
        ggplot2: \$(R -e "suppressMessages(library(ggplot2)); sessionInfo()" | grep ggplot2 | awk '{print \$2}')
        dplyr: \$(R -e "suppressMessages(library(dplyr)); sessionInfo()" | grep dplyr | awk '{print \$2}')
        cowplot: \$(R -e "suppressMessages(library(cowplot)); sessionInfo()" | grep cowplot | awk '{print \$2}')
        genomicranges: \$(R -e "suppressMessages(library(GenomicRanges)); sessionInfo()" | grep GenomicRanges | awk '{print \$2}')
    END_VERSIONS
    """
} 