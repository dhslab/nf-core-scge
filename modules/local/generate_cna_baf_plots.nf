process GENERATE_CNA_BAF_PLOTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-clerbase:250719'

    publishDir "${params.outdir}/pipeline_info/cna_baf_plots/${meta.id}", mode: 'copy', pattern: '*.png'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*cna_plot.png"), emit: cna_plot
    tuple val(meta), path("*baf_plot.png"), emit: baf_plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail
    baf_path=\$(find -L . -name "*baf.bedgraph.gz")
    cnv_path=\$(find -L . -name "*.tn.tsv.gz")
    if [ -n "\$baf_path" ] && [ -n "\$cnv_path" ] && [ -e "\$baf_path" ] && [ -e "\$cnv_path" ]; then
        generate_cna_baf_plots.R \\
            ${meta.id} \\
            \$baf_path \\
            \$cnv_path
    fi

    ls -lR

    cat <<-'END_VERSIONS' > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed 's/.*version \\([0-9.]*\\).*/\\1/')
        ggplot2: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('ggplot2')))")
        dplyr: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('dplyr')))")
        cowplot: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('cowplot')))")
        genomicranges: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('GenomicRanges')))")
    END_VERSIONS
    """
} 