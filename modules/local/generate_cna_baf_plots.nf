process GENERATE_CNA_BAF_PLOTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/dhslab/docker-clerbase:250719'

    publishDir "${params.outdir}/pipeline_info/cna_baf_plots/${meta.id}", mode: 'copy', pattern: '*.png'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*")

    output:
    tuple val(meta), path("*cna_plot.png"), emit: cna_plot
    tuple val(meta), path("*baf_plot.png"), emit: baf_plot
    path "versions.yml", emit: versions

    script:
    def inputs = [
        meta.id,
        dragen_files.find{ it ==~ /.*\.(baf.bedgraph.gz)$/ } ?: "",
        dragen_files.find{ it ==~ /.*\.(tn.tsv.gz)$/ } ?: ""
    ].join(' ').trim()

    """
    generate_cna_baf_plots.R ${inputs}

    cat <<-'END_VERSIONS' > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed 's/.*version \\([0-9.]*\\).*/\\1/')
        ggplot2: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('ggplot2')))")
        dplyr: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('dplyr')))")
        cowplot: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('cowplot')))")
        genomicranges: \$(R --vanilla --quiet -e "cat(as.character(packageVersion('GenomicRanges')))")
    END_VERSIONS
    """

    stub:
    def inputs = [
        ${meta.id},
        dragen_files.find{ it ==~ /.*\.(baf.bedgraph.gz)$/ } ?: "",
        dragen_files.find{ it ==~ /.*\.(tn.tsv.gz)$/ } ?: ""
    ].join(' ').trim()

    """

    touch ${meta.id}.cna_plot.png
    touch ${meta.id}.baf_plot.png


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