process  MAKE_CIRCOS_PLOT {
    tag "$meta.id"
    label 'process_low'
    container "staphb/circos"

    input:
    tuple val(meta), path(circos_input)

    output:
    tuple val(meta), path("*png"), emit: circos_plot

    script:
    """
    cp ${projectDir}/bin/circos.conf . && \\
    sed -i "s|__FILE_PLACEHOLDER__|$circos_input|" circos.conf && \\
    /circos-0.69-9/bin/circos -conf circos.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circos: \$(echo \$(circos -v 2>&1) | sed 's/circos.*v //; s/ .*\$//')
    END_VERSIONS
    """
}