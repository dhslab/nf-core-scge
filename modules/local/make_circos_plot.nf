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
    if [ -s "$circos_input" ]; then
        cp ${projectDir}/bin/circos.conf . && \\
        sed -i "s|__FILE_PLACEHOLDER__|$circos_input|" circos.conf && \\
        /circos-0.69-9/bin/circos -conf circos.conf
    else
        touch "${meta.id}_circos.png"
    fi

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}