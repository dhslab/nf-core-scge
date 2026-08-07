process  MAKE_CIRCOS_PLOT {
    tag "$meta.id"
    label 'process_low'
    container "staphb/circos"

    input:
    tuple val(meta), path(circos_input)

    output:
    tuple val(meta), path("${meta.id}.transgene_insertions_circos.png"), emit: plot

    // Invoke circos from PATH, not an absolute versioned path. The container tag above floats,
    // so a hardcoded /circos-<version>/bin/circos breaks silently the first time upstream
    // publishes a new release -- 0.69-9 -> 0.69-10 did exactly that (exit 127).
    script:
    """
    if [ -s "$circos_input" ]; then
        cp ${projectDir}/bin/circos.conf . && \\
        sed -i "s|__FILE_PLACEHOLDER__|$circos_input|" circos.conf && \\
        circos -conf circos.conf
        mv *.png ${meta.id}.transgene_insertions_circos.png
    else
        touch "${meta.id}.transgene_insertions_circos.png"
    fi

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.transgene_insertions_circos.png"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}