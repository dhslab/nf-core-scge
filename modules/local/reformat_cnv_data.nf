process REFORMAT_CNV_DATA { 
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils:240129"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.cnv_visualization.tsv")    , emit: out_file
    path("versions.yml")                                         , emit: versions

    script:
    """
    set -euo pipefail
    if [ -s "${meta.id}.tumor.baf.bedgraph.gz" ] && [ -s "${meta.id}.tn.tsv.gz" ]; then
        cnv_visualization.py -o ${meta.id}.cnv_visualization.tsv "${meta.id}.tumor.baf.bedgraph.gz" "${meta.id}.tn.tsv.gz"
    else
        echo -e "sample\tmessage" > ${meta.id}.cnv_visualization.tsv
        echo -e "${meta.id}\tCNA/BAF inputs unavailable" >> ${meta.id}.cnv_visualization.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}