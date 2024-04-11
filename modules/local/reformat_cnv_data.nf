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
    cnv_visualization.py -o ${meta.id}.cnv_visualization.tsv ${meta.id}.tumor.baf.bedgraph.gz ${meta.id}.tn.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cnv_visualization.py --version)
    END_VERSIONS
    """

}