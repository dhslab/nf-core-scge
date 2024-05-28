process GET_TRANSGENE_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.transgene_out.tsv")
    path "versions.yml"                                     , emit: versions

    script:
    """
    getTransgeneJunctions.v2.py -x 3130,5930 PLVM_CD19_CARv4_cd34 ${meta.id}_tumor.cram > ${meta.id}.transgene_out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}