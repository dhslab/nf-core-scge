process CONVERT_MGI_SAMPLEMAP {
//    tag "${task.ext.prefix.id}"
    label 'process_low'

    container 'docker.io/gregorysprenger/pandas-excel:v2.2.2'

    input:
    tuple val(meta), path(samplemap)

    output:
    tuple val(meta), path("${meta.id}.fastq_list.csv")   , emit: fastq_list
    path("versions.yml")                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    convert_mgi_samplemap.py ${samplemap} ${meta.id}.fastq_list.csv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    convert_mgi_samplemap.py ${samplemap} ${meta.id}.fastq_list.csv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
