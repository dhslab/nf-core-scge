process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'
    container "ghcr.io/dhslab/docker-cleutils:240129"

    input:
    path samplesheet
    val data_path

    output:
    path 'samplesheet.valid.csv'    , emit: csv
    path "versions.yml"             , emit: versions

    script:
    def mgi = params.mgi == true ? "-m " : ""
    def dir = data_path != "" ? "-d ${data_path} " : ""
    """
    samplesheet_check.py ${mgi}${dir}$samplesheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def mgi = params.mgi == true ? "-m " : ""
    def dir = data_path != "" ? "-d ${data_path} " : ""
    """
    samplesheet_check.py ${mgi}${dir}$samplesheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
