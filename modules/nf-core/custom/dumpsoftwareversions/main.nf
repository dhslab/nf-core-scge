process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'
    container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
