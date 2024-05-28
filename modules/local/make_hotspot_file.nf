process MAKE_HOTSPOT_FILE {
    tag "${info[0].id}"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(info), path(csv_file), path (bed_file)

    output:
    tuple val(info), path ("*.vcf")  , emit: hotspot_vcf
    path "versions.yml" , emit: versions

    script:
    if( bed_file.name == 'NO_FILE.bed' && csv_file.name == 'NO_FILE.csv' ) {
        """
        touch NO_FILE.vcf
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
    else {
        def bed = bed_file.name != 'NO_FILE.bed' ? "--bed $bed_file" : ''
        def csv = csv_file.name != 'NO_FILE.csv' ? "--csv $csv_file" : ''
        """
        make_hotspot_file.py $bed $csv --window $params.window

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}
