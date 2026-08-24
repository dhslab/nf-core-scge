process MAKE_HOTSPOT_VCF {
    tag "${id}"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(id), path(editing_targets)
    path(bed_file)
    path(reference)

    output:
    tuple val(id), path("${id}.hotspots.vcf")         , emit: hotspot_vcf
    path "versions.yml"                               , emit: versions

    script:
    def args = [
        bed_file.find{ it ==~ /.*\.(bed)$/ }?.with{ "--bed $it" }                                  ?: "",
        editing_targets.find{ it ==~ /.*\.(vcf|vcf.gz)$/ }?.with{ "--targets $it" }                       ?: "",
        params.hotspot_window_size ? "--window ${params.hotspot_window_size}"                       : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" }                          ?: ""
    ].join(' ').trim()
    """
    make_hotspot_vcf.py ${args} --outfile ${id}.hotspots.vcf

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${id}.hotspots.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
