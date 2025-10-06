process MAKE_HOTSPOT_VCF {
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    path(bed_file)
    path(reference)

    output:
    path("hotspots.vcf")  , emit: hotspot_vcf
    path "versions.yml"   , emit: versions

    script:
    def args = [
        bed_file.find{ it ==~ /.*\.(bed)$/ }?.with{ "--bed $it" }                                  ?: "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" }                          ?: ""
    ].join(' ').trim()
    """
    make_hotspot_vcf.py ${args} --outfile hotspots.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        make_hotspot_vcf: \$(make_hotspot_vcf.py --version | awk '{print \$2}')
    END_VERSIONS
    """

}
