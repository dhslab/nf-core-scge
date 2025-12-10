process VEP_TO_TSV {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-baseimage:latest"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(input)
    val(type)

    output:
    tuple val(meta), path("*.tsv"), emit: vep_tsv
    path "versions.yml", emit: versions

    script:
    def args = [
     type.toLowerCase() == "vcf" ? "-i 0 1 -v" : "",
     type.toLowerCase() == "sv" ? "-i 1 -s" : "",
     type.toLowerCase() == "cnv" ? "-i 0 -s" : ""
    ].join(' ').trim()
    def vcf = input.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def output = vcf.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')

    """
    export PATH=/usr/local/bin:\$PATH
    vep2table.py $args $vcf -o $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def args = [
     type.toLowerCase() == "vcf" ? "-i 0 1 -v" : "",
     type.toLowerCase() == "sv" ? "-i 1 -s" : "",
     type.toLowerCase() == "cnv" ? "-i 0 -s" : ""
    ].join(' ').trim()
    def vcf = input.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def output = vcf.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')

    """
    touch "${output}"

    cat <<-END_CMDS > "${task.process}_cmds.txt"
    vep2table.py $args $vcf -o $output
    END_CMDS

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}