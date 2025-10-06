process VEP_TO_TSV {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), val(type), path(input)

    output:
    tuple val(meta), path("*.tsv"), emit: vep_tsv
    path "versions.yml", emit: versions

    script:
    def args =
    (type == "vcf" ? "-i 1 -v" :
    type == "cnv" ? "-i 0 -s" :
    type == "sv" ? "-i 1 -s" : "")
    def output = input.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')
    """
    vep2table.py $args $input -o $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}