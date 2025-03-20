// python3 vep2table.py -i 0 -s out.Eta_C33_2.cnv.vcf -o cnv.tsv
// python3 vep2table.py -i 1 -s out.Eta_C33_2.sv.vcf -o sv.tsv
// python3 vep2table.py -i 1 -v out.Eta_C33_2.vcf -o vcf.tsv
process VEP_TO_TSV {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), val(type), path(input)

    output:
    tuple val(meta), path("$output")
    path "versions.yml",    emit: versions

    script:
    def args =
    (type == "vcf" ? "-i 1 -v" : 
    type == "cnv" ? "-i 0 -s" : 
    type == "sv" ? "-i 1 -s" : "")
    output = input.getName().replaceAll(/\.vcf$/, ".tsv")
    """
    vep2table.py $args $input -o $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}