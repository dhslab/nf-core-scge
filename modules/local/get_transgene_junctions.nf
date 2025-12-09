process GET_TRANSGENE_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(meta), path(dragen_files)
    path(fasta)

    output:
    tuple val(meta), path("${meta.id}.transgene_out.tsv"), emit: transgene_file
    path "versions.yml", emit: versions

    script:
    def input = [
        fasta.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--reference $it" } ?: "",
        params.transgene_match_coordinates ? "-x ${params.transgene_match_coordinates}" : "",
        params.transgene_name ? "${params.transgene_name}" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "$it" } ?: "",
    ].join(' ').trim()

    """

    getTransgeneJunctions.py ${input} > ${meta.id}.transgene_out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def input = [
        fasta ? "--reference ${fasta}" : "",
        params.transgene_match_coordinates ? "-x ${params.transgene_match_coordinates}" : "",
        transgene_name ? "${transgene_name}" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "$it" } ?: "",
    ].join(' ').trim()

    """
    touch ${meta.id}.transgene_out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}