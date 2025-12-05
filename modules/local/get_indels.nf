process GET_INDELS {
    tag "$meta.id"
    label 'process_high'
    label 'final_output'
    container "ghcr.io/dhslab/docker-scge:latest"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*"), path(hotspot_file)
    path(crispr_model)

    output:
    tuple val(meta), path("${meta.id}.indels.txt"), emit: indels_file
    tuple val(meta), path("${meta.id}.ml_results.txt"), emit: ml_results
    tuple val(meta), path("${meta.id}.fp_filtered.txt"), emit: fp_log
    path "versions.yml",    emit: versions

    script:
    def inputs = [
        crispr_model ? "--enable-crispr-prediction --crispr-model ${crispr_model} --crispr-threshold 0.7" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "--edited-bam $it" } ?: "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }?.with{ "--control-bam $it" } ?: "",
        hotspot_file ? "--target-file ${hotspot_file}" : ""
    ].join(' ').trim()

    """
    extract_variant_reads_ML.py ${inputs} --filter-off-target-fp \\
        --fp-log ${meta.id}.fp_filtered.txt -v -o ${meta.id}.indels.txt

    # Extract ML results into a separate file, preserving the header
    cut -f 17-19 ${meta.id}.indels.txt > ${meta.id}.ml_results.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def inputs = [
        crispr_model ? "--enable-crispr-prediction --crispr-model ${crispr_model} --crispr-threshold 0.7" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "--edited-bam $it" } ?: "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }?.with{ "--control-bam $it" } ?: "",
        hotspot_file ? "--target-file ${hotspot_file}" : ""
    ].join(' ').trim()
    
    """
    touch ${meta.id}.fp_filtered.txt
    touch ${meta.id}.indels.txt
    touch ${meta.id}.ml_results.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}