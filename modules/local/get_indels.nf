process GET_INDELS {
    tag "$meta.id"
    label 'process_highmem'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*"), path(hotspot_file)
    path(crispr_model)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.offtarget_analysis.tsv"), emit: indels_file
    tuple val(meta), path("${meta.id}.ml_results.txt"), emit: ml_results
    tuple val(meta), path("${meta.id}.fp_filtered.txt"), emit: fp_log
    path "versions.yml",    emit: versions

    script:
    def inputs = [
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
        crispr_model ? "--enable-crispr-prediction --crispr-model ${crispr_model} --crispr-threshold 0.7" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "--edited-bam $it" } ?: "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }?.with{ "--control-bam $it" } ?: "",
        hotspot_file ? "--target-file ${hotspot_file}" : ""
    ].join(' ').trim()

    """    
    find_edited_reads.py ${inputs} -u ${meta.id}.unevaluable_reads.txt --filter-off-target-fp \\
        --fp-log ${meta.id}.fp_filtered.txt -v -o ${meta.id}.offtarget_analysis.tsv

    # Extract ML results into a separate file, preserving the header
    cut -f 17-19 ${meta.id}.offtarget_analysis.tsv > ${meta.id}.ml_results.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def inputs = [
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
        crispr_model ? "--enable-crispr-prediction --crispr-model ${crispr_model} --crispr-threshold 0.7" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "--edited-bam $it" } ?: "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }?.with{ "--control-bam $it" } ?: "",
        hotspot_file ? "--target-file ${hotspot_file}" : ""
    ].join(' ').trim()
    
    """
    touch ${meta.id}.offtarget_indels.tsv
    touch ${meta.id}.fp_filtered.txt
    touch ${meta.id}.ml_results.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}