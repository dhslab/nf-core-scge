process GET_INDELS {
    tag "$meta.id"
    label 'process_highmem'
    label 'final_output'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*"), path(hotspot_file)
    path(crispr_model)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.offtarget_analysis.tsv"), emit: indels_file
    tuple val(meta), path("${meta.id}.offtarget_edits.vcf"), emit: indels_vcf
    tuple val(meta), path("${meta.id}.tagged.bam*"), optional: true, emit: tagged_bam
    //tuple val(meta), path("${meta.id}.ml_results.txt"), emit: ml_results
    //tuple val(meta), path("${meta.id}.fp_filtered.txt"), emit: fp_log
    path "versions.yml",    emit: versions

    script:
    def inputs = [
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: "",
//        crispr_model ? "--enable-crispr-prediction --crispr-model ${crispr_model} --crispr-threshold 0.7" : "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.max{ it.toString().length() }?.with{ "--edited-bam $it" } ?: "",
        dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }?.with{ "--control-bam $it" } ?: "",
        hotspot_file ? "--target-file ${hotspot_file}" : ""
    ].join(' ').trim()

    // Optional review aid: a window-restricted BAM whose reads carry an XC tag naming
    // the per-read call, for colouring the pileup in IGV. Off unless explicitly asked.
    def tagged_bam = params.offtarget_tagged_bam ? "--tagged-bam-out ${meta.id}.tagged.bam" : ""

    """
    find_edited_reads.py ${inputs} ${tagged_bam} -u ${meta.id}.unevaluable_reads.txt --vcf-out ${meta.id}.offtarget_edits.vcf -o ${meta.id}.offtarget_analysis.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def tagged_bam = params.offtarget_tagged_bam ? "touch ${meta.id}.tagged.bam ${meta.id}.tagged.bam.bai" : ""
    """
    touch ${meta.id}.offtarget_analysis.tsv
    touch ${meta.id}.offtarget_edits.vcf
    ${tagged_bam}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}