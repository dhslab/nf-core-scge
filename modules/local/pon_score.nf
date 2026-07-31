// PON_SCORE — score a run's own UNEDITED control against that run's own target file.
//
// This is what makes rule 4 of the review filter transfer to a brand-new guide RNA. A panel of
// normals is a list of genomic positions, so a PoN built for one guide shares almost no
// coordinates with a different guide's Cas-OFFinder panel and silently filters nothing.
// Scoring the control at the SAME target file the samples were called against produces a PoN
// with 100% site coverage by construction, for any guide, with no extra sequencing: the
// unedited control CRAM is already required by the pipeline.
//
// The query here is the control itself, so every call it makes is a false positive by
// construction -- germline or a recurrent alignment artifact at a homology-predicted site.
//
//   * -x 999999999 disables control subtraction. We want the raw evidence, not a filtered call.
//   * NO -u flag: that log is ~0.5-1 TB/sample and is what filled scratch2 during the CART run.
process PON_SCORE {
    tag "$meta.id"
    label 'process_highmem'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*"), path(hotspot_file)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.pon_scored.tsv"), emit: scored
    path "versions.yml", emit: versions

    script:
    def fasta = reference.find{ it ==~ /.*\.(fasta|fa)$/ }
    // Same control-CRAM selection as GET_INDELS (shortest filename), so the PoN is built from
    // exactly the sample GET_INDELS treats as the matched normal.
    def cram  = dragen_files.findAll{ it ==~ /.*\.(cram)$/ }.min{ it.toString().length() }
    """
    find_edited_reads.py \\
        --fasta ${fasta} \\
        --edited-bam ${cram} \\
        --control-bam ${cram} \\
        --target-file ${hotspot_file} \\
        -x 999999999 \\
        -o ${meta.id}.pon_scored.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pon_scored.tsv
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


// BUILD_PON — collapse the scored controls into the blacklist review_filter.py consumes.
process BUILD_PON {
    tag "cohort"
    label 'process_single'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path scored_tsvs
    path existing_pon

    output:
    path "offtarget_pon.tsv", emit: pon
    path "versions.yml",      emit: versions

    script:
    // Carrying a previous PoN forward accumulates recurrent bad regions across runs. Sites are
    // keyed by position, so this only compounds where panels genuinely overlap -- which is
    // precisely the repetitive sequence worth remembering.
    def merge = existing_pon.name != 'NO_FILE' ? "--merge ${existing_pon}" : ''
    """
    build_offtarget_pon.py ${scored_tsvs} ${merge} \\
        --min-reads ${params.pon_min_reads} \\
        --min-donors ${params.pon_min_donors} \\
        --min-vaf ${params.pon_min_vaf} \\
        -o offtarget_pon.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // versions.yml must carry real content even in a stub: nf-core's processVersionsFromYAML
    // does yaml.load(f).collectEntries{...}, and an empty file loads as null -> NPE.
    """
    touch offtarget_pon.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
