process PARTITION_ALIGNMENT_FILE {
    tag "${meta.id}"
    label "process_medium"

    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
    tuple val(meta), path(alignment_file)
    path(reference)

    output:
    tuple val(meta), path("*.cram"), emit: cram_files
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cram_reference  = reference         ? "--reference ${reference.min{ it.toString().length() }}"    : ""
    def cram_file       = alignment_file    ? "${alignment_file.min{ it.toString().length() }}"           : ""

    """
    samtools split \\
        -f '%!.cram' \\
        --write-index \\
        -@ ${task.cpus} \\
        ${cram_reference} \\
        ${cram_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    """

    stub:
    """
    find \\
        "${projectDir}/assets/stub/alignment_files"/*.{5,6}.cram \\
        -type f \\
        -name "*.cram" \\
        | while read -r file; do
            cp "\$file" .
        done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
