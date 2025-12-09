process CONVERT_ALIGNMENT_FILE_TO_FASTQ {
    tag "${meta.id}"
    label "process_high"

    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
    tuple val(meta), path(alignment_file)
    path(reference)

    output:
    tuple val(meta), path("*R1*.fastq.gz"), path("*R2*.fastq.gz"), emit: fastqs
    path("versions.yml")                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cram_reference  = reference         ? "--reference ${reference.min{ it.toString().length() }}"    : ""
    def cram_file       = alignment_file    ? "${alignment_file.min{ it.toString().length() }}"           : ""
    def read_group      = cram_file         ? cram_file - '.cram'                                         : "${meta.id}"
    def lane            = cram_file         ? read_group.split('\\.')[-1]                                 : "?"

    """
    samtools collate \\
        -@ ${task.cpus} \\
        ${cram_reference} \\
        -u \\
        -O \\
        ${alignment_file} \\
        | samtools fastq \\
            -@ ${task.cpus} \\
            -1 "${meta.id}.${read_group}.R1.fastq.gz" \\
            -2 "${meta.id}.${read_group}.R2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def cram_reference  = reference         ? "--reference ${reference.min{ it.toString().length() }}"    : ""
    def cram_file       = alignment_file    ? "${alignment_file.min{ it.toString().length() }}"           : ""
    def read_group      = cram_file         ? cram_file - '.cram'                                         : "${meta.id}"
    def lane            = cram_file         ? read_group.split('\\.')[-1]                                 : "?"

    """
    find \\
        "${projectDir}/assets/stub/alignment_files/" \\
        -name "${meta.id}.${read_group}*.fastq.gz" \\
        | sort \\
        | head -n 2 \\
        | while read -r line; do
            cp \$line \$PWD
        done

    cat <<-END_CMDS > "${meta.id}_cmds.txt"
    samtools collate \\
        -@ ${task.cpus} \\
        ${cram_reference} \\
        -u \\
        -O \\
        ${alignment_file} \\
        | samtools fastq \\
            -@ ${task.cpus} \\
            -1 "${meta.id}.R1.fastq.gz" \\
            -2 "${meta.id}.R2.fastq.gz"
    END_CMDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
