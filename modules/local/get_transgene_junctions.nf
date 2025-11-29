process GET_TRANSGENE_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-cleutils"

    input:
    tuple val(meta), path(dragen_dir), val(transgene_name), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.transgene_out.tsv"), emit: transgene_file
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    TUMOR_CRAM=\$(ls ${dragen_dir}/*_tumor.cram 2>/dev/null | head -n1 || echo "")
    if [ -z "\$TUMOR_CRAM" ]; then
        TUMOR_CRAM=\$(ls ${dragen_dir}/*.cram 2>/dev/null | head -n1 || echo "")
    fi
    
    if [ -z "\$TUMOR_CRAM" ]; then
        echo "Error: No tumor CRAM file found." >&2
        exit 1
    fi

    getTransgeneJunctions.py -x 3130,5930 ${transgene_name} \$TUMOR_CRAM --reference ${fasta} > ${meta.id}.transgene_out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}