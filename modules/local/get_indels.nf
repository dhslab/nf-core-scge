process GET_INDELS {
    tag "$meta.id"
    label 'process_high'
    label 'final_output'
    container "ghcr.io/dhslab/docker-scge:latest"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(files), path(hotspot_file, stageAs: 'hotspots.csv')
    path crispr_model

    output:
    tuple val(meta), path("${meta.id}.indels.txt"), emit: indels_file
    tuple val(meta), path("${meta.id}.ml_results.txt"), emit: ml_results
    path "versions.yml",    emit: versions

    script:
    """
    set -euo pipefail

    # Determine tumor and control CRAMs from staged files
    TUMOR_CRAM=""
    CONTROL_CRAM=""

    if ls dragen/*_tumor.cram >/dev/null 2>&1; then
        TUMOR_CRAM=\$(ls dragen/*_tumor.cram | head -n1)
    fi

    if ls dragen/*_normal.cram >/dev/null 2>&1; then
        CONTROL_CRAM=\$(ls dragen/*_normal.cram | head -n1)
    fi

    # Fallbacks: try explicit control CRAM staged from "unedited" directory
    if [ -z "\${CONTROL_CRAM}" ]; then
        if ls dragen/*unedited*.cram >/dev/null 2>&1; then
            CONTROL_CRAM=\$(ls dragen/*unedited*.cram | head -n1)
        fi
    fi

    # Fallbacks: try by meta id or any CRAMs present
    if [ -z "\${TUMOR_CRAM}" ] && ls dragen/${meta.id}.cram >/dev/null 2>&1; then
        TUMOR_CRAM=dragen/${meta.id}.cram
    fi
    if [ -z "\${CONTROL_CRAM}" ] && ls dragen/${meta.id}.cram >/dev/null 2>&1; then
        CONTROL_CRAM=dragen/${meta.id}.cram
    fi

    if [ -z "\${TUMOR_CRAM}" ] || [ ! -s "\${TUMOR_CRAM}" ]; then
        # pick the first CRAM as tumor
        TUMOR_CRAM=\$(ls dragen/*.cram | head -n1)
    fi
    if [ -z "\${CONTROL_CRAM}" ] || [ ! -s "\${CONTROL_CRAM}" ]; then
        # pick the next CRAM as control if available; else use tumor
        CONTROL_CRAM=\$(ls dragen/*.cram | sed -n '2p')
        CONTROL_CRAM=\${CONTROL_CRAM:-\$TUMOR_CRAM}
    fi

    # Ensure hotspot file is present as hotspots.csv; if missing, create an empty template with expected columns
    BED_FILE="hotspots.csv"
    if [ ! -s "\${BED_FILE}" ]; then
        echo "[GET_INDELS] hotspots.csv not found; creating empty template" >&2
        echo "Chromosome,Start,End,On_target,Source,DNA_Sequence,PAM,Strand,Mismatch,Bulge_Type,Bulge_Size" > \${BED_FILE}
    fi

    extract_variant_reads_ML.py \
        --target-file \${BED_FILE} \
        --edited-bam \${TUMOR_CRAM} \
        --control-bam \${CONTROL_CRAM} \
        --enable-crispr-prediction \
        --crispr-model ${crispr_model} \
        --crispr-threshold 0.7 \
        --filter-off-target-fp \
        -v \
        -o ${meta.id}.indels.txt

    # Extract ML results into a separate file, preserving the header
    cut -f 17-19 ${meta.id}.indels.txt > ${meta.id}.ml_results.txt

    cat <<-'END_VERSIONS' > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}