process ANNOTATE_OFFTARGETS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep:release_105" // Use the same VEP container

    input:
    tuple val(meta), path(indels_file)

    output:
    tuple val(meta), path("${meta.id}.indels.annotated.tsv"), emit: annotated_indels
    path "versions.yml", emit: versions

    script:
    def pattern_chr = /'(#chrom|chr|chromosome)'/
    def pattern_start = /'(start|pos)'/
    def pattern_end = /'(end)'/
    """
    # Check if indels file is empty or has only a header
    if [ ! -s "${indels_file}" ] || [ \$(tail -n +2 "${indels_file}" | wc -l) -eq 0 ]; then
        if [ -s "${indels_file}" ]; then
            head -n 1 "${indels_file}" | tr -d '\\r\\n' > "${meta.id}.indels.annotated.tsv"
            echo -e '\\tAnnotation' >> "${meta.id}.indels.annotated.tsv"
        else
            touch "${meta.id}.indels.annotated.tsv"
        fi
        echo -e '"${task.process}":\\n  vep: N/A' > versions.yml
        exit 0
    fi

    HEADER=\$(head -n 1 ${indels_file} | sed 's/\\r\$//')
    CHRO_COL=\$(echo "\$HEADER" | tr '\\t' '\\n' | grep -n -i -E ${pattern_chr} | head -n 1 | cut -d: -f1)
    START_COL=\$(echo "\$HEADER" | tr '\\t' '\\n' | grep -n -i -E ${pattern_start} | head -n 1 | cut -d: -f1)
    END_COL=\$(echo "\$HEADER" | tr '\\t' '\\n' | grep -n -i -E ${pattern_end} | head -n 1 | cut -d: -f1)

    if [ -z "\$CHRO_COL" ] || [ -z "\$START_COL" ] || [ -z "\$END_COL" ]; then
        echo "Error: Could not find all required coordinate columns (chr, start, end) in ${indels_file}" >&2
        echo "Header was: \$HEADER" >&2
        exit 1
    fi

    tail -n +2 "${indels_file}" | awk -v c=\$CHRO_COL -v s=\$START_COL -v e=\$END_COL -F'\\t' 'BEGIN {OFS="\\t"} {print \$c, \$s, \$e, NR, "+"}' > vep_input.tsv

    VEP_OUTPUT="vep_output.tsv"
    /opt/vep/src/ensembl-vep/vep \\
        --offline \\
        --cache \\
        --dir ${params.vep_cache} \\
        --fasta ${params.fasta} \\
        --symbol \\
        --per_gene \\
        --tab \\
        --fields Location,SYMBOL \\
        --format region \\
        -i vep_input.tsv \\
        -o \${VEP_OUTPUT}

    if [ -s "\${VEP_OUTPUT}" ]; then
        grep -v '^##' \${VEP_OUTPUT} | cut -f 2 | sed 's/-/intergenic/g' > vep_symbols.txt
    else
        touch vep_symbols.txt
    fi

    NUM_INDELS=\$(tail -n +2 "${indels_file}" | wc -l)
    NUM_VEP_RESULTS=\$(cat vep_symbols.txt | wc -l)
    for i in \$(seq \$NUM_VEP_RESULTS \$((NUM_INDELS - 1)) ); do
        echo "intergenic" >> vep_symbols.txt
    done

    (echo "Annotation"; cat vep_symbols.txt) > annotations.txt

    paste "${indels_file}" annotations.txt > "${meta.id}.indels.annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\\s*//g')
    END_VERSIONS
    """
}


