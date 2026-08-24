// CRISPRME — enumerate a guide's off-target sites (mismatches + DNA/RNA bulges) against a
// PREBUILT CRISPRme/CRISPRitz index. Emits the CRISPRme integrated-results TSV that
// bin/combine_offtarget_results.py parses. Off by default (params.run_crisprme); the index is a
// heavy one-time offline asset built by bin/build_crisprme_index.sh — this module only consumes it.
//
// The index dir (params.crisprme_index_dir) must contain, as built by build_crisprme_index.sh:
//   Genome/                              per-chromosome unzipped FASTAs
//   genome_library/<PAM>_<bMax>_Genome/  the TST index (bMax must equal bDNA+bRNA of the search)
//   <len>bp-<PAM>-<nuclease>.txt          the PAM file (its name encodes the nuclease)
// CRISPRme looks for genome_library relative to the CWD, so we symlink the staged index into the
// task work dir before running complete-search — that reuses the prebuilt index instead of
// rebuilding it (verified: prebuilt .bin files are left untouched).
process CRISPRME {
    tag "${meta.id}"
    label 'process_high'
    // Lab wrapper over pinellolab/crisprme that adds procps (`ps`) — Nextflow's task wrapper needs
    // it under -euo pipefail. See containers/docker-crisprme/Dockerfile. Build/push before enabling.
    container "ghcr.io/dhslab/docker-crisprme:latest"
    // CRISPRme writes scratch files (e.g. vuoto.txt) into its own read-only install dir under
    // Apptainer; an ephemeral writable overlay lets those succeed without persisting anything.
    containerOptions '--writable-tmpfs'

    input:
    tuple val(meta), val(spacer), val(pam)
    path index_dir

    output:
    tuple val(meta), path("${meta.id}.crisprme.tsv"), emit: hits
    path "versions.yml",                              emit: versions

    script:
    def mm  = params.offtarget_mismatches
    def bul = params.offtarget_bulges
    """
    # CRISPRme ships in a conda env its entrypoint activates; Nextflow bypasses the entrypoint,
    # so put its tools on PATH.
    export PATH=/opt/conda/bin:\$PATH

    # Expose the prebuilt index in the CWD so complete-search reuses it (it resolves
    # genome_library relative to the working directory) instead of rebuilding.
    ln -s ${index_dir}/Genome Genome
    ln -s ${index_dir}/genome_library genome_library
    PAMFILE=\$(basename \$(ls ${index_dir}/*bp-*-*.txt | head -1))
    cp ${index_dir}/\$PAMFILE .

    # guide file: spacer padded with N over the PAM positions (CRISPRme's crRNA format)
    printf '%s%s\\n' "${spacer}" "\$(printf 'N%.0s' \$(seq 1 ${pam.length()}))" > guide.txt

    crisprme.py complete-search \\
        --genome Genome \\
        --pam \$PAMFILE \\
        --guide guide.txt \\
        --mm ${mm} \\
        --bDNA ${bul} \\
        --bRNA ${bul} \\
        --output ${meta.id}_crisprme \\
        --thread ${task.cpus}

    # the parseable target table is the integrated-results TSV (schema consumed by the combiner)
    cp Results/${meta.id}_crisprme/*_integrated_results.tsv ${meta.id}.crisprme.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crisprme: \$(crisprme.py --version 2>&1 | tail -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'Spacer+PAM\\tChromosome\\tStart\\tStrand\\tx\\tDNA\\tx\\tPAM\\tMismatches\\tBulges\\tx\\tx\\tx\\tx\\tx\\tBulge_type\\n' > ${meta.id}.crisprme.tsv
    touch versions.yml
    """
}
