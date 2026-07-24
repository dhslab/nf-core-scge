// CRISPRME — enumerate a guide's off-target sites (mismatches + bulges, optionally
// variant-aware) against a PREBUILT CRISPRme index. Emits the CRISPRme targets TSV that
// bin/combine_offtarget_results.py parses. Off by default (params.run_crisprme); the index
// is a heavy one-time offline asset built by bin/build_crisprme_index.sh — this module
// only consumes it (params.crisprme_index_dir), never builds it.
//
// NOTE: the exact complete-search argument layout must be validated against the built
// index before the first real run (this is the fast-follow to the Cas-OFFinder path).
process CRISPRME {
    tag "${meta.id}"
    label 'process_high'
    container "docker.io/pinellolab/crisprme:latest"

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
    # PAM file + guide file expected by crisprme complete-search
    printf '%s%s %d\\n' "\$(printf 'N%.0s' \$(seq 1 ${spacer.length()}))" "${pam}" ${pam.length()} > pam.txt
    printf '%s%s\\n' "${spacer}" "\$(printf 'N%.0s' \$(seq 1 ${pam.length()}))" > guide.txt

    crisprme.py complete-search \\
        --genome ${index_dir}/Genome \\
        --thread ${task.cpus} \\
        --bmax ${bul} \\
        --mm ${mm} \\
        --pam pam.txt \\
        --guide guide.txt \\
        --output ${meta.id}_crisprme

    # complete-search writes its best-hits table into the output dir; expose the canonical name
    cp \$(find ${meta.id}_crisprme -name '*.best.txt' | head -1) ${meta.id}.crisprme.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crisprme: \$(crisprme.py --version 2>&1 | tail -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'crisprme_header\\n' > ${meta.id}.crisprme.tsv
    touch versions.yml
    """
}
