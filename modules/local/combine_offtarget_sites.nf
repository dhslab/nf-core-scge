// COMBINE_OFFTARGET_SITES — merge one guide's per-tool off-target predictions
// (Cas-OFFinder bulge output, CRISPRme, optionally IDT) into the canonical
// `<guide>.targets.csv` schema the OFFTARGET arm consumes as a target_file.
// Any source may be omitted by passing assets/NO_FILE; at least one must be real.
process COMBINE_OFFTARGET_SITES {
    tag "${meta.id}"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    tuple val(meta), path(casoffinder), path(crisprme), path(idt)

    output:
    tuple val(meta), path("${meta.id}.targets.csv"), emit: sites
    path "versions.yml",                             emit: versions

    script:
    // Optional slots are passed as distinct NO_* sentinel files (distinct names avoid a
    // staging collision); a real tool output never starts with 'NO_'.
    def real = { f -> f && !f.name.startsWith('NO_') }
    def cas_arg = real(casoffinder) ? "--casoffinder ${casoffinder}" : ''
    def cme_arg = real(crisprme)    ? "--crisprme ${crisprme}"       : ''
    def idt_arg = real(idt)         ? "--idt ${idt}"                 : ''
    """
    python ${projectDir}/bin/combine_offtarget_results.py \\
        ${cas_arg} \\
        ${cme_arg} \\
        ${idt_arg} \\
        --pam ${params.offtarget_pam} \\
        -o ${meta.id}.targets.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.targets.csv versions.yml"
}
