// CASOFFINDER — enumerate a guide's genome-wide off-target sites with cas-offinder-bulge
// (mismatches + DNA/RNA bulges). Emits the native bulge-format table
// (Bulge type, crRNA, DNA, Chromosome, Position, Direction, Mismatches, Bulge Size) that
// bin/combine_offtarget_results.py parses. CPU device ('C'), so no GPU is required.
process CASOFFINDER {
    tag "${meta.id}"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-casoffinder-bulge:latest"

    input:
    tuple val(meta), val(spacer), val(pam)
    path twobit

    output:
    tuple val(meta), path("${meta.id}.casoffinder.txt"), emit: hits
    path "versions.yml",                                 emit: versions

    script:
    def mm  = params.offtarget_mismatches
    def bul = params.offtarget_bulges
    """
    # cas-offinder-bulge input: genome, then "<pattern> <DNA bulge> <RNA bulge>", then
    # "<query> <mismatches>". Pattern masks the protospacer with N and appends the PAM;
    # the query is the spacer with N placeholders for the PAM bases.
    SPACER=${spacer}
    PAM=${pam}
    PATTERN=\$(printf 'N%.0s' \$(seq 1 \${#SPACER}))\${PAM}
    QUERY=\${SPACER}\$(printf 'N%.0s' \$(seq 1 \${#PAM}))

    printf '%s\\n%s %s %s\\n%s %s\\n' \\
        "\$(readlink -f ${twobit})" "\$PATTERN" "${bul}" "${bul}" "\$QUERY" "${mm}" > casoffinder_input.txt

    cas-offinder-bulge casoffinder_input.txt C ${meta.id}.casoffinder.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cas-offinder: \$(cas-offinder 2>&1 | grep -oE 'Cas-OFFinder v[0-9.a-z]+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'Bulge type\\tcrRNA\\tDNA\\tChromosome\\tPosition\\tDirection\\tMismatches\\tBulge Size\\n' > ${meta.id}.casoffinder.txt
    touch versions.yml
    """
}
