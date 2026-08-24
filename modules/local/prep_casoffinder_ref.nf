// PREP_CASOFFINDER_REF — build the .2bit Cas-OFFinder searches against, from the
// pipeline reference FASTA. Runs once per pipeline (faToTwoBit is cheap/deterministic);
// skip it entirely by setting params.casoffinder_2bit to a prebuilt .2bit.
process PREP_CASOFFINDER_REF {
    tag "${fasta.baseName}"
    label 'process_low'
    container "ghcr.io/dhslab/docker-casoffinder-bulge:latest"

    input:
    path fasta

    output:
    path "${fasta.baseName}.2bit", emit: twobit
    path "versions.yml",           emit: versions

    script:
    """
    faToTwoBit ${fasta} ${fasta.baseName}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faToTwoBit: \$(faToTwoBit 2>&1 | grep -oE 'faToTwoBit v[0-9.]+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    "touch ${fasta.baseName}.2bit versions.yml"
}
