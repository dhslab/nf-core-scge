// WGS_WORKLIST — genome-wide, homology-free off-target discovery + shape ranking.
// Wraps bin/worklist_from_vcf.py: reads each DRAGEN somatic VCF + tumor CRAM listed in
// the cram map, applies the depth-augmented shape model, cross-checks the matched normal.
// Cohort-level (the script loops all samples in the map). CRAMs/VCFs are read by absolute
// path from bind-mounted storage (as elsewhere in this pipeline), so they are not staged.
process WGS_WORKLIST {
    tag "cohort"
    label 'process_medium'
    container "ghcr.io/dhslab/docker-scge:latest"

    input:
    path cram_map
    path model
    val  reference        // absolute FASTA path on mounted storage

    output:
    path "wgs_offtarget_worklist_genomewide.csv", emit: worklist
    path "snapshots/*.png",  optional: true,       emit: snapshots
    path "versions.yml",                           emit: versions

    script:
    def homology = params.offtarget_homology_table ? "--homology-table '${params.offtarget_homology_table}'" : ""
    def snaps    = params.offtarget_snapshots ? "--snapshot-dir snapshots" : ""
    """
    worklist_from_vcf.py \\
        --cram-list ${cram_map} \\
        --ref ${reference} \\
        --model ${model} \\
        --min-af ${params.offtarget_min_af} \\
        --min-span ${params.offtarget_min_span} \\
        --top ${params.offtarget_top} \\
        ${homology} ${snaps} \\
        --out wgs_offtarget_worklist_genomewide.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "mkdir -p snapshots; touch wgs_offtarget_worklist_genomewide.csv versions.yml"
}
