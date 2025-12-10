process ANNOTATE_OFFTARGETS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep_release113:250810" // Use the same VEP container

    input:
    tuple val(meta), path(targetfile)
    path(vep_cache)
    path(fasta)

    output:
    tuple val(meta), path("${meta.id}.indels.annotated.tsv"), emit: targetfile
    path "versions.yml", emit: versions

    script:
    // Use params directly for vep_cache to ensure absolute path is used
    // Handle both params.vepcache (absolute path) and staged vep_cache (relative path)
    def vep_cache_dir = params.vepcache 
        ? params.vepcache.toString().replaceAll(/\/$/, '')  // Remove trailing slash if present
        : (vep_cache ? vep_cache.toString() : "")
    
    // Find fasta file from the collection
    def fasta_file = fasta.find{ it ==~ /.*\.(fasta|fa)$/ }?.toString() ?: ""
    
    def vep_args = [
        targetfile                                          ? "-i ${targetfile}"    : "",
        vep_cache_dir                                       ? "--dir ${vep_cache_dir}"   : "",
        fasta_file                                          ? "--fasta ${fasta_file}" : ""
    ].findAll{ it }.join(' ')

    """
    # Extract the last column header from the target file (contains the keys for the Info column)
    # Assumes the target file has a header line and is tab-separated (or compatible)
    LAST_COL_HEADER=\$(head -n 1 ${targetfile} | awk -F'\\t' '{print \$NF}' | tr -d '\\r')

    # Use add_vep2targetfile.pl to format VEP output into a compatible TSV
    # The helper script now handles header creation and column reformatting using the extracted key string.
    /opt/vep/src/ensembl-vep/vep \\
            --offline \\
            --cache \\
            ${vep_args} \\
            --force --symbol --term SO --per_gene --fields "Location,SYMBOL,DISTANCE,INTRON,EXON" --numbers -o stdout \\
            | add_vep2targetfile.pl "\$LAST_COL_HEADER" > "${meta.id}.indels.annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """
}


