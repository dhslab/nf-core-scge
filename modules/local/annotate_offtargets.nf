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
    # Extract the ENTIRE header from the target file (used as keys for the Info column)
    # This header should be comma-separated and will be parsed by extract_variant_reads_ML.py
    LAST_COL_HEADER=\$(head -n 1 ${targetfile} | tr -d '\\r')

    # Detect delimiter (comma or tab)
    if head -n 1 ${targetfile} | grep -q ','; then
        DELIM=","
    else
        DELIM="\\t"
    fi

    # Convert target CSV to VEP's ensembl input format with proper allele representation
    # VEP ensembl format: chromosome start end allele strand [identifier]
    # For target site annotation, we use N/N as a dummy SNP (any nucleotide)
    # CSV format: Source,DNA_Sequence,PAM,Chromosome,Strand,Start,...
    awk -F"\$DELIM" 'NR==1 {
        for (i=1; i<=NF; i++) {
            gsub(/^[ \\t]+|[ \\t]+\$/, "", \$i)  # trim whitespace
            if (\$i ~ /^[Cc]hromosome\$|^[Cc]hrom\$/) chr_col=i
            if (\$i ~ /^[Ss]tart\$/) start_col=i
            if (\$i ~ /^[Ss]trand( [Dd]irection)?\$/) strand_col=i
            if (\$i ~ /^[Oo]n.?[Tt]arget\$/) ontarget_col=i
        }
        next
    }
    chr_col && start_col {
        chr=\$chr_col
        start=\$start_col
        end=start  # For point annotation, end = start
        strand = strand_col ? \$strand_col : "+"
        if (strand == "") strand = "+"
        ontarget = ontarget_col ? \$ontarget_col : "0"
        # Build identifier from all columns (comma-separated to match header)
        # This will be parsed by extract_variant_reads_ML.py as key-value pairs
        id=""
        for (i=1; i<=NF; i++) {
            if (i>1) id=id","
            id=id\$i
        }
        # VEP ensembl format: chr start end allele strand identifier
        # Use A/T as dummy SNP allele - VEP will accept this and annotate the position
        print chr, start, end, "A/T", strand, id
    }' ${targetfile} > vep_input_unsorted.txt

    # Sort by chromosome and position (VEP requires sorted input)
    sort -k1,1V -k2,2n vep_input_unsorted.txt > vep_input.txt

    # Check if vep_input.txt has data
    if [ ! -s vep_input.txt ]; then
        echo "Warning: No valid records found in target file for VEP annotation"
        # Create empty output with header
        echo -e "#chromosome\\tstart\\tend\\tLocation\\tAllele\\tGene\\tFeature\\tFeature_type\\tConsequence\\tcDNA_position\\tCDS_position\\tProtein_position\\tAmino_acids\\tCodons\\tExisting_variation\\tExtra\\t\$LAST_COL_HEADER" > "${meta.id}.indels.annotated.tsv"
    else
        # Run VEP on converted input
    /opt/vep/src/ensembl-vep/vep \\
            --offline \\
            --cache \\
                -i vep_input.txt \\
                --dir ${vep_cache_dir} \\
                --fasta ${fasta_file} \\
            --force --symbol --term SO --per_gene --fields "Location,SYMBOL,DISTANCE,INTRON,EXON" --numbers -o stdout \\
            | add_vep2targetfile.pl "\$LAST_COL_HEADER" > "${meta.id}.indels.annotated.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """
}


