process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), val(type), path(input_vcf)

    output:
    tuple val(meta), val(type), path("$annotate_vcf_output"), emit: annotated_vcf
    path "versions.yml",    emit: versions

    script:
    annotate_vcf_output = "${meta.id}" + 
    (type == "vcf" ? ".hard_filtered.annotated.vcf" : 
    type == "cnv" ? ".cnv.annotated.vcf" : 
    type == "sv" ? ".sv.annotated.vcf" : "")
    """
    /usr/bin/perl \\
   -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
   --format vcf \\
   --vcf --fasta ${params.fasta} \\
   --hgvs \\
   --symbol \\
   --term SO \\
   --flag_pick \\
   --custom ${params.assay_inputs.cytobands},cytobands,bed \\
   -o $annotate_vcf_output \\
   -i $input_vcf \\
   --offline \\
   --cache \\
   --max_af --dir ${params.vep_cache}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        vep: \$(/opt/vep/src/ensembl-vep/vep --help 2>&1 | grep "ensembl-vep" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//')
    END_VERSIONS
    """

}