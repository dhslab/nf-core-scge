process ANNOTATE_CNV_VARIANTS {
    tag "${meta.id}"
    label "process_low"
    container "ghcr.io/dhslab/docker-vep_release113:250810"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_files)
    path(reference)
    path(vep_cache)
    path(cytobands)

    output:
    tuple val(meta), path("*.cnv.annotated.vcf.gz*"), emit: vcf
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def vcf = dragen_files.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def vep_gene_args = [
        vep_cache                                 ? "--dir ${vep_cache}"   : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: ""
    ].join(' ').trim()

    def bcftools_args = cytobands ? "${cytobands.min{ it.toString().length() }}" : ""
    """
    set -eo pipefail

    gunzip -c ${vcf} \\
        | awk -v FS="\t" -v OFS="\t" '{ if(\$5=="<DEL>,<DUP>"){ \$5="<CNV>"; } print; }' \\
        | bgzip -c > vep_input.vcf.gz

    /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --offline \\
        -i vep_input.vcf.gz \\
        --format vcf \\
        --fields SYMBOL \\
        ${vep_gene_args} \\
        -o STDOUT \\
        --max_sv_size 300000000 \\
        --vcf_info_field VEPGENES \\
    | bcftools annotate \\
        -a "${bcftools_args}" \\
        -c CHROM,BEG,END,INFO/Cytobands,- \\
        -H '##INFO=<ID=Cytobands,Number=.,Type=String,Description="Cytobands">' \\
        -l Cytobands:append \\
        | awk -v FS="\t" -v OFS="\t" '{ if(\$5=="<CNV>"){ \$5="<DEL>,<DUP>"; } print; }' \\
        | bgzip -c > "${meta.id}.cnv.annotated.vcf.gz"

    tabix -p vcf "${meta.id}.cnv.annotated.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | awk -F ': ' '{print \$NF}')
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def vcf = dragen_files.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def vep_gene_args = [
        vep_cache                                 ? "--dir ${vep_cache}"   : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" } ?: ""
    ].join(' ').trim()

    def bcftools_args = cytobands ? "${cytobands.min{ it.toString().length() }}" : ""
    """
    set -eo pipefail

    touch \\
        "${meta.id}.cnv.annotated.vcf.gz" \\
        "${meta.id}.cnv.annotated.vcf.gz.tbi"

    cat <<-END_CMDS > "${meta.id}_cmds.txt"
    gunzip -c ${vcf} \\
        | awk -v FS="\t" -v OFS="\t" '{ if(5=="<DEL>,<DUP>"){5="<CNV>"; } print; }' \\
        | bgzip -c > vep_input.vcf.gz

    /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --offline \\
        -i vep_input.vcf.gz \\
        --format vcf \\
        --fields SYMBOL \\
        ${vep_gene_args} \\
        -o STDOUT \\
        --max_sv_size 300000000 \\
        --vcf_info_field VEPGENES \\
    | bcftools annotate \\
        -a "${bcftools_args}" \\
        -c CHROM,BEG,END,INFO/Cytobands,- \\
        -H '##INFO=<ID=Cytobands,Number=.,Type=String,Description="Cytobands">' \\
        -l Cytobands:append \\
        | awk -v FS="\t" -v OFS="\t" '{ if(5=="<CNV>"){ 5="<DEL>,<DUP>"; } print; }' \\
        | bgzip -c > "${meta.id}.cnv.annotated.vcf.gz"

    tabix -p vcf "${meta.id}.cnv.annotated.vcf.gz"
    END_CMDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | awk -F ': ' '{print \$NF}')
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
