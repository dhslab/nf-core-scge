process ANNOTATE_SV_VARIANTS {
    tag "${meta.id}"
    label "process_medium"
    container "ghcr.io/dhslab/docker-vep_release113:250810"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_files, stageAs: "dragen_files/*")
    path(reference)
    path(vep_cache)
    path(cytobands)

    output:
    tuple val(meta), path("*.sv.annotated.vcf.gz*"), emit: vcf
    path("versions.yml")                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def vep_args = [
        vep_cache                                 ? "--dir ${vep_cache}"                            : "",
        params.sv_annotation_distance             ? "--distance ${params.sv_annotation_distance}"   : "",
        params.max_filter_sv_length               ? "--max_sv_size ${params.max_filter_sv_length}"  : "",
        task.cpus > 1                             ? "--fork ${task.cpus}"                           : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" }                          ?: ""
    ].join(' ').trim()

    def bcftools_cytobands     = cytobands           ? "${cytobands.min{ it.toString().length() }}"           : ""
    """
    dragen_sv_file=\$(find -L dragen_files/ -type f -name "*.sv*.vcf.gz" | head -n 1)

    if [ -e dragen_files/*.dux4.vcf.gz ]; then
        bcftools concat \\
            -a "\${dragen_sv_file}" dragen_files/*.dux4.vcf.gz \\
            | bcftools sort \\
                --write-index \\
                -Oz \\
                -o "${meta.id}.merged.vcf.gz"
    else
        cp \${dragen_sv_file} "${meta.id}.merged.vcf.gz"

        if [[ -f "\${dragen_sv_file}.tbi" ]]; then
            cp \${dragen_sv_file}.tbi "${meta.id}.merged.vcf.gz.tbi"
        else
            exit 1
        fi
    fi

    bcftools view -i 'SVTYPE!="BND"' "${meta.id}.merged.vcf.gz" \\
    | /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --symbol \\
        --term SO \\
        --offline \\
        ${vep_args} \\
        --flag_pick \\
        --format vcf \\
        --dir_plugins /opt/lib/perl/VEP/Plugins/ \\
        --plugin TranscriptCoordinates,START,END \\
        --compress_output bgzip \\
        -o vep_cnvs_out.vcf.gz

    tabix -p vcf vep_cnvs_out.vcf.gz

    bcftools view -i 'SVTYPE=="BND"' "${meta.id}.merged.vcf.gz" \\
    | /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --symbol \\
        --term SO \\
        --offline \\
        ${vep_args} \\
        --flag_pick \\
        --format vcf \\
        --dir_plugins /opt/lib/perl/VEP/Plugins/ \\
        --plugin TranscriptCoordinates,START,END \\
        --compress_output bgzip \\
        -o vep_bnds_out.vcf.gz

    tabix -p vcf vep_bnds_out.vcf.gz

    bcftools concat -a vep_cnvs_out.vcf.gz vep_bnds_out.vcf.gz \\
    | bcftools sort \\
    | bcftools annotate -a ${bcftools_cytobands} \\
        -c CHROM,BEG,END,INFO/Cytobands,- \\
        -H '##INFO=<ID=Cytobands,Number=.,Type=String,Description="Cytobands">' \\
        -l Cytobands:append \\
        -Oz \\
        -o "${meta.id}.sv.annotated.vcf.gz"

    tabix -p vcf "${meta.id}.sv.annotated.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | awk -F ': ' '{print \$NF}')
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def vep_args = [
        vep_cache                                 ? "--dir ${vep_cache}"                            : "",
        params.sv_annotation_distance             ? "--distance ${params.sv_annotation_distance}"   : "",
        params.max_filter_sv_length               ? "--max_sv_size ${params.max_filter_sv_length}"  : "",
        task.cpus > 1                             ? "--fork ${task.cpus}"                           : "",
        reference.find{ it ==~ /.*\.(fasta|fa)$/ }?.with{ "--fasta $it" }                          ?: ""
    ].join(' ').trim()

    def bcftools_cytobands     = cytobands           ? "${cytobands.min{ it.toString().length() }}"           : ""
    """
    set -eo pipefail

    touch \\
        "${meta.id}.sv.annotated.vcf.gz" \\
        "${meta.id}.sv.annotated.vcf.gz.tbi"

    cat <<-END_CMDS > "${meta.id}_cmds.txt"
    dragen_sv_file=\$(find -L dragen_files/ -type f -name "*.sv*.vcf.gz" | head -n 1)

    if [ -e dragen_files/*.dux4.vcf.gz ]; then
        bcftools concat \\
            -a "{dragen_sv_file}" dragen_files/*.dux4.vcf.gz \\
            | bcftools sort \\
                --write-index \\
                -Oz \\
                -o "${meta.id}.merged.vcf.gz"
    else
        cp {dragen_sv_file} "${meta.id}.merged.vcf.gz"

        if [[ -f "{dragen_sv_file}.tbi" ]]; then
            cp {dragen_sv_file}.tbi "${meta.id}.merged.vcf.gz.tbi"
        else
            exit 1
        fi
    fi

    bcftools view -i 'SVTYPE!="BND"' "${meta.id}.merged.vcf.gz" \\
    | /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --symbol \\
        --term SO \\
        --offline \\
        ${vep_args} \\
        --flag_pick \\
        --format vcf \\
        --dir_plugins /opt/lib/perl/VEP/Plugins/ \\
        --plugin TranscriptCoordinates,START,END \\
        --compress_output bgzip \\
        -o vep_cnvs_out.vcf.gz

    tabix -p vcf vep_cnvs_out.vcf.gz

    bcftools view -i 'SVTYPE=="BND"' "${meta.id}.merged.vcf.gz" \\
    | /opt/vep/src/ensembl-vep/vep \\
        --vcf \\
        --cache \\
        --symbol \\
        --term SO \\
        --offline \\
        ${vep_args} \\
        --flag_pick \\
        --format vcf \\
        --dir_plugins /opt/lib/perl/VEP/Plugins/ \\
        --plugin TranscriptCoordinates,START,END \\
        --compress_output bgzip \\
        -o vep_bnds_out.vcf.gz

    tabix -p vcf vep_bnds_out.vcf.gz

    bcftools concat -a vep_cnvs_out.vcf.gz vep_bnds_out.vcf.gz \\
    | bcftools sort \\
    | bcftools annotate -a ${bcftools_cytobands} \\
        -c CHROM,BEG,END,INFO/Cytobands,- \\
        -H '##INFO=<ID=Cytobands,Number=.,Type=String,Description="Cytobands">' \\
        -l Cytobands:append \\
        -Oz \\
        -o "${meta.id}.sv_vep.vcf.gz"

    tabix -p vcf "${meta.id}.sv_vep.vcf.gz"
    END_CMDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | awk -F ': ' '{print \$NF}')
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
