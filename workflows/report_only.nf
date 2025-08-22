#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Include the reporting subworkflow from your existing pipeline
include { VEP_TO_TSV                  } from '../modules/local/vep_to_tsv.nf'
include { GENERATE_CNA_BAF_PLOTS      } from '../modules/local/generate_cna_baf_plots.nf'
include { COMPILE_REPORT_JSON         } from '../modules/local/compile_report_json.nf'
include { MAKE_CIRCOS_PLOT            } from '../modules/local/make_circos_plot.nf'
include { TRANSFORM_TRANSGENE         } from '../modules/local/transform_transgene.nf'
include { ANNOTATE_OFFTARGETS         } from '../modules/local/annotate_offtargets.nf'
include { BND_FROM_INDELS_TO_VCF      } from '../modules/local/bnd_from_indels_to_vcf.nf'
include { ANNOTATE_VARIANTS           } from '../modules/local/annotate_variants.nf'
include { MAKE_SCGE_REPORT            } from '../subworkflows/local/make_scge_report.nf'


workflow ANALYSIS_ONLY {
    ch_mastersheet = Channel.fromPath(params.input)

    ch_samples = ch_mastersheet
        .splitCsv(header:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id
            [ meta, file(row.dragen_path) ]
        }

    ch_fasta_reference = file(params.fasta)

    ANNOTATE_VARIANTS(ch_samples, ch_fasta_reference)

    VEP_TO_TSV(ANNOTATE_VARIANTS.out.vcf.map { meta, vcf_files -> [meta, "vcf", vcf_files] })

    GENERATE_CNA_BAF_PLOTS(ch_samples)

    // At this point, you will need to gather all the required inputs for COMPILE_REPORT_JSON.
    // This will likely involve several more steps of data processing and channel manipulation.
    // For now, I will create a placeholder channel.
    ch_for_compile = ANNOTATE_VARIANTS.out.vcf
        .join(VEP_TO_TSV.out.vep_tsv)
        .join(GENERATE_CNA_BAF_PLOTS.out.cna_plot)

    COMPILE_REPORT_JSON(ch_for_compile)

    MAKE_SCGE_REPORT(COMPILE_REPORT_JSON.out.json)
} 