/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_VARIANTS           } from '../../modules/local/annotate_variants.nf'
include { ANNOTATE_TRANSGENE_VARIANTS } from '../../modules/local/annotate_transgene.nf'
include { GET_INDELS                  } from '../../modules/local/get_indels.nf'
include { GET_TRANSGENE_JUNCTIONS     } from '../../modules/local/get_transgene_junctions.nf'
include { REFORMAT_CNV_DATA           } from '../../modules/local/reformat_cnv_data.nf'
include { VEP_TO_TSV                  } from '../../modules/local/vep_to_tsv.nf'
include { GENERATE_CNA_BAF_PLOTS      } from '../../modules/local/generate_cna_baf_plots.nf'
include { COMPILE_REPORT_JSON         } from '../../modules/local/compile_report_json.nf'
include { MAKE_CIRCOS_PLOT            } from '../../modules/local/make_circos_plot.nf'
include { TRANSFORM_TRANSGENE         } from '../../modules/local/transform_transgene.nf'
include { ANNOTATE_OFFTARGETS         } from '../../modules/local/annotate_offtargets.nf'
include { BND_FROM_INDELS_TO_VCF      } from '../../modules/local/bnd_from_indels_to_vcf.nf'
include { MAKE_SCGE_REPORT            } from './make_scge_report.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRANSGENE_TO_VCF            } from '../../modules/local/transgene_to_vcf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Vep cache
ch_vep_cache = params.vep_cache
    ? Channel.fromPath(params.vep_cache, checkIfExists: true)
    : Channel.empty()

// Gene regions
ch_editing_targets = params.editing_targets ?
        Channel.fromPath("${params.editing_targets}", checkIfExists: true) :
        Channel.empty()

// Transgene name
ch_transgene_name = params.transgene ?
        Channel.from("${params.transgene}") 
        : Channel.empty()

ch_transgene_fasta = params.transgene_fasta ?
        Channel.fromPath("${params.transgene_fasta}", checkIfExists: true)
        : Channel.empty()

ch_hotspot_file = params.hotspot_file ?
        Channel.fromPath("${params.hotspot_file}", checkIfExists: true)
        : Channel.empty()

/*
========================================================================================
    SUBWORKFLOW TO ANALYZE DATA
========================================================================================
*/

workflow SCGE_ANALYSIS {

    take:
    ch_dragen_files // channel: meta
    ch_hotspot_file // channel: path(hotspot_file)
    ch_fasta_reference // channel: path(fasta_reference)
    ch_vep_cache // channel: path(vep_cache)
    ch_crispr_model // channel: path(crispr_model)

    main:
    ch_versions = Channel.empty()

    //
    // Main analysis
    //
    ch_annotate_variants_input = ch_dragen_files
        .combine(ch_fasta_reference)
        .combine(ch_vep_cache)

    ANNOTATE_VARIANTS(ch_annotate_variants_input)
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    VEP_TO_TSV(ANNOTATE_VARIANTS.out.vcf.map { meta, vcf -> [meta, "vcf", vcf] })
    ch_versions = ch_versions.mix(VEP_TO_TSV.out.versions)

    ch_get_indels_input = ch_dragen_files
        .combine(ch_hotspot_file)
        .combine(ch_crispr_model)

    GET_INDELS(ch_get_indels_input)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

    ANNOTATE_OFFTARGETS(GET_INDELS.out.indels_file)
    ch_versions = ch_versions.mix(ANNOTATE_OFFTARGETS.out.versions)

    ch_get_transgene_junctions_input = ch_dragen_files
        .combine(ch_transgene_name)
        .combine(ch_transgene_fasta)

    GET_TRANSGENE_JUNCTIONS(ch_get_transgene_junctions_input)
    ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

    TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

    MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)

    ch_transgene_fasta = Channel.fromPath(params.transgene_fasta)

    TRANSGENE_TO_VCF( GET_TRANSGENE_JUNCTIONS.out.transgene_file )
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ch_annotate_transgene_variants_input = ANNOTATE_VARIANTS.out.vcf
        .join(TRANSGENE_TO_VCF.out.transgene_vcf)

    ANNOTATE_TRANSGENE_VARIANTS(ch_annotate_transgene_variants_input)
    ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_VARIANTS.out.versions)

    BND_FROM_INDELS_TO_VCF (
        GET_INDELS.out.indels_file
            .join(VEP_TO_TSV.out.vep_tsv)
            .map { meta, indels_file, vep_tsv -> [meta, indels_file] }
    )

    //
    // Generate plots
    //
    GENERATE_CNA_BAF_PLOTS(ch_dragen_files)
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)


    //
    // Collate outputs
    //
    ch_coverage_files = ch_dragen_files.map { meta, dragen_path ->
        def tumor_cov_file = file(dragen_path).listFiles().find { it.name.endsWith('.wgs_overall_mean_cov_tumor.csv') } ?: file("${baseDir}/assets/empty_tumor_coverage.txt")
        def normal_cov_file = file(dragen_path).listFiles().find { it.name.endsWith('.wgs_overall_mean_cov_normal.csv') } ?: file("${baseDir}/assets/empty_normal_coverage.txt")
        return [meta.id, tumor_cov_file, normal_cov_file]
    }

    def ch_plots = GENERATE_CNA_BAF_PLOTS.out.cna_plot
        .join(GENERATE_CNA_BAF_PLOTS.out.baf_plot)
        .map { meta, cna, baf -> [meta.id, meta, cna, baf] }

    def ch_circos = MAKE_CIRCOS_PLOT.out.circos_plot
        .map { meta, circos -> [meta.id, circos] }
    ch_circos.view { "Circos: $it" }


    def ch_annotated_transgene = ANNOTATE_TRANSGENE_VARIANTS.out.annotated_transgene_variants
        .map { meta, transgene -> [meta.id, transgene] }
    ch_annotated_transgene.view { "Annotated Transgene: $it" }

    def ch_vep_tsv = VEP_TO_TSV.out.vep_tsv
        .map { meta, tsv -> [meta.id, tsv] }
    ch_vep_tsv.view { "VEP TSV: $it" }

    def ch_indels = ANNOTATE_OFFTARGETS.out.annotated_indels
        .map { meta, indels -> [meta.id, indels] }
    ch_indels.view { "Indels: $it" }

    ch_plots
        .join(ch_circos, by: 0)
        .join(ch_annotated_transgene, by: 0)
        .join(ch_vep_tsv, by: 0)
        .join(ch_indels, by: 0)
        .join(ch_coverage_files, by: 0)
        .map { id, meta, cna, baf, circos, transgene, tsv, indels, tumor_cov, normal_cov ->
            def timestamp = new Date().getTime()
            [meta, cna, baf, circos, transgene, tsv, indels, tumor_cov, normal_cov, timestamp]
        }
        .set { ch_compile_report_input }

    COMPILE_REPORT_JSON(ch_compile_report_input)
    ch_versions = ch_versions.mix(COMPILE_REPORT_JSON.out.versions)

    def ch_plots_for_report = GENERATE_CNA_BAF_PLOTS.out.cna_plot
        .join(GENERATE_CNA_BAF_PLOTS.out.baf_plot)
        .map { meta, cna, baf -> [meta.id, cna, baf] }

    def ch_report_input = COMPILE_REPORT_JSON.out.json
        .map { meta, json -> [meta.id, meta, json] }
        .join(ch_plots_for_report)
        .map { id, meta, json, cna, baf -> [meta, json, cna, baf] }

    MAKE_SCGE_REPORT(ch_report_input)

    CUSTOM_DUMPSOFTWAREVERSIONS (
         ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    emit:
    versions = ch_versions  // channel: [ path(file) ]

}

