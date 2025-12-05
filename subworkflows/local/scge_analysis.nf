/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_VARIANTS             } from '../../modules/local/annotate_variants.nf'
include { VEP_TO_TSV as VARIANTS_TO_TSV } from '../../modules/local/vep_to_tsv.nf'
include { ANNOTATE_SV_VARIANTS          } from '../../modules/local/annotate_sv_variants.nf'
include { VEP_TO_TSV as SV_TO_TSV       } from '../../modules/local/vep_to_tsv.nf'
include { ANNOTATE_CNV_VARIANTS         } from '../../modules/local/annotate_cnv_variants.nf'
include { VEP_TO_TSV as CNV_TO_TSV      } from '../../modules/local/vep_to_tsv'
include { ANNOTATE_OFFTARGETS         } from '../../modules/local/annotate_offtargets.nf'
include { GET_INDELS                  } from '../../modules/local/get_indels.nf'
include { GET_TRANSGENE_JUNCTIONS     } from '../../modules/local/get_transgene_junctions.nf'
include { TRANSGENE_TO_VCF            } from '../../modules/local/transgene_to_vcf'
include { ANNOTATE_TRANSGENE_VARIANTS } from '../../modules/local/annotate_transgene.nf'
include { TRANSFORM_TRANSGENE         } from '../../modules/local/transform_transgene.nf'
include { MAKE_CIRCOS_PLOT            } from '../../modules/local/make_circos_plot.nf'
include { REFORMAT_CNV_DATA           } from '../../modules/local/reformat_cnv_data.nf'
include { GENERATE_CNA_BAF_PLOTS      } from '../../modules/local/generate_cna_baf_plots.nf'
include { COMPILE_REPORT_JSON         } from '../../modules/local/compile_report_json.nf'
include { BND_FROM_INDELS_TO_VCF      } from '../../modules/local/bnd_from_indels_to_vcf.nf'
include { MAKE_SCGE_REPORT            } from './make_scge_report.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Vep cache
ch_vepcache = params.vepcache
    ? Channel.fromPath(params.vepcache, checkIfExists: true)
    : Channel.empty()

ch_cytobands = params.cytobands
    ? Channel.fromPath("${params.cytobands}*", checkIfExists: true).collect()
    : []

ch_crispr_model = params.crispr_model ?
    Channel.fromPath("${params.crispr_model}", checkIfExists: true) 
    : []

/*
========================================================================================
    SUBWORKFLOW TO ANALYZE DATA
========================================================================================
*/

workflow SCGE_ANALYSIS {

    take:
    ch_analysis_samples // channel: meta
    ch_fasta_reference // channel: path(fasta_reference)

    main:
    ch_versions = Channel.empty()

    //
    // Main analysis
    //

    // Extract dragenfiles and target files from analysis samples
    ch_dragen_files = ch_analysis_samples
        .map { meta, dragenfiles, targetfile -> [ meta, dragenfiles ]}

    ANNOTATE_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    VARIANTS_TO_TSV (ANNOTATE_VARIANTS.out.vcf,Channel.value('vcf'))
    ch_versions = ch_versions.mix(VARIANTS_TO_TSV.out.versions)

    ANNOTATE_SV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    SV_TO_TSV (ANNOTATE_SV_VARIANTS.out.vcf,Channel.value('sv'))
    ch_versions = ch_versions.mix(SV_TO_TSV.out.versions)

    ANNOTATE_CNV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_CNV_VARIANTS.out.versions)

    CNV_TO_TSV (ANNOTATE_CNV_VARIANTS.out.vcf,Channel.value('cnv'))
    ch_versions = ch_versions.mix(CNV_TO_TSV.out.versions)

    ANNOTATE_OFFTARGETS(
        ch_analysis_samples.map{ meta, dragenfiles, targetfile -> [meta, targetfile] },
        ch_vepcache,
        ch_fasta_reference
    )
    ch_versions = ch_versions.mix(ANNOTATE_OFFTARGETS.out.versions)

    GET_INDELS(ch_dragen_files.join(ANNOTATE_OFFTARGETS.out.targetfile),
            ch_crispr_model)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

    GET_TRANSGENE_JUNCTIONS(ch_dragen_files,
                            ch_fasta_reference)
    ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

    TRANSGENE_TO_VCF(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ANNOTATE_TRANSGENE_VARIANTS(
        TRANSGENE_TO_VCF.out.transgene_vcf,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_VARIANTS.out.versions)

    TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

    MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)

/*
    BND_FROM_INDELS_TO_VCF (
        GET_INDELS.out.indels_file
            .join(VEP_TO_TSV.out.vep_tsv)
            .map { meta, indels_file, vep_tsv -> [meta, indels_file] }
    )
*/
    //
    // Generate plots
    //
    GENERATE_CNA_BAF_PLOTS(ch_dragen_files)
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)

    //
    // Collate outputs
    //
    /*
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
    
    */
     
    CUSTOM_DUMPSOFTWAREVERSIONS (
         ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    emit:
    versions = ch_versions  // channel: [ path(file) ]

}

