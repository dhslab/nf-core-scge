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
include { ANNOTATE_OFFTARGETS           } from '../../modules/local/annotate_offtargets.nf'
include { GET_INDELS                    } from '../../modules/local/get_indels.nf'
include { GET_TRANSGENE_JUNCTIONS       } from '../../modules/local/get_transgene_junctions.nf'
include { TRANSGENE_TO_VCF              } from '../../modules/local/transgene_to_vcf'
include { ANNOTATE_TRANSGENE_JUNCTIONS } from '../../modules/local/annotate_transgene_junctions.nf'
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
    ? Channel.fromPath(params.vepcache, type: 'dir', checkIfExists: true)
    : Channel.empty()

ch_cytobands = params.cytobands
    ? Channel.fromPath("${params.cytobands}*", checkIfExists: true).collect()
    : []

ch_crispr_model = params.crispr_model ?
    Channel.fromPath("${params.crispr_model}", checkIfExists: true) 
    : []

ch_transgene_name = params.transgene_name && params.transgene_name != false && params.transgene_name != null
    ? Channel.value(params.transgene_name) : Channel.empty()


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
    ch_report_inputs = Channel.empty() // channel: meta, [ file1, file2, etc ]

    //
    // Main analysis
    //

    // Extract dragenfiles and target files from analysis samples
    ch_dragen_files = ch_analysis_samples
        .map { meta, dragenfiles, targetfile -> [ meta, dragenfiles ]}

    // Get coverage files for report
    ch_report_inputs = ch_report_inputs.mix(
            ch_dragen_files.map { meta, dragen_path ->
                def tumor_cov  = dragen_path.find { it.name.endsWith('.wgs_overall_mean_cov_tumor.csv') }
                def normal_cov = dragen_path.find { it.name.endsWith('.wgs_overall_mean_cov_normal.csv') }            
                return [ meta, [ file(tumor_cov,checkIfExists: true), file(normal_cov,checkIfExists: true) ] ]
            }
        )
 
    // Annotate small variants
    ANNOTATE_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    VARIANTS_TO_TSV (ANNOTATE_VARIANTS.out.vcf,Channel.value('vcf'))
    ch_report_inputs = ch_report_inputs.mix(VARIANTS_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(VARIANTS_TO_TSV.out.versions)

    // Annotate structural variants
    ANNOTATE_SV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    SV_TO_TSV (ANNOTATE_SV_VARIANTS.out.vcf,Channel.value('sv'))
    ch_report_inputs = ch_report_inputs.mix(SV_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(SV_TO_TSV.out.versions)

    // Annotate copy number variants
    ANNOTATE_CNV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_CNV_VARIANTS.out.versions)

    CNV_TO_TSV (ANNOTATE_CNV_VARIANTS.out.vcf,Channel.value('cnv'))
    ch_report_inputs = ch_report_inputs.mix(CNV_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(CNV_TO_TSV.out.versions)

    // Annotate and analyze off-target sites
    ANNOTATE_OFFTARGETS(
        ch_analysis_samples.map{ meta, dragenfiles, targetfile -> [meta, targetfile] },
        ch_vepcache,
        ch_fasta_reference
    )
    ch_versions = ch_versions.mix(ANNOTATE_OFFTARGETS.out.versions)

    GET_INDELS(
        ch_dragen_files.join(ANNOTATE_OFFTARGETS.out.targetfile),
        ch_crispr_model,
        ch_fasta_reference
    )

    ch_report_inputs = ch_report_inputs.mix(GET_INDELS.out.indels_file)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

    BND_FROM_INDELS_TO_VCF (
        GET_INDELS.out.indels_file
    )
    ch_report_inputs = ch_report_inputs.mix(BND_FROM_INDELS_TO_VCF.out.vcf)
    ch_versions = ch_versions.mix(BND_FROM_INDELS_TO_VCF.out.versions)

    // Get transgene junctions
    GET_TRANSGENE_JUNCTIONS(ch_dragen_files,
                            ch_fasta_reference,
                            ch_transgene_name
    )
    ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

    TRANSGENE_TO_VCF(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ANNOTATE_TRANSGENE_JUNCTIONS(
        TRANSGENE_TO_VCF.out.transgene_vcf,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_report_inputs = ch_report_inputs.mix(ANNOTATE_TRANSGENE_JUNCTIONS.out.annotated_junctions)
    ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_JUNCTIONS.out.versions)

    // Make Circos plot
    TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

    MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)
    // Add to report inputs
    ch_report_inputs = ch_report_inputs.mix(MAKE_CIRCOS_PLOT.out.plot)

    //
    // Generate CNA plots
    //
    GENERATE_CNA_BAF_PLOTS(ch_dragen_files)
    ch_report_inputs = ch_report_inputs.mix(GENERATE_CNA_BAF_PLOTS.out.plots) 
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)

    // Make report JSON
    COMPILE_REPORT_JSON(
        ch_report_inputs
        .groupTuple()
        .map{ it -> [ it[0], it[1].flatten() ] },
        Channel.value("${new Date().getTime()}"))
    ch_versions = ch_versions.mix(COMPILE_REPORT_JSON.out.versions)

/*
    // Render report
    MAKE_SCGE_REPORT(
        COMPILE_REPORT_JSON.out.json
        .join(GENERATE_CNA_BAF_PLOTS.out.plots)
    )
*/
    ch_versions.unique().collectFile(name: 'collated_versions.yml').view()
//    CUSTOM_DUMPSOFTWAREVERSIONS (
//         ch_versions.unique().collectFile(name: 'collated_versions.yml').view()
//    )

    emit:
    versions = ch_versions  // channel: [ path(file) ]

}

