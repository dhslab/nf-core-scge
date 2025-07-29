/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SOMATIC_INPUT_CHECK } from '../subworkflows/local/somatic_input_check.nf'
include { MAKE_SCGE_REPORT } from '../subworkflows/local/make_scge_report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MAKE_HOTSPOT_FILE           } from '../modules/local/make_hotspot_file.nf'
include { DRAGEN_SCGE                 } from '../modules/local/dragen_scge.nf'
include { ANNOTATE_VARIANTS           } from '../modules/local/annotate_variants.nf'
include { ANNOTATE_TRANSGENE_VARIANTS } from '../modules/local/annotate_transgene.nf'
include { GET_INDELS                  } from '../modules/local/get_indels.nf'
include { GET_TRANSGENE_JUNCTIONS     } from '../modules/local/get_transgene_junctions.nf'
include { REFORMAT_CNV_DATA           } from '../modules/local/reformat_cnv_data.nf'
include { ANNOTATE_VCF                } from '../modules/local/annotate_vcf.nf'
include { VEP_TO_TSV                  } from '../modules/local/vep_to_tsv.nf'
include { GENERATE_CNA_BAF_PLOTS      } from '../modules/local/generate_cna_baf_plots.nf'
include { COMPILE_REPORT_JSON         } from '../modules/local/compile_report_json.nf'
include { MAKE_CIRCOS_PLOT            } from '../modules/local/make_circos_plot.nf'
include { TRANSFORM_TRANSGENE         } from '../modules/local/transform_transgene.nf'

def stageFileset(Map filePathMap) {
    def basePathMap = [:]
    def filePathsList = []

    filePathMap.each { key, value ->
        if (value != null) {
            def filepath = file(value)
            if (filepath.exists()) {
                // Add basename and key to the map
                basePathMap[key] = value.split('/')[-1]
                // Add file path to the list
                filePathsList << filepath
            } else {
                println "Warning: File at '${value}' for key '${key}' does not exist."
            }
        }
    }
    return [basePathMap, filePathsList]
}

// If MGI samplesheet is used, we need to set the 
// data path because only files are given. This sets the 
// data path to the samplesheet directory, or the data_path parameter.
def data_path = ""
def mastersheet = params.input
if (params.mgi == true) {
    data_path = new File(params.input).parentFile.absolutePath
} else if (params.data_path != null){
    data_path  = params.data_path
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCGE {

    ch_versions = Channel.empty()
    ch_dragen_outputs = Channel.empty()
    hotspot_bed = params.hotspot_bed ? Channel.fromPath(params.hotspot_bed) : Channel.fromPath("$projectDir/assets/NO_FILE.bed")

    SOMATIC_INPUT_CHECK(Channel.fromPath(mastersheet), data_path)

    ch_input_data = SOMATIC_INPUT_CHECK.out.input_data
    ch_hotspots = ch_input_data.map{ meta, file -> return [meta[0]['id'], file] }

    if (params.run_dragen == true) {

        hotspot_input = ch_input_data.combine(hotspot_bed)
        MAKE_HOTSPOT_FILE(hotspot_input)
        ch_input_data = MAKE_HOTSPOT_FILE.out.hotspot_vcf
            .map{ info, hotspot_vcf ->
            def newinfo = []
            newinfo = info + [hotspot_vcf]
            newinfo
            }

        ch_dragen_outputs = ch_dragen_outputs.mix(SOMATIC_INPUT_CHECK.out.dragen_outputs)
        ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))
        ch_assay_inputs = Channel.value(stageFileset(params.assay_inputs))

        DRAGEN_SCGE (ch_input_data, ch_dragen_inputs)
        ch_versions = ch_versions.mix(DRAGEN_SCGE.out.versions)
        ch_dragen_outputs = ch_dragen_outputs.mix(DRAGEN_SCGE.out.dragen_output)
    } else {
        ch_dragen_outputs = ch_dragen_outputs.mix(SOMATIC_INPUT_CHECK.out.dragen_outputs)
        ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))
        ch_assay_inputs = Channel.value(stageFileset(params.assay_inputs))
    }

    dragen_baf = ch_dragen_outputs.map{ meta, files -> [meta, files.find { it.endsWith(".tumor.baf.bedgraph.gz") }] }
    dragen_cnv = ch_dragen_outputs.map{ meta, files -> [meta, files.find { it.endsWith(".tn.tsv.gz") }] }

    GENERATE_CNA_BAF_PLOTS(dragen_baf.join(dragen_cnv))
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)

    if (params.run_analysis == true) {
        
        ANNOTATE_VARIANTS (ch_dragen_outputs, ch_assay_inputs, params.fasta)
        ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)
        
        get_indels_input = ch_dragen_outputs
            .map{ meta -> [meta[0]['id'], meta] }
            .combine(ch_hotspots, by: 0)
            .map{id, meta, hotspot_file -> [meta[0], meta[1], hotspot_file]}
        GET_INDELS (get_indels_input)
        ch_versions = ch_versions.mix(GET_INDELS.out.versions)

        ch_dragen_outputs.dump(tag: 'ch_dragen_outputs')
        if (params.transgene_analysis == true) {
            GET_TRANSGENE_JUNCTIONS (ch_dragen_outputs)
            ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

            TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
            ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

            MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)
            ch_versions = ch_versions.mix(MAKE_CIRCOS_PLOT.out.versions)
            ch_circos_plot = MAKE_CIRCOS_PLOT.out.circos_png

            annotate_transgene_input = ch_dragen_outputs.join(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
            ANNOTATE_TRANSGENE_VARIANTS (annotate_transgene_input)
            ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_VARIANTS.out.versions)
        } else {
            ch_circos_plot = ch_dragen_outputs.map { meta, files -> [meta, "NO_FILE.png"] }
        }

        REFORMAT_CNV_DATA (ch_dragen_outputs)
        ch_versions = ch_versions.mix(REFORMAT_CNV_DATA.out.versions)

        annotate_vcf_input = ch_dragen_outputs.flatMap{ meta, files -> 
            def cnv = files.find { it.endsWith("${meta.id}.cnv.vcf.gz") }
            def sv = files.find { it.endsWith("${meta.id}.sv.vcf.gz") }
            def vcf = files.find { it.endsWith("${meta.id}.hard-filtered.vcf.gz") }
            return [[meta, "cnv", cnv], [meta, "sv", sv], [meta, "vcf", vcf]] }
        ANNOTATE_VCF(annotate_vcf_input)
        ch_versions = ch_versions.mix(ANNOTATE_VCF.out.versions)

        VEP_TO_TSV(ANNOTATE_VCF.out.annotated_vcf)
        ch_versions = ch_versions.mix(VEP_TO_TSV.out.versions)

        // Combine all inputs for the report
        report_inputs = GENERATE_CNA_BAF_PLOTS.out.cna_plot
            .join(GENERATE_CNA_BAF_PLOTS.out.baf_plot)
            .join(ch_circos_plot)
            .join(VEP_TO_TSV.out.vep_tsv)
            .join(GET_INDELS.out.indels_file)

        COMPILE_REPORT_JSON(report_inputs)
        ch_versions = ch_versions.mix(COMPILE_REPORT_JSON.out.versions)

        MAKE_SCGE_REPORT(COMPILE_REPORT_JSON.out.json)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
    // MODULE: MultiQC
    
    // workflow_summary    = WorkflowDragenmultiworkflow.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowDragenmultiworkflow.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
