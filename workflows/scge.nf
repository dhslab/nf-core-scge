/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

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
include { ANNOTATE_OFFTARGETS         } from '../modules/local/annotate_offtargets.nf'
include { BND_FROM_INDELS_TO_VCF      } from '../modules/local/bnd_from_indels_to_vcf.nf'
include { TRANSGENE_TO_VCF            } from '../modules/local/transgene_to_vcf.nf'

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

def generateMetaFromCsv(csv_file) {
    def lines = csv_file.text.readLines()
    def headers = lines[0].split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
    return lines.drop(1).collect { line ->
        def fields = line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
        [headers, fields].transpose().collectEntries { k, v -> v ? [(k): v] : [:] }
    }.findAll { it }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCGE {

    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def summary_params = paramsSummaryMap(workflow)

    // Print parameter summary log to screen
    log.info logo + paramsSummaryLog(workflow) + citation

    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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

    def multiqc_report = []

    ch_versions = Channel.empty()
    ch_dragen_outputs = Channel.empty()

    //
    // *** Main inputs can be made into channels here. ***
    //
    // - reference fasta file (and fai)
    // - hotspot file (now a parameter and not in the input sheet)
    // - VEP cache
    // FastA reference
    ch_fasta_reference = Channel.fromPath(params.fasta, checkIfExists: true)
    // Vep cache
    ch_vep_cache = params.vep_cache ?
        Channel.fromPath(params.vep_cache, checkIfExists: true).collect() :
        Channel.empty()
    // Gene regions
    ch_editing_targets = params.editing_targets ?
        Channel.fromPath("${params.editing_targets}", checkIfExists: true) :
        Channel.empty()
    // Gene regions
    ch_transgene_name = params.transgene_name ?
        Channel.from("${params.transgene_name}") :
        Channel.empty()
    // input mastersheet
    ch_mastersheet = params.input ?
        Channel.fromPath("${params.input}", checkIfExists: true) :
        Channel.empty()

    //
    // This input check needs to be overhauled. See assets/stub/sample_mastersheet.csv
    //

    // Channel of meta data for alignment samples
    ch_samples = ch_mastersheet
                        .map { generateMetaFromCsv(it) }
                        .flatten()

    // This pseudo code at this point, but its supposed to get the dragen_path out of the dict
    // and return a tuple of (dict,path)
    ch_dragen_output = ch_samples
    .map { dict ->
        if (dict.dragen_path && dict.dragen_path != '') {
            def path = file(dict.dragen_path)
            if (path.exists()) {
                return [dict, path]
            } else {
                log.warn "DRAGEN path does not exist for sample ${dict.id}: ${dict.dragen_path}"
                return null
            }
        } else {
            log.warn "No DRAGEN path specified for sample ${dict.id}"
            return null
        }
    }
    .filter { it != null }

    // this stages all dragen output and then returns the images
    GENERATE_CNA_BAF_PLOTS(ch_dragen_output)
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)

    // get indels
    // Build input for GET_INDELS by joining dragen outputs with per-sample hotspot file
    ch_hotspot_file = ch_samples
    .map { dict ->
        if (dict.hotspot_file && dict.hotspot_file != '') {
            def path = file(dict.hotspot_file)
            if (path.exists()) {
                return [dict.id, path]
            } else {
                log.warn "Hotspot file does not exist for sample ${dict.id}: ${dict.hotspot_file}"
                return [dict.id, null]
            }
        } else {
            log.warn "No hotspot file specified for sample ${dict.id}"
            return [dict.id, null]
        }
    }
    .filter { it[1] != null }


    ch_dragen_output_for_join = ch_dragen_output.map { meta, files -> [meta.id, meta, files] }

    ch_get_indels_input = ch_dragen_output_for_join
        .join(ch_hotspot_file)
        .map { id, meta, files, hotspot_file -> [meta, files, hotspot_file] }

    GET_INDELS(ch_get_indels_input)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

    ANNOTATE_OFFTARGETS(GET_INDELS.out.indels_file)
    ch_versions = ch_versions.mix(ANNOTATE_OFFTARGETS.out.versions)

    ch_vcf_for_annotation = ch_dragen_output.map { meta, dragen_path ->
        def vcf_file = file("${dragen_path}/${meta.id}.hard-filtered.vcf.gz")
        if (vcf_file.exists()) {
            return [meta, vcf_file]
        } else {
            log.warn "VCF file not found for sample ${meta.id}: ${vcf_file}"
            return null
        }
    }.filter { it != null }

    ANNOTATE_VARIANTS (ch_vcf_for_annotation, ch_fasta_reference.first(), ch_vep_cache.first())
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    ANNOTATE_VARIANTS.out.vcf.view()
    VEP_TO_TSV(ANNOTATE_VARIANTS.out.vcf.map { meta, vcf -> [meta, "vcf", vcf] })
        ch_versions = ch_versions.mix(VEP_TO_TSV.out.versions)

    GET_TRANSGENE_JUNCTIONS(ch_dragen_output)
    ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

    TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

    MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)

    ch_transgene_fasta = Channel.fromPath(params.transgene_fasta)

    TRANSGENE_TO_VCF(
        GET_TRANSGENE_JUNCTIONS.out.transgene_file,
        ch_transgene_fasta
    )
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ch_annotate_transgene_variants_input = ch_dragen_output.join(TRANSGENE_TO_VCF.out.vcf)

    ANNOTATE_TRANSGENE_VARIANTS(ch_annotate_transgene_variants_input)
    ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_VARIANTS.out.versions)

    // Temporarily disabled to avoid early failure when no versions are present
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    ch_coverage_files = ch_dragen_output.map { meta, dragen_path ->
        def tumor_cov_file = file("${dragen_path}/${meta.id}.wgs_overall_mean_cov_tumor.csv")
        def normal_cov_file = file("${dragen_path}/${meta.id}.wgs_overall_mean_cov_normal.csv")
        if (tumor_cov_file.exists() && normal_cov_file.exists()) {
            return [meta.id, tumor_cov_file, normal_cov_file]
        } else {
            if (!tumor_cov_file.exists()) log.warn "Tumor coverage metrics file not found for sample ${meta.id}: ${tumor_cov_file}"
            if (!normal_cov_file.exists()) log.warn "Normal coverage metrics file not found for sample ${meta.id}: ${normal_cov_file}"
            return [meta.id, null, null]
        }
    }
    .filter { it[1] != null && it[2] != null }

    def ch_plots = GENERATE_CNA_BAF_PLOTS.out.cna_plot
        .join(GENERATE_CNA_BAF_PLOTS.out.baf_plot)
        .map { meta, cna, baf -> [meta.id, meta, cna, baf] }
    ch_plots.view { "Plots: $it" }

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

    // MODULE: MultiQC
    workflow_summary    = WorkflowScge.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowScge.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
