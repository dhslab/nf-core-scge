
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_scge_pipeline'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { PARSE_INPUT_SAMPLESHEET     } from '../modules/local/parse_input_samplesheet.nf'
include { GATHER_ALIGNMENT_SAMPLES    } from '../subworkflows/local/gather_alignment_samples.nf'
include { MAKE_HOTSPOT_VCF            } from '../modules/local/make_hotspot_vcf.nf'
include { DRAGEN_SCGE                 } from '../modules/local/dragen_scge.nf'
include { SCGE_ANALYSIS               } from '../subworkflows/local/scge_analysis.nf'

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
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// DRAGEN reference directory
ch_reference_dir = params.refdir
    ? Channel.fromPath(params.refdir, type: 'dir', checkIfExists: true).collect()
    : Channel.empty()

// DRAGEN adapter sequences for read 1
ch_adapter1_file = params.adapter1
    ? Channel.fromPath(params.adapter1, checkIfExists: true).collect()
    : []

// DRAGEN adapter sequences for read 2
ch_adapter2_file = params.adapter2
    ? Channel.fromPath(params.adapter2, checkIfExists: true).collect()
    : []

// DRAGEN intermediate directory
if (params.intermediate_dir?.toString()?.startsWith('/staging')) {
    ch_intermediate_dir = Channel.of(params.intermediate_dir).map{ [ it, [] ] }.collect()
} else if (params.intermediate_dir) {
    ch_intermediate_dir = Channel.fromPath(params.intermediate_dir).map{ [ [], it ] }.collect()
} else {
    ch_intermediate_dir = [ [], [] ]
}

// Gene hotspots
ch_hotspot_bed = params.hotspot_bed
    ? Channel.fromPath("${params.hotspot_bed}", checkIfExists: true).collect()
    : []
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
    ch_crispr_model = params.crispr_model ?
        Channel.fromPath("${params.crispr_model}", checkIfExists: true) :
        Channel.empty()

// SNV systematic noise BED file
ch_snv_noisefile = params.snv_noisefile
    ? Channel.fromPath(params.snv_noisefile, checkIfExists: true).collect()
    : []

// SV systematic noise BED file
ch_sv_noisefile = params.sv_noisefile
    ? Channel.fromPath(params.sv_noisefile, checkIfExists: true).collect()
    : []

// High confidence CNV VCF file
ch_cnv_population_vcf = params.cnv_population_vcf
    ? Channel.fromPath(params.cnv_population_vcf, checkIfExists: true).collect()
    : []

// CRAM reference file
ch_cram_reference = params.cram_reference
    ? Channel.fromPath("${params.cram_reference}*", checkIfExists: true).collect()
    : []

ch_fasta_reference = params.fasta
    ? Channel.fromPath("${params.fasta}*", checkIfExists: true).collect()
    : Channel.empty()

// Vep cache
ch_vep_cache = params.vep_cache
    ? Channel.fromPath(params.vep_cache, checkIfExists: true).collect()
    : Channel.empty()
    GET_INDELS(ch_get_indels_input, ch_crispr_model)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

/*
~~~~~~~~~~~~~~~~~~
MultiQC parameters
~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

ch_multiqc_custom_config = params.multiqc_config 
    ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) 
    : Channel.empty()
    ANNOTATE_VARIANTS (ch_vcf_for_annotation, ch_fasta_reference, ch_vep_cache)
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

ch_multiqc_logo = params.ch_multiqc_logo
    ? Channel.fromPath( params.multiqc_logo, checkIfExists: true )
    : Channel.empty()

ch_multiqc_custom_methods_description = params.multiqc_methods_description 
    ? file(params.multiqc_methods_description, checkIfExists: true) 
    : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCGE {

    take:
    ch_input_samplesheet  // channel: [ path(file) ]
    TRANSGENE_TO_VCF(
        GET_TRANSGENE_JUNCTIONS.out.transgene_file
    )
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ch_annotate_transgene_variants_input = ch_dragen_output
        .map { meta, files -> [meta.id, meta, files] }
        .join(TRANSGENE_TO_VCF.out.transgene_vcf.map { meta, vcf -> [meta.id, vcf] })
        .map { id, meta, files, vcf -> [meta, files, vcf] }

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

    main:
    ch_versions = Channel.empty()
    ch_alignment_samples = Channel.empty()
    ch_dragen_output = Channel.empty()
    ch_dragen_usage = Channel.empty()

    //
    // dump samplesheet channel
    ch_input_samplesheet.dump(tag:'mastersheet')

    //
    // MODULE: Parse input samplesheet to format samples for processing.
    //         Output of this process are csv files for samples that need to be aligned
    //         and samples that need to be analyzed
    //
    PARSE_INPUT_SAMPLESHEET (
        ch_input_samplesheet
    )
    ch_versions = ch_versions.mix(PARSE_INPUT_SAMPLESHEET.out.versions)

    PARSE_INPUT_SAMPLESHEET.out.samples_to_align.dump(tag:'alignmentsamples')
    PARSE_INPUT_SAMPLESHEET.out.samples_to_analyze.dump(tag:'analysissamples')

    ch_dragen_output = ch_dragen_output.mix(
        PARSE_INPUT_SAMPLESHEET.out.samples_to_analyze
            .map{ generateMetaFromCsv(it) }
            .flatten()
            .filter{ it.dragen_path }
            .map{ [ it, it.dragen_path ? file("${it.dragen_path}/*") : [] ] }
    )

    //
    // SUBWORKFLOW: Gather alignment samples
    //
    GATHER_ALIGNMENT_SAMPLES (
        PARSE_INPUT_SAMPLESHEET.out.samples_to_align,
        ch_cram_reference
    )
    ch_versions          = ch_versions.mix(GATHER_ALIGNMENT_SAMPLES.out.versions)
    ch_alignment_samples = ch_alignment_samples.mix(GATHER_ALIGNMENT_SAMPLES.out.samples)

    ch_alignment_samples.dump(tag: 'ch_alignment_samples', pretty: true)

    //
    // Run dragen
    //
    if (params.run_alignment) {

        //
        // Make hotspot vcf file
        //
        MAKE_HOTSPOT_VCF (
            ch_hotspot_bed,
            ch_fasta_reference
        )

        DRAGEN_SCGE (
                ch_alignment_samples,
                ch_intermediate_dir,
                ch_reference_dir,
                ch_adapter1_file,
                ch_adapter2_file,
                ch_cram_reference,
                ch_sv_noisefile,
                ch_snv_noisefile,
                ch_cnv_population_vcf,
                MAKE_HOTSPOT_VCF.out.hotspot_vcf
        )
        ch_versions     = ch_versions.mix(DRAGEN_SCGE.out.versions)
        ch_dragen_usage = ch_dragen_usage.mix(DRAGEN_SCGE.out.usage)
        ch_dragen_output = ch_dragen_output.mix(DRAGEN_SCGE.out.dragen_output)

        // Output DRAGEN usage information
        ch_dragen_usage.map{
                            def meta = it.getSimpleName().split("_usage")[0]
                            def data = it.text.split("\\: ").join('\t')
                            return "Accession\tLicense Type\tUsage\n${meta}\t${data}"
                        }
                        .collectFile(
                            name      : "DRAGEN_usage.tsv",
                            keepHeader: true,
                            storeDir  : "${params.outdir}/pipeline_info"
                        )

    }

    if (params.run_analysis) {
        SCGE_ANALYSIS (
            ch_dragen_output
        )
        ch_versions = ch_versions.mix(SCGE_ANALYSIS.out.versions)
    }

    //
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name    : 'software_versions.yml',
            sort    : true,
            newLine : true
        )
        .set{ ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_collated_versions
                        .mix(
                            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
                        )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()  // channel: [ path(file) ]
    versions       = ch_versions                  // channel: [ path(file) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
