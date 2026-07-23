
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

include { PARSE_INPUT_SAMPLESHEET     } from '../modules/local/parse_input_samplesheet'
include { GATHER_ALIGNMENT_SAMPLES    } from '../subworkflows/local/gather_alignment_samples.nf'
include { PREPARE_SOMATIC_FASTQS      } from '../subworkflows/local/gather_alignment_samples.nf'

include { MAKE_HOTSPOT_VCF            } from '../modules/local/make_hotspot_vcf.nf'
include { DRAGEN_SCGE                 } from '../modules/local/dragen_scge.nf'
include { SCGE_ANALYSIS               } from '../subworkflows/local/scge_analysis.nf'
// Note: TRANSGENE_TO_VCF is included and invoked inside the SCGE_ANALYSIS subworkflow
// (subworkflows/local/scge_analysis.nf), not here.

def generateMetaFromCsv(csv_string) {
    def lines = csv_string.readLines()
    def headers = lines[0].split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
    return lines.drop(1).collect { line ->
        def fields = line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
        [headers, fields].transpose().collectEntries { k, v -> v ? [(k): v] : [:] }
    }.findAll { it }
}

def check_reference_contig(fasta, contig) {
    if (!contig || contig == null || contig == false){
        return true
    }
    def fai_file = new File("${fasta}.fai")
    if (!fai_file.exists()) {
        error "ERROR: FASTA index file not found: ${fai_file}"
    }
    def found = false
    fai_file.eachLine { line ->
        def current_contig = line.split('\t')[0]        
        if (current_contig == contig) {
            found = true
        }
    }
    if (!found) {
        error "ERROR: Contig '${contig}' not found in reference index: ${fai_file}"
    }
    return true
}

def check_dragen_hash_contig(dragen_ref, contig) {
    // if no contig is passed, then continue
    if (!contig || contig == null || contig == false){
        return true
    }
    def cfg_file = new File("${dragen_ref}/hash_table.cfg")
    if (!cfg_file.exists()) {
        error "ERROR: DRAGEN hash table config not found: ${cfg_file}"
    }
    def found = false
    def target_string = "'${contig}'"
    cfg_file.eachLine { line ->
        if (line.contains(target_string)) {
            found = true
        }
    }
    if (!found) {
        error "ERROR: Contig '${contig}' not found in DRAGEN reference config: ${cfg_file}"
    }
    return true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// input mastersheet
ch_mastersheet = params.input ?
    Channel.fromPath("${params.input}", checkIfExists: true) :
    Channel.empty()

// DRAGEN reference directory
ch_reference_dir = params.refdir && check_dragen_hash_contig(params.refdir, params.transgene_name)
    ? Channel.fromPath(params.refdir, type: 'dir', checkIfExists: true).collect()
    : Channel.empty()

// FastA reference
ch_fasta_reference = params.fasta && check_reference_contig(params.fasta, params.transgene_name)
    ? Channel.fromPath("${params.fasta}*", checkIfExists: true).collect()
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

// SNV systematic noise BED file
ch_snv_noisefile = params.snv_noisefile
    ? Channel.fromPath(params.snv_noisefile, checkIfExists: true).collect()
    : []

// SV systematic noise BED file
ch_sv_noisefile = params.sv_noisefile
    ? Channel.fromPath(params.sv_noisefile, checkIfExists: true).collect()
    : []

// CRAM reference file
ch_cram_reference = params.cram_reference
    ? Channel.fromPath("${params.cram_reference}*", checkIfExists: true).collect()
    : []

// Gene regions / target file for analysis
ch_param_target_file = params.target_file ?
    Channel.fromPath("${params.target_file}", checkIfExists: true).first() 
    : []

// Nirvana path
ch_nirvana_path = params.nirvana_path
    ? Channel.fromPath("${params.nirvana_path}", checkIfExists: true).collect()
    : []

/*
~~~~~~~~~~~~~~~~~~
MultiQC parameters
~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

ch_multiqc_custom_config = params.multiqc_config 
    ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) 
    : Channel.empty()

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

    main:
    ch_versions = Channel.empty()
    ch_demux_output = Channel.empty()
    ch_alignment_samples = Channel.empty()
    ch_target_files = Channel.empty()   
    ch_dragen_output = Channel.empty()
    ch_dragen_usage = Channel.empty()

    //
    // MODULE: Parse input samplesheet to format samples for processing.
    //         Output of this process are csv files for samples that need to be aligned
    //         and samples that need to be analyzed
    //
    PARSE_INPUT_SAMPLESHEET (
        ch_input_samplesheet
    )
    ch_versions = ch_versions.mix(PARSE_INPUT_SAMPLESHEET.out.versions)

    // Get dragen outputs and add target files to analyze
    ch_dragen_output = ch_dragen_output.mix(
        PARSE_INPUT_SAMPLESHEET.out.samples_to_analyze
                .map{ generateMetaFromCsv(it) }
                .flatten()
                .filter{ it.dragen_path }
                .map{ meta -> 
                    def dragenfiles = file("${meta.dragen_path}/*")
                    // Use global target_file param if provided, otherwise use meta.target_file from CSV
                    def targetfile = params.target_file ? 
                        file(params.target_file, checkIfExists: true) : 
                        (meta.target_file ? file(meta.target_file, checkIfExists: true) : [])
                        [ meta, dragenfiles, targetfile ] 
                }
    )

    // Now get samples to align
    ch_sample_meta = PARSE_INPUT_SAMPLESHEET.out.samples_to_align
        .map{ generateMetaFromCsv(it) }
        .flatten()
        
    // get editing target file as separate channel
    ch_target_files = ch_target_files.mix(
        ch_sample_meta
            .filter{ it.sample_type == "tumor" }
            .map { meta -> 
                def targetfile = params.target_file ? 
                        file(params.target_file, checkIfExists: true) : 
                        (meta.target_file ? file(meta.target_file, checkIfExists: true) : [])
                if (targetfile!=[]){ 
                    [ meta.id, targetfile ]
                } else {
                    error "NO Target file provided."
                }
            }
    )

    // Get reads/fastqlists to align
    GATHER_ALIGNMENT_SAMPLES (
        ch_sample_meta,
        ch_demux_output.ifEmpty([]),
        ch_cram_reference
    )
    ch_versions = ch_versions.mix(GATHER_ALIGNMENT_SAMPLES.out.versions)

    PREPARE_SOMATIC_FASTQS(GATHER_ALIGNMENT_SAMPLES.out.samples)

    // Make hotspot VCF file, which includes gene bed file and 
    // Nominated off-target sites in the target_file. These are keyed by id
    MAKE_HOTSPOT_VCF(
        ch_target_files,
        ch_hotspot_bed,
        ch_fasta_reference
    )
    ch_versions = ch_versions.mix(MAKE_HOTSPOT_VCF.out.versions)
    ch_hotspot_vcf = MAKE_HOTSPOT_VCF.out.hotspot_vcf

    // Join alignment samples with hotspot VCF
    ch_alignment_samples = PREPARE_SOMATIC_FASTQS.out.samples
        .map{ meta, reads, fastqlist -> [ meta.id, meta, reads, fastqlist ] }
        .join(
            MAKE_HOTSPOT_VCF.out.hotspot_vcf
        )
        .map{ id, meta, reads, fastqlist, hotspot_vcf -> [ meta, reads, fastqlist, hotspot_vcf ] }
    
    ch_alignment_samples.dump(tag:'alignment_samples',pretty:true)

    if (params.run_alignment) {
        DRAGEN_SCGE (
            ch_alignment_samples,
            ch_intermediate_dir,
            ch_reference_dir,
            ch_adapter1_file,
            ch_adapter2_file,
            ch_sv_noisefile,
            ch_snv_noisefile,
            ch_nirvana_path
        )
        ch_versions     = ch_versions.mix(DRAGEN_SCGE.out.versions)
        ch_dragen_usage = ch_dragen_usage.mix(DRAGEN_SCGE.out.usage)
        ch_dragen_output = ch_dragen_output.mix(
                DRAGEN_SCGE.out.dragen_output
                .map { meta, dragenfiles -> [ meta.id, meta, dragenfiles ] }
                .join(ch_target_files)
                .map { id, meta, dragenfiles, targetfile -> [ meta, dragenfiles, targetfile ] }
        )
    }

    ch_dragen_output.dump(tag:'dragen_output',pretty:true)

    if (params.run_analysis) {

        SCGE_ANALYSIS(
            ch_dragen_output,
            ch_fasta_reference
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
