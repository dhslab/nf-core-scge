#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dhslab/scge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/dhslab/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCGE                    } from './workflows/scge'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scge_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scge_pipeline'
include { VALIDATE_PARAMS         } from './subworkflows/local/utils_nfcore_scge_pipeline'
include { OFFTARGET_WORKFLOW      } from './workflows/offtarget'
include { TRAIN_WORKFLOW          } from './workflows/train'
include { GENERATE_HOTSPOTS       } from './subworkflows/local/generate_hotspots'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Validate params against nextflow_schema.json, print the run summary,
    // and build the input channel. This is what makes --help and parameter validation
    // work; INPUT_CHECK is invoked inside it, so it must not also be called here.
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    ch_versions = ch_versions.mix(PIPELINE_INITIALISATION.out.versions)

    SCGE (
        PIPELINE_INITIALISATION.out.input
    )
    ch_versions = ch_versions.mix(SCGE.out.versions)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SCGE.out.multiqc_report
    )
}

//
// WORKFLOW: Unified CRISPR Off-Target Workflow.
// Run with:  nextflow run . -entry OFFTARGET -profile ris --input offtarget_samplesheet.csv --outdir ./results
// When -entry OFFTARGET is given, the default SCGE workflow above does not run.
//
workflow OFFTARGET {
    VALIDATE_PARAMS (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        "nextflow run ${workflow.manifest.name} -entry OFFTARGET -profile ris2,apptainer --input offtarget_samplesheet.csv --outdir <OUTDIR>"
    )
    OFFTARGET_WORKFLOW()
}

//
// WORKFLOW: Off-Target Model Trainer (offline retrain loop).
// Run with:  nextflow run . -entry TRAIN -profile ris2,apptainer --input training.tsv --outdir ./results
// Fits a new wgs_shape_model.pkl from a BUILD_TRAINING_TABLE output; deploy it via
// -entry OFFTARGET --offtarget_shape_model <new.pkl>.
//
workflow TRAIN {
    VALIDATE_PARAMS (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        "nextflow run ${workflow.manifest.name} -entry TRAIN -profile ris2,apptainer --input training.tsv --outdir <OUTDIR>"
    )
    TRAIN_WORKFLOW()
}

//
// WORKFLOW: Auto-Hotspot Finder — gRNA -> predicted off-target sites (target_file).
// Run with:  nextflow run . -entry HOTSPOTS -profile ris2,apptainer --input grna_samplesheet.csv --outdir ./results
// Emits, per guide, <guide>.targets.csv (readable) and <guide>.targets.vcf (feed as the
// OFFTARGET arm's per-row target_file). Cas-OFFinder by default; add --run_crisprme with a
// prebuilt --crisprme_index_dir to also include CRISPRme.
//
workflow HOTSPOTS {
    VALIDATE_PARAMS (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        "nextflow run ${workflow.manifest.name} -entry HOTSPOTS -profile ris2,apptainer --input grna_samplesheet.csv --outdir <OUTDIR>"
    )

    if (!params.input) { error "HOTSPOTS: provide --input grna_samplesheet.csv (columns: guide_id,spacer[,pam,idt])" }

    ch_guides = Channel.fromPath(params.input, checkIfExists: true)
        | splitCsv(header: true)
        | map { row ->
            if (!row.guide_id?.trim() || !row.spacer?.trim()) {
                error "HOTSPOTS: every --input row needs non-empty 'guide_id' and 'spacer' (got: ${row})"
            }
            def pam = (row.pam?.trim()) ?: params.offtarget_pam
            def idt = (row.idt?.trim()) ? file(row.idt.trim(), checkIfExists: true)
                                        : file("${projectDir}/assets/NO_IDT")
            tuple([id: row.guide_id.trim()], row.spacer.trim().toUpperCase(), pam, idt)
        }

    GENERATE_HOTSPOTS(ch_guides)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
