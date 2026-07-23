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
include { INPUT_CHECK             } from './subworkflows/local/input_check'
include { OFFTARGET_WORKFLOW      } from './workflows/offtarget'
include { TRAIN_WORKFLOW          } from './workflows/train'

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

    INPUT_CHECK(params.input)
        .set { ch_input }

    ch_versions = ch_versions.mix(ch_input.versions)

    SCGE (
        ch_input.input
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
    OFFTARGET_WORKFLOW()
}

//
// WORKFLOW: Off-Target Model Trainer (offline retrain loop).
// Run with:  nextflow run . -entry TRAIN -profile ris2,apptainer --input training.tsv --outdir ./results
// Fits a new wgs_shape_model.pkl from a BUILD_TRAINING_TABLE output; deploy it via
// -entry OFFTARGET --offtarget_shape_model <new.pkl>.
//
workflow TRAIN {
    TRAIN_WORKFLOW()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
