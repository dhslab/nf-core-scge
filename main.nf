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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
