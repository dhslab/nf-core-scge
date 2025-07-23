#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Include the reporting subworkflow from your existing pipeline
include { MAKE_SCGE_REPORT } from './subworkflows/local/make_scge_report.nf'

params.samplesheet = 'custom_samplesheet.csv'
params.outdir = 'report_only'

workflow {
    // Read the custom samplesheet
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id]
            def transgene_file = file(row.transgene_file)
            def indels_file = file(row.indels_file)
            return [meta, transgene_file, indels_file]
        }
        .set { input_ch }

    // Split into separate channels for the reporting workflow
    input_ch
        .map { meta, transgene, indels -> [meta, indels] }
        .set { indels_ch }

    input_ch
        .map { meta, transgene, indels -> [meta, transgene] }
        .set { transgene_ch }

    // Run only the reporting step
    MAKE_SCGE_REPORT(indels_ch, transgene_ch)
} 