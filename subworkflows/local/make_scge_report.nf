// make scge report

include { RENDER_SCGE_REPORT } from '../../modules/local/render_scge_report.nf'
include { GENERATE_EXCEL_REPORT } from '../../modules/local/generate_excel_report.nf'

workflow MAKE_SCGE_REPORT {
    take:
    ch_report_input

    main:
    RENDER_SCGE_REPORT(ch_report_input, params.scge_report_qmd)

    // Join inputs:
    // ch_report_input: [meta, json, cna, baf, circos]
    // RENDER.out.indel_plot: [meta, indel_png]
    // RENDER.out.off_targets_plot: [meta, off_target_png]
    
    ch_excel_input = ch_report_input
        .join(RENDER_SCGE_REPORT.out.indel_plot)
        .join(RENDER_SCGE_REPORT.out.off_targets_plot)
        // Result: [meta, json, cna, baf, circos, indel, off_target]
        
    GENERATE_EXCEL_REPORT(ch_excel_input)

    emit:
    scge_report = RENDER_SCGE_REPORT.out.scge_report
    excel_report = GENERATE_EXCEL_REPORT.out.excel_report
}