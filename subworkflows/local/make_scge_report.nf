// make scge report

include { RENDER_SCGE_REPORT } from '../../modules/local/render_scge_report.nf'

workflow MAKE_SCGE_REPORT {
    take:
    ch_report_input

    main:
    RENDER_SCGE_REPORT(ch_report_input, params.scge_report_qmd)
    
}