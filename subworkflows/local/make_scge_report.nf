// make scge report

include { MAKE_QUARTO_REPORT } from '../../modules/local/make_quarto_report.nf'

workflow MAKE_SCGE_REPORT {
    take:
    ch_report_json

    main:
    MAKE_QUARTO_REPORT(ch_report_json, params.scge_report_qmd)
    
}