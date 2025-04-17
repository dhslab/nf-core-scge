// make circos plots
// make scge report

include { MAKE_CIRCOS_PLOT   } from '../../modules/local/make_circos_plot.nf'
include { MAKE_QUARTO_REPORT } from '../../modules/local/make_quarto_report.nf'

workflow MAKE_SCGE_REPORT {
    take:
    ch_indels
    ch_transgene 

    main:
    get_circos_input(ch_transgene)

    MAKE_CIRCOS_PLOT(get_circos_input.out.circos_input)

    quarto_report_input = MAKE_CIRCOS_PLOT.out.circos_png.join(ch_indels).dump()
    MAKE_QUARTO_REPORT(quarto_report_input, params.scge_report_qmd)
    
}

process get_circos_input {
    label 'process_single'
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(meta), path(transgene_file)

    output:
    tuple val(meta), path("${meta.id}.circos_input.tsv"), emit: circos_input

    script:
    """
    transform_transgene.py --input $transgene_file --output ${meta.id}.circos_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}