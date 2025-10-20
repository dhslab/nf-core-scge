/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_XLSX_TO_CSV } from '../../modules/local/convert_xlsx_to_csv'

/*
========================================================================================
    SUBWORKFLOW TO CHECK INPUTS
========================================================================================
*/

workflow INPUT_CHECK {

    take:
    input          //  string: Path to input samplesheet

    main:
    ch_versions = Channel.empty()

    /*
    ================================================================================
                        Process input MGI samplesheet
    ================================================================================
    */

    if (input) {
        // Verify input samplesheet has a file extension in [xlsx,csv,tsv]
        if (hasExtension(input, 'xlsx')) {
            CONVERT_XLSX_TO_CSV (
                Channel.fromPath(input, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(CONVERT_XLSX_TO_CSV.out.versions)

            ch_input = CONVERT_XLSX_TO_CSV.out.csv

        } else if (hasExtension(input, 'csv') || hasExtension(input, 'tsv')) {
            ch_input = Channel.fromPath(input, checkIfExists: true)
        } else {
            error("Input samplesheet input does not end in `.{xlsx,csv,tsv}`!")
        }
    } else {
        ch_input = Channel.empty()
    }

    emit:
    input    = ch_input     // channel: [ path(file) ]
    versions = ch_versions  // channel: [ path(file) ]

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Get file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Parse FastQ list
def parseFastqList(file) {
    def separator = file.toString().endsWith("tsv") ? '\t' : ','
    def lines = file.readLines()
    def headers = lines.first().split(separator)
    lines.drop(1).collect{ line ->
        [headers, line.split(separator)].transpose().collectEntries{ it }
    }
}
