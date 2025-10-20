/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_MGI_SAMPLEMAP                          } from '../../modules/local/convert_mgi_samplemap.nf'
include { CREATE_FASTQ_LIST as CREATE_TUMOR_FASTQ_LIST   } from '../../modules/local/create_fastq_list.nf'
include { CREATE_FASTQ_LIST as CREATE_NORMAL_FASTQ_LIST  } from '../../modules/local/create_fastq_list.nf'

// Parse CSV files from 'ch_samples_to_align' for each sample into 'meta'
def generateMetaFromCsv(csv_file) {
    def lines = csv_file.text.readLines()
    def headers = lines[0].split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')

    return lines.drop(1).collect { line ->
        def fields = line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
        [headers, fields].transpose().collectEntries { k, v -> v ? [(k): v] : [:] }
    }.findAll { it }
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

/*
========================================================================================
    SUBWORKFLOW TO GATHER ALIGNMENT SAMPLES
========================================================================================
*/

workflow GATHER_ALIGNMENT_SAMPLES {

    take:
    ch_samples_to_align  // channel: [ path(file) ]
    ch_cram_reference    // channel: [ path(file) ]

    main:
    ch_versions          = Channel.empty()
    ch_fastqlists        = Channel.empty()
    ch_gathered_fastqs   = Channel.empty() // channel: [ val(meta), path(read1), path(read2), path(runinfo.xml) ]

    // This is the output of this subworkflow and is the format for DRAGEN alignment
    ch_gathered_samples  = Channel.empty() // [ val(sample_info), path(tumor_reads), path(tumor_fastq_list), path(normal_reads), path(normal_fastq_list) ]

    // Channel of meta data for alignment samples
    ch_sample_alignment_meta = ch_samples_to_align
                                .map{ generateMetaFromCsv(it) }
                                .flatten()

    ch_sample_alignment_meta.dump(tag: 'gather_alignment_samples:ch_sample_alignment_meta', pretty: true)

    //
    // Get CRAM/BAM files that can be directly realigned
    //
    ch_gathered_samples = ch_gathered_samples.mix(
                            ch_sample_alignment_meta
                                .map { meta -> 
                                    if (meta.edited_cram && meta.control_cram){
                                        return [ meta, [], [], file("${meta.edited_cram}*"), file("${meta.control_cram}*") ]
                                    } else if (meta.edited_bam && meta.control_bam){
                                        return [ meta, [], [], file("${meta.edited_bam}*"), file("${meta.control_bam}*") ]
                                    }
                                }
                            )

    //
    // Convert MGI samplemap to fastq_list
    //
    CONVERT_MGI_SAMPLEMAP (
        ch_sample_alignment_meta
            .filter{ it.mgi_samplemap && it.edited_id && it.control_id }
            .map { meta -> [ meta, file(meta.mgi_samplemap, checkIfExists: true) ] }
    )

    ch_fastqlists = ch_fastqlists.mix(
        CONVERT_MGI_SAMPLEMAP.out.fastq_list
            .map { meta, fastq_list -> 
                meta.fastq_list = fastq_list
                return meta
        },
        ch_sample_alignment_meta.filter{ it.fastq_list && it.edited_id && it.control_id }
    )

    ch_fastqlists.dump(tag: 'gather_alignment_samples:ch_fastqlists', pretty: true)

    //
    // Collect tumor fastq_list, and runinfo.
    //
    ch_gathered_tumor_fastqs = ch_fastqlists
        .flatMap{ meta -> 
            def requiredColumns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
            def fastq_list = file(meta.fastq_list, checkIfExists: true)                    
            def data = parseFastqList(fastq_list)
            data = data.findAll{ it.RGSM == meta.edited_id }                    
            data.collect{
                if (!it.keySet().containsAll(requiredColumns)) {
                    error("Missing required columns in input FastQ list!")
                }
                def R1 = file(it['Read1File'], checkIfExists: true)
                def R2 = file(it['Read2File'], checkIfExists: true)                    
                [ meta, meta.edited_id, R1, R2 ]
            }
        }
        .filter{ it!= [] }
    
    ch_gathered_tumor_fastqs.dump(tag: 'gather_alignment_samples:ch_gathered_tumor_fastqs', pretty: true)

    //
    // MODULE: Create normal fastq_list with local/staged fastq paths and appropriate metadata
    //
    CREATE_TUMOR_FASTQ_LIST (
        ch_gathered_tumor_fastqs
    )
    ch_versions = ch_versions.mix(CREATE_TUMOR_FASTQ_LIST.out.versions)

    //
    // Collect normal fastq_list, and runinfo.
    //
    ch_gathered_normal_fastqs = ch_fastqlists
        .flatMap{ meta -> 
            def requiredColumns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
            def fastq_list = file(meta.fastq_list, checkIfExists: true)                    
            def data = parseFastqList(fastq_list)
            data = data.findAll{ it.RGSM == meta.control_id }                    
            data.collect{
                if (!it.keySet().containsAll(requiredColumns)) {
                    error("Missing required columns in input FastQ list!")
                }
                def R1 = file(it['Read1File'], checkIfExists: true)
                def R2 = file(it['Read2File'], checkIfExists: true)                    
                [ meta, meta.control_id, R1, R2 ]
            }
        }
        .filter{ it!= [] }
    
    ch_gathered_normal_fastqs.dump(tag: 'gather_alignment_samples:ch_gathered_normal_fastqs', pretty: true)

    //
    // MODULE: Create normal fastq_list with local/staged fastq paths and appropriate metadata
    //
    CREATE_NORMAL_FASTQ_LIST (
        ch_gathered_normal_fastqs
    )
    ch_versions = ch_versions.mix(CREATE_NORMAL_FASTQ_LIST.out.versions)

    //
    // SUBWORKFLOW: Verify and parse tumor fastq_list files and staged reads
    //
    ch_gathered_tumor_samples = ch_sample_alignment_meta
        .join( 
            ch_gathered_tumor_fastqs
                .filter{ it != [] }
                .map{ meta, id, read1, read2 -> [ meta, [ read1, read2 ] ] }
                .groupTuple()
                .map { meta, reads -> [ meta, reads.flatten() ] }
        )
        .combine(
            CREATE_TUMOR_FASTQ_LIST.out.fastq_list  // NOTE (DHS): This creates a fastq_list will all samples in this run that is reused.
                .map{ meta, fastq_list ->
                    def data = parseFastqList(fastq_list)
                    data.each{
                        if (it) {
                            it['Read1File'] = "tumor_fastq_files/${it['Read1File'].split('/')[-1]}"
                            it['Read2File'] = "tumor_fastq_files/${it['Read2File'].split('/')[-1]}"
                        }
                    }

                    if (data) {
                        def header = data[0].keySet().join(',')
                        def content = data.collect { it.values().join(',') }.join('\n')

                        [ header, content ]
                    } else {
                        []
                    }
                }
                .flatten()
                .collectFile(
                    name   : "tumor_fastq_list.csv",
                    newLine: true,
                    sort   : 'index'
                )
        )
        .map{ meta, reads, fastq_list -> [ meta, reads.flatten(), fastq_list ] }

    //
    // SUBWORKFLOW: Verify and parse normal fastq_list files and staged reads
    //
    ch_gathered_normal_samples = ch_sample_alignment_meta
        .join( 
            ch_gathered_normal_fastqs
                .filter{ it != [] }
                .map{ meta, id, read1, read2 -> [ meta, [ read1, read2 ] ] }
                .groupTuple()
                .map { meta, reads -> [ meta, reads.flatten() ] }
        )
        .combine(
            CREATE_NORMAL_FASTQ_LIST.out.fastq_list  // NOTE (DHS): This creates a fastq_list will all samples in this run that is reused.
                .map{ meta, fastq_list ->
                    def data = parseFastqList(fastq_list)
                    data.each{
                        if (it) {
                            it['Read1File'] = "normal_fastq_files/${it['Read1File'].split('/')[-1]}"
                            it['Read2File'] = "normal_fastq_files/${it['Read2File'].split('/')[-1]}"
                        }
                    }

                    if (data) {
                        def header = data[0].keySet().join(',')
                        def content = data.collect { it.values().join(',') }.join('\n')

                        [ header, content ]
                    } else {
                        []
                    }
                }
                .flatten()
                .collectFile(
                    name   : "normal_fastq_list.csv",
                    newLine: true,
                    sort   : 'index'
                )
        )
        .map{ meta, reads, fastq_list -> [ meta, reads.flatten(), fastq_list ] }

    ch_gathered_samples = ch_gathered_tumor_samples.join( ch_gathered_normal_samples )

    ch_gathered_samples.dump(tag: 'gather_alignment_samples:ch_gathered_samples', pretty: true)

    emit:
    samples  = ch_gathered_samples  // channel: [ val(meta), path(reads), path(fastq_list), path(alignment_file) ]
    versions = ch_versions           // channel: [ path(file) ]

}
