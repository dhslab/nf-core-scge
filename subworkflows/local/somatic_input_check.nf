include { SAMPLESHEET_CHECK              } from '../../modules/local/samplesheet_check.nf'
include { MAKE_FASTQLIST                 } from '../../modules/local/make_fastqlist.nf'

workflow SOMATIC_INPUT_CHECK {
    take:
    master_samplesheet
    data_path

    main:

    ch_mastersheet        = Channel.empty()
    ch_input_data         = Channel.empty()
    ch_dragen_outputs     = Channel.empty()

    // Runs a python script that parses the sample sheet and adds key metadata, 
    // including index sequences, flowcell, and lane. If fastq_list.csv files are passed,
    // these are parsed also and read1/read2 pairs are returned. 
    // Also, if tumor/normal and case ids are supplied, mark the 
    // rows in the sheet as tumor or normal and assign a case id.
    // The output is then channelified.
    SAMPLESHEET_CHECK ( master_samplesheet, data_path )
    .csv
    .splitCsv ( header:true, sep:',' )
    .map { create_master_samplesheet(it) }
    .map { row -> [ row.uid ?: row.id, row ] }
    .groupTuple()
    .map { key, rows ->
        def tumor_rows = rows.findAll { it.sample_type == 'tumor' }
        def normal_rows = rows.findAll { it.sample_type == 'normal' }

        // We expect one tumor and one normal, but can handle multiple rows (e.g. for multiple fastqs)
        if (!tumor_rows.isEmpty() && !normal_rows.isEmpty()) {
            def tumor_meta = tumor_rows.first().subMap('id', 'assay', 'dragen_path', 'uid')
            def normal_meta = normal_rows.first().subMap('id', 'assay', 'dragen_path', 'uid')
            tumor_meta.normal_id = normal_meta.id // Keep track of the paired normal
            
            // Collect all files from the respective DRAGEN output paths
            def tumor_files = tumor_rows.collect { it.dragen_path }.unique().findAll { it != null }.collect { file(it).listFiles().collect{it.toString()} }.flatten()
            def normal_files = normal_rows.collect { it.dragen_path }.unique().findAll { it != null }.collect { file(it).listFiles().collect{it.toString()} }.flatten()
            
            def all_files = []
            all_files.addAll(tumor_files)
            all_files.addAll(normal_files)

            return [ tumor_meta, all_files ]
        } else {
            // Handle cases with only precomputed dragen_path and no tumor/normal pairs
            def meta = rows.first().subMap('id', 'assay', 'dragen_path')
            def files = rows.collect{it.dragen_path}.unique().findAll { it != null }.collect{ file(it).listFiles().collect{it.toString()} }.flatten()
            return [ meta, files ]
        }
    }
    .set { ch_dragen_outputs }


    // The following logic for fastqs, crams, and bams is for DRAGEN runs.
    // Since we are running with --run_dragen false, we can simplify this subworkflow
    // to primarily focus on emitting the ch_dragen_outputs channel correctly.
    // The ch_input_data and ch_hotspots are still needed for analysis mode.
    ch_mastersheet = SAMPLESHEET_CHECK.out.csv
        .splitCsv( header:true, sep:',' )
        .map { create_master_samplesheet(it) }
        .filter { it.sample_type == 'tumor' || it.dragen_path != null } // Only process tumors or samples with dragen_path


    ch_input_data = ch_mastersheet.map{ row ->
        def meta = row.subMap('id', 'assay', 'uid', 'sample_type')
        // When not running dragen, we don't need files, just the metadata
        [meta, []]
    }

    ch_hotspots = ch_mastersheet
        .map{ row ->
            def hs = params.hotspot_csv ?: row.hotspot_file ?: "$projectDir/assets/NO_FILE.csv"
            def key = row.uid ?: row.id
            [key, hs]
        }
        .unique()


    // This join is probably not needed in analysis only mode, but we keep it for now.
    ch_input_data = ch_input_data
        .map{ meta, files ->
            def key = meta.uid ?: meta.id
            [key, [meta, files]]
        }
        .join(ch_hotspots, by: 0)
        .map{ id, info, hotspot ->
            // The structure for ch_input_data is [meta, file_list]
            // In analysis mode, the file list is empty, but we add the hotspot file
            [info[0], [info[1], hotspot]]
        }


    emit:
    dragen_outputs = ch_dragen_outputs
    input_data = ch_input_data
    hotspots = ch_hotspots
}

def create_master_samplesheet(LinkedHashMap row) {

    def meta = [:]
    meta.id             = row.id
    meta.uid            = row.uid ?: null
    meta.sample_type    = row.sample_type ?: null
    meta.sample_id      = row.sample_id ?: null
    meta.assay          = row.assay ?: null
    meta.i7index        = row.i7index ?: null
    meta.i5index        = row.i5index ?: null
    meta.flowcell       = row.flowcell ?: null
    meta.lane           = row.lane ?: null
    meta.dragen_path    = null
    meta.read1          = null
    meta.read2          = null
    meta.cram           = null
    meta.bam            = null
    meta.hotspot_file   = row.hotspot_file ?: null

    if (row.fastq_list) {
        if (!file(row.fastq_list).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> fastq_list does not exist!\n${row.fastq_list}"
        }
        meta.fastq_list = file(row.fastq_list)
    }

    if (row.dragen_path) {
        if (!file(row.dragen_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> dragen_path does not exist!\n${row.dragen_path}"
        }
        meta.dragen_path = file(row.dragen_path)
    }

    if (row.read1) {
        if (!file(row.read1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 file does not exist!\n${row.read1}"
        }
        meta.read1 = file(row.read1)
    }
    if (row.read2) {
        if (!file(row.read2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 file does not exist!\n${row.read2}"
        }
        meta.read2 = file(row.read2)
    }

    if (row.cram) {
        if (!file(row.cram).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Cram file does not exist!\n${row.cram}"
        }
        meta.cram = file(row.cram)
    }

    if (row.bam) {
        if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
        }
        meta.bam = file(row.bam)
    }
    return meta
}
