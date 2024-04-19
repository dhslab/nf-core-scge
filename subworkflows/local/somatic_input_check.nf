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
    .map { meta -> 
        if (meta.id == params.tumorid){
            meta.sample_type = 'tumor'
            meta.uid = params.id
        } else if (meta.id == params.normalid){
            meta.sample_type = 'normal'
            meta.uid = params.id
        }   
        return meta 
    }
    .filter { it.sample_type != null || it.dragen_path != null }
    .set { ch_mastersheet }

    ch_mastersheet.dump()
    
    ch_mastersheet
    .map { meta -> 
        if (meta.read1 != null && meta.read2 != null){
            def new_meta = [:]
            new_meta['id'] = meta.uid
            new_meta['assay'] = meta.assay
            def sample_id = ""
            if (meta.sample_id != null){
                sample_id = '.' + meta.sample_id
            }
            def rgid = meta.flowcell + sample_id + '.' + meta.i7index + '.' + meta.i5index + '.' + meta.lane 
            def rglb = meta.id + '.' + meta.i7index + '.' + meta.i5index
            tumor_id = ""
            normal_id = ""
            if (meta.sample_type == 'tumor'){
                tumor_id = meta.id
            } else if (meta.sample_type == 'normal'){
                normal_id = meta.id
            }

            [ new_meta, tumor_id, normal_id, [ rgid, meta.id, rglb, meta.lane, file(meta.read1), file(meta.read2) ] ]
        }
    }
    .groupTuple()
    .map { meta, tumor, normal, fqlist -> 
            new_meta = meta.subMap('id', 'assay')
            new_meta['tumor'] = tumor.findAll { it != '' }.unique()[0]
            new_meta['normal'] = normal.findAll { it != '' }.unique()[0]

            def fileList = ['RGID,RGSM,RGLB,Lane,Read1File,Read2File']
            def read1 = []
            def read2 = []

            // Create data rows
            for (int i = 0; i < fqlist.size(); i++) {
                def row = fqlist[i]
                read1 << file(row[4])
                read2 << file(row[5])
                fileList << [ row[0], row[1], row[2], row[3], row[4].toString().split('/')[-1], row[5].toString().split('/')[-1] ].join(',')
            }
            return [ new_meta, fileList.join('\n'), read1, read2 ]
    }
    .filter { it[0].tumor != "" && it[0].normal != "" }
    .set { ch_fastqs }

    ch_fastqs.dump()

    // Put read1 and read2 files into separate channels.
    ch_fastqs
    .map { meta, fqlist, read1, read2 -> [ meta, read1 ] }
    .transpose()
    .set { ch_read1 }

    ch_fastqs
    .map { meta, fqlist, read1, read2 -> [ meta, read2 ] }
    .transpose()
    .set { ch_read2 }

    // Make fastq_list file from string.
    ch_fastqs
    .map { meta, fqlist, read1, read2 -> [ meta, fqlist ] } | MAKE_FASTQLIST

    // Concatenate read1, read2, and fastqlist channels and group by meta.
    MAKE_FASTQLIST.out.fastq_list
    .concat(ch_read1,ch_read2)
    .transpose()
    .groupTuple()
    .map { meta, files -> 
        return [ meta, 'fastq', files ]
    }
    .set { ch_input_data }

    // Organize cram files into a channel.
    // this saves the cram file names as values and the files as a list of files
    ch_mastersheet
    .map { meta -> 
        if (meta.cram != null){
            def new_meta = [:]
            new_meta['id'] = meta.uid
            new_meta['assay'] = meta.assay
            if (meta.sample_type == 'tumor'){
                return [ new_meta, file(meta.cram).getName(), '', [ file(meta.cram), file(meta.cram + '.crai') ] ]
            } else if (meta.sample_type == 'normal'){
                return [ new_meta, '', file(meta.cram).getName(), [ file(meta.cram), file(meta.cram + '.crai') ] ]
            }
        }
    }
    .groupTuple()
    .map { meta, tumor, normal, crams ->
        meta['tumor'] = tumor.findAll { it != '' }.unique()[0]
        meta['normal'] = normal.findAll { it != '' }.unique()[0]
        return [ meta, 'cram', crams ]
    }
    .filter { it.tumor != "" && it.normal != "" }
    .set { ch_cram }

    ch_input_data = ch_input_data.mix(ch_cram)

    // Organize bam files into a channel.
    ch_mastersheet
    .map { meta -> 
        if (meta.bam != null){
            def new_meta = [:]
            new_meta['id'] = meta.uid
            new_meta['assay'] = meta.assay
            if (meta.sample_type == 'tumor'){
                return [ new_meta, file(meta.bam).getName(), '', [ file(meta.bam), file(meta.bam + '.crai') ] ]
            } else if (meta.sample_type == 'normal'){
                return [ new_meta, '', file(meta.bam).getName(), [ file(meta.bam), file(meta.bam + '.crai') ] ]
            }
        }
    }
    .groupTuple()
    .map { meta, tumor, normal, bams ->
        meta['tumor'] = tumor.findAll { it != '' }.unique()[0]
        meta['normal'] = normal.findAll { it != '' }.unique()[0]
        return [ meta, 'bam', bams ]
    }
    .filter { it.tumor != "" && it.normal != "" }
    .set { ch_bam }

    ch_input_data = ch_input_data.mix(ch_bam)

    ch_mastersheet
    .map { meta -> 
        if (meta.dragen_path != null){
            def new_meta = meta.subMap('id','assay','dragen_path')
            return [ new_meta, file(meta.dragen_path).listFiles() ]
        }
    }
    .set { ch_dragen_outputs }

    emit:
    dragen_outputs = ch_dragen_outputs
    input_data = ch_input_data
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
