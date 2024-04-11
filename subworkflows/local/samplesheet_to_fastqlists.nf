//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow SAMPLESHEET_TO_FASTQLISTS {
    take:
    input // file: /path/to/samplesheet.csv

    main:
    make_fastqlists(input)

    normal_fastq_info = make_fastqlists.out.normal_fastqlist.splitCsv(header: true)
        .map{row -> def meta = row.subMap('RGID', 'RGLB', 'RGSM')
        [row.'RGSM', meta]
        }

    tumor_fastq_info = make_fastqlists.out.tumor_fastqlist.splitCsv(header: true)
        .map{row -> def meta = row.subMap('RGID', 'RGLB', 'RGSM')
        [row.'RGSM', meta]
        }

    normals_id = Channel.fromPath(input).splitCsv(header: true).filter{ it.type == 'normal' }
        .map{row -> def meta = row.subMap('sample')
        [row.id, meta]
        }
    
    tumors_id = Channel.fromPath(input).splitCsv(header: true).filter{ it.type == 'tumor' }
        .map{row -> def meta = row.subMap('sample')
        [row.id, meta]
        }
    
    // order normals and tumors
    normals = normals_id.join(normal_fastq_info).combine(make_fastqlists.out.normal_fastqlist)
        .map{id, sample, fastq_info, fastqlist ->
        [sample.sample, fastq_info, fastqlist]
        }
        
    tumors = tumors_id.join(tumor_fastq_info).combine(make_fastqlists.out.tumor_fastqlist)
        .map{id, sample, fastq_info, fastqlist ->
        [sample.sample, fastq_info, fastqlist]
        }
        
    dragen_input = normals.combine(tumors, by: 0).map{sample, normal_meta, normal_fastqlist, tumor_meta, tumor_fastqlist -> [normal_meta, normal_fastqlist, tumor_meta, tumor_fastqlist]}
    dragen_input.dump()

    emit:
    dragen_input
    versions = make_fastqlists.out.versions

}

process make_fastqlists {
    label 'process_single'
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path input

    output:
    path "tumor_fastqlist.csv", emit: tumor_fastqlist
    path "normal_fastqlist.csv", emit: normal_fastqlist
    path('versions.yml'),          emit: versions

    script:
    """
    make_fastqlists.py $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}