/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_VARIANTS             } from '../../modules/local/annotate_variants.nf'
include { VEP_TO_TSV as VARIANTS_TO_TSV } from '../../modules/local/vep_to_tsv.nf'
include { ANNOTATE_SV_VARIANTS          } from '../../modules/local/annotate_sv_variants.nf'
include { VEP_TO_TSV as SV_TO_TSV       } from '../../modules/local/vep_to_tsv.nf'
include { ANNOTATE_CNV_VARIANTS         } from '../../modules/local/annotate_cnv_variants.nf'
include { VEP_TO_TSV as CNV_TO_TSV      } from '../../modules/local/vep_to_tsv'
include { ANNOTATE_OFFTARGETS           } from '../../modules/local/annotate_offtargets.nf'
include { GET_INDELS                    } from '../../modules/local/get_indels.nf'
include { REVIEW_FILTER                 } from '../../modules/local/review_filter.nf'
include { REVIEW_FILTER_BND             } from '../../modules/local/review_filter_bnd.nf'
include { REVIEW_SNAPSHOTS             } from '../../modules/local/review_snapshots.nf'
include { REVIEW_BND_SNAPSHOTS         } from '../../modules/local/review_bnd_snapshots.nf'
include { GET_TRANSGENE_JUNCTIONS       } from '../../modules/local/get_transgene_junctions.nf'
include { TRANSGENE_TO_VCF              } from '../../modules/local/transgene_to_vcf'
include { ANNOTATE_TRANSGENE_JUNCTIONS } from '../../modules/local/annotate_transgene_junctions.nf'
include { TRANSFORM_TRANSGENE         } from '../../modules/local/transform_transgene.nf'
include { MAKE_CIRCOS_PLOT            } from '../../modules/local/make_circos_plot.nf'
include { REFORMAT_CNV_DATA           } from '../../modules/local/reformat_cnv_data.nf'
include { GENERATE_CNA_BAF_PLOTS      } from '../../modules/local/generate_cna_baf_plots.nf'
include { COMPILE_REPORT_JSON         } from '../../modules/local/compile_report_json.nf'
include { BND_FROM_INDELS_TO_VCF      } from '../../modules/local/bnd_from_indels_to_vcf.nf'
include { MAKE_SCGE_REPORT            } from './make_scge_report.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Vep cache
ch_vepcache = params.vepcache
    ? Channel.fromPath(params.vepcache, type: 'dir', checkIfExists: true).collect()
    : Channel.empty()

ch_cytobands = params.cytobands
    ? Channel.fromPath("${params.cytobands}*", checkIfExists: true).collect()
    : []

ch_crispr_model = params.crispr_model ?
    Channel.fromPath("${params.crispr_model}", checkIfExists: true).collect() 
    : []

ch_transgene_name = params.transgene_name && params.transgene_name != false && params.transgene_name != null
    ? Channel.value(params.transgene_name) : Channel.empty()


/*
========================================================================================
    SUBWORKFLOW TO ANALYZE DATA
========================================================================================
*/

workflow SCGE_ANALYSIS {

    take:
    ch_analysis_samples // channel: meta
    ch_fasta_reference // channel: path(fasta_reference)

    main:
    ch_versions = Channel.empty()
    ch_report_inputs = Channel.empty() // channel: meta, [ file1, file2, etc ]

    //
    // Main analysis
    //

    // Extract dragenfiles and target files from analysis samples
    ch_dragen_files = ch_analysis_samples
        .map { meta, dragenfiles, targetfile -> [ meta, dragenfiles ]}

    // Get coverage files for report
    ch_report_inputs = ch_report_inputs.mix(
            ch_dragen_files.map { meta, dragen_path ->
                def tumor_cov  = dragen_path.find { it.name.endsWith('.wgs_overall_mean_cov_tumor.csv') }
                def normal_cov = dragen_path.find { it.name.endsWith('.wgs_overall_mean_cov_normal.csv') }            
                return [ meta, [ file(tumor_cov,checkIfExists: true), file(normal_cov,checkIfExists: true) ] ]
            }
        )
 
    // Annotate small variants
    ANNOTATE_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    VARIANTS_TO_TSV (ANNOTATE_VARIANTS.out.vcf,Channel.value('vcf'))
    ch_report_inputs = ch_report_inputs.mix(VARIANTS_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(VARIANTS_TO_TSV.out.versions)

    // Annotate structural variants
    ANNOTATE_SV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    SV_TO_TSV (ANNOTATE_SV_VARIANTS.out.vcf,Channel.value('sv'))
    ch_report_inputs = ch_report_inputs.mix(SV_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(SV_TO_TSV.out.versions)

    // Annotate copy number variants
    ANNOTATE_CNV_VARIANTS (
        ch_dragen_files,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_versions = ch_versions.mix(ANNOTATE_CNV_VARIANTS.out.versions)

    CNV_TO_TSV (ANNOTATE_CNV_VARIANTS.out.vcf,Channel.value('cnv'))
    ch_report_inputs = ch_report_inputs.mix(CNV_TO_TSV.out.tsv)
    ch_versions = ch_versions.mix(CNV_TO_TSV.out.versions)

    // Annotate and analyze off-target sites
    ANNOTATE_OFFTARGETS(
        ch_analysis_samples.map{ meta, dragenfiles, targetfile -> [meta, targetfile] },
        ch_vepcache,
        ch_fasta_reference
    )
    ch_versions = ch_versions.mix(ANNOTATE_OFFTARGETS.out.versions)

    GET_INDELS(
        ch_dragen_files.join(ANNOTATE_OFFTARGETS.out.targetfile),
        ch_crispr_model,
        ch_fasta_reference
    )

    ch_report_inputs = ch_report_inputs.mix(GET_INDELS.out.indels_file)
    ch_versions = ch_versions.mix(GET_INDELS.out.versions)

    // Collapse the per-sample call tables into one short review queue. Every sample is staged
    // into a SINGLE call on purpose: the cross-guide fallback for rule 4 can only see recurrence
    // across whatever is passed together, so a per-sample invocation would silently disable it.
    // With review_noise_model set that fallback is not used at all -- the beta-binomial scores
    // each sample against its OWN control, which is what makes a single-sample submission work.
    if (params.review_filter) {
        ch_no_file = Channel.fromPath("${projectDir}/assets/NO_FILE")

        // Rule 5, repeat context. Needs no controls, no cohort and no guide information, so it
        // is the one rule that works on day one for a guide never run before.
        ch_repeats = params.review_repeat_beds
            ? Channel.fromPath(params.review_repeat_beds.tokenize(','), checkIfExists: true).collect()
            : Channel.value([])

        // Rule 6, external DRAGEN noise panel. Like rule 5 it needs nothing from this run, but it
        // is a single optional file, so it uses the NO_FILE sentinel rather than an empty list.
        ch_snv_noise = params.review_snv_noise
            ? Channel.fromPath(params.review_snv_noise, checkIfExists: true)
            : ch_no_file

        // .collect() yields a value channel, so the same staged list feeds both review processes.
        ch_analysis_tsvs = GET_INDELS.out.indels_file.map{ meta, tsv -> tsv }.collect()

        REVIEW_FILTER(
            ch_analysis_tsvs,
            ch_repeats,
            ch_snv_noise
        )
        ch_versions = ch_versions.mix(REVIEW_FILTER.out.versions)

        // Review packet inputs, shared by the indel and breakend renderers. CRAM paths are
        // passed as a map rather than staged -- the queue names arbitrary samples, and staging
        // every cohort CRAM to draw a few dozen pictures would copy TBs.
        def want_bnd_snaps = params.review_filter_bnd && params.review_bnd_snapshots
        if (params.review_snapshots || want_bnd_snaps) {
            ch_review_cram_map = ch_dragen_files
                .map{ meta, dragenfiles ->
                    def crams = dragenfiles.findAll{ it ==~ /.*\.(cram)$/ }
                    def ed = crams.max{ it.toString().length() }
                    def ct = crams.min{ it.toString().length() }
                    "${meta.id}\t${ed}\t${ct}\n"
                }
                .collectFile(name: 'review_cram_map.tsv', sort: true)
                .first()

            // No .first() here: ch_fasta_reference is built with .collect(), so it is
            // already a value channel and .first() would only earn a warning.
            ch_review_fasta = ch_fasta_reference
                .map{ it.find{ f -> f ==~ /.*\.(fasta|fa)$/ } }
        }
        // .first() on the cram map: collectFile() yields a QUEUE channel, and it now feeds
        // two processes. A value channel is re-readable; a queue channel would let whichever
        // renderer ran first consume the only item and leave the other waiting forever.

        // The same triage applied to breakends rather than indels. Separate process because the
        // rules differ: a junction has two ends, so promiscuity replaces length diversity and the
        // noise panel is a BEDPE rather than a BED.
        if (params.review_filter_bnd) {
            ch_sv_noise = params.review_sv_noise
                ? Channel.fromPath(params.review_sv_noise, checkIfExists: true)
                : ch_no_file

            REVIEW_FILTER_BND(
                ch_analysis_tsvs,
                ch_sv_noise
            )
            ch_versions = ch_versions.mix(REVIEW_FILTER_BND.out.versions)

            // One figure per junction, not per queue row -- see the module header.
            if (params.review_bnd_snapshots) {
                REVIEW_BND_SNAPSHOTS(
                    REVIEW_FILTER_BND.out.queue,
                    ch_review_cram_map,
                    ch_review_fasta
                )
                ch_versions = ch_versions.mix(REVIEW_BND_SNAPSHOTS.out.versions)
            }
        }

        // Review packet: an IGV-style pileup per surviving site, edited over matched control.
        if (params.review_snapshots) {
            REVIEW_SNAPSHOTS(
                REVIEW_FILTER.out.queue,
                ch_review_cram_map,
                ch_review_fasta
            )
            ch_versions = ch_versions.mix(REVIEW_SNAPSHOTS.out.versions)
        }
    }

    BND_FROM_INDELS_TO_VCF (
        GET_INDELS.out.indels_file
    )
    ch_report_inputs = ch_report_inputs.mix(BND_FROM_INDELS_TO_VCF.out.vcf)
    ch_versions = ch_versions.mix(BND_FROM_INDELS_TO_VCF.out.versions)

    // Get transgene junctions
    GET_TRANSGENE_JUNCTIONS(ch_dragen_files,
                            ch_fasta_reference,
                            ch_transgene_name
    )
    ch_versions = ch_versions.mix(GET_TRANSGENE_JUNCTIONS.out.versions)

    TRANSGENE_TO_VCF(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSGENE_TO_VCF.out.versions)

    ANNOTATE_TRANSGENE_JUNCTIONS(
        TRANSGENE_TO_VCF.out.transgene_vcf,
        ch_fasta_reference,
        ch_vepcache,
        ch_cytobands
    )
    ch_report_inputs = ch_report_inputs.mix(ANNOTATE_TRANSGENE_JUNCTIONS.out.annotated_junctions)
    ch_versions = ch_versions.mix(ANNOTATE_TRANSGENE_JUNCTIONS.out.versions)

    // Make Circos plot
    TRANSFORM_TRANSGENE(GET_TRANSGENE_JUNCTIONS.out.transgene_file)
    ch_versions = ch_versions.mix(TRANSFORM_TRANSGENE.out.versions)

    MAKE_CIRCOS_PLOT(TRANSFORM_TRANSGENE.out.circos_input)
    // Add to report inputs
    ch_report_inputs = ch_report_inputs.mix(MAKE_CIRCOS_PLOT.out.plot)

    //
    // Generate CNA plots
    //
    GENERATE_CNA_BAF_PLOTS(ch_dragen_files)
    ch_report_inputs = ch_report_inputs.mix(GENERATE_CNA_BAF_PLOTS.out.plots) 
    ch_versions = ch_versions.mix(GENERATE_CNA_BAF_PLOTS.out.versions)

    // Make report JSON
    COMPILE_REPORT_JSON(
        ch_report_inputs
        .groupTuple()
        .map{ it -> [ it[0], it[1].flatten() ] }.view(),
        Channel.value("${new Date().getTime()}"))
    ch_versions = ch_versions.mix(COMPILE_REPORT_JSON.out.versions)


    // Render report
    MAKE_SCGE_REPORT(
        COMPILE_REPORT_JSON.out.json
        .join(GENERATE_CNA_BAF_PLOTS.out.plots)
        .join(MAKE_CIRCOS_PLOT.out.plot)
    )

    emit:
    versions = ch_versions  // channel: [ path(file) ]
    scge_report = MAKE_SCGE_REPORT.out.scge_report

}

