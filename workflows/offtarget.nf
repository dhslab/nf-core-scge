/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified CRISPR Off-Target Workflow  (entry: -entry OFFTARGET)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Self-contained. Leaves the default SCGE workflow untouched.

    samplesheet (sample,datatype{ecs|wgs},guide,edited_cram,control_cram,target_file,vcf)
       ├─ ecs rows ─► ECS_INDELS  @ hotspots ─────────────► ECS truth (indel_fraction = VAF)
       └─ wgs rows ─► WGS_WORKLIST ─► PON_OFFTARGET_FILTER ─► genome-wide homology-free worklist
                          │
       paired: HOTSPOT_TO_TABLE ─► SCORE_HOTSPOTS ─► BUILD_TRAINING_TABLE ─► training.tsv
                                                     └► RECALL_VS_VAF ─► recall_vs_vaf.{csv,png}
       all   : RECONCILE_OFFTARGET_REPORT ─► offtarget_report.csv  (is_hotspot / ecs_confirmed)

    Notes for the first real run:
      * WGS edited_cram must be the DRAGEN '<base>_tumor.cram'; the matched normal
        '<base>.cram' and '<base>.hard-filtered.vcf.gz' must sit beside it (the scripts
        derive them by naming convention). CRAMs/VCFs are read by absolute path from
        bind-mounted storage, so paths in the samplesheet must be absolute.
      * The paired training/recall arm is new code — validate its join on the first run
        (coordinate off-by-one between ECS and WGS is the usual culprit for an empty join).
*/

include { WGS_WORKLIST              } from '../modules/local/wgs_worklist.nf'
include { PON_OFFTARGET_FILTER      } from '../modules/local/pon_offtarget_filter.nf'
include { ECS_INDELS                } from '../modules/local/ecs_indels.nf'
include { HOTSPOT_TO_TABLE          } from '../modules/local/hotspot_to_table.nf'
include { SCORE_HOTSPOTS            } from '../modules/local/score_hotspots.nf'
include { BUILD_TRAINING_TABLE      } from '../modules/local/build_training_table.nf'
include { RECALL_VS_VAF             } from '../modules/local/recall_vs_vaf.nf'
include { RECONCILE_OFFTARGET_REPORT } from '../modules/local/reconcile_offtarget_report.nf'

workflow OFFTARGET_WORKFLOW {

    if (!params.input) { error "OFFTARGET: provide --input samplesheet.csv" }
    ch_ss = file(params.input, checkIfExists: true)

    // ---- preflight: detect run mode from the sheet and validate datatypes up front, so a
    // mistyped or single-arm sheet fails with a clear message (or runs the right arms) instead
    // of silently sending empty inputs downstream. Read synchronously — it's a small local file.
    def rows = ch_ss.readLines().findAll { it.trim() }
    def hdr  = rows[0].split(',', -1).collect { it.trim().toLowerCase() }
    def dti  = hdr.indexOf('datatype')
    if (dti < 0) { error "OFFTARGET: --input has no 'datatype' column" }
    def dts  = rows.drop(1).collect { (it.split(',', -1)[dti] ?: '').trim().toLowerCase() }
    def bad  = dts.findAll { it && it != 'ecs' && it != 'wgs' }.unique()
    if (bad) { error "OFFTARGET: unknown datatype(s) ${bad} in --input (expected ecs|wgs)" }
    def n_ecs = dts.count { it == 'ecs' }
    def n_wgs = dts.count { it == 'wgs' }
    if (n_ecs + n_wgs == 0) { error "OFFTARGET: no ecs/wgs rows in --input" }
    def mode = (n_ecs && n_wgs) ? 'paired' : (n_wgs ? 'wgs_only' : 'ecs_only')
    log.info "OFFTARGET mode: ${mode}  (${n_ecs} ecs, ${n_wgs} wgs rows)"

    ch_rows = Channel.fromPath(ch_ss)
        | splitCsv(header: true)
        | branch { row ->
            ecs: (row.datatype ?: '').toLowerCase() == 'ecs'
            wgs: (row.datatype ?: '').toLowerCase() == 'wgs'
        }

    // ---- ECS arm: truth at hotspots ----
    ch_ecs_truth_files = Channel.empty()
    if (n_ecs > 0) {
        ch_ecs_in = ch_rows.ecs.map { row ->
            tuple([id: row.sample, guide: row.guide],
                  row.edited_cram, row.control_cram, row.target_file)
        }
        ECS_INDELS(ch_ecs_in, params.fasta)
        ch_ecs_truth_files = ECS_INDELS.out.indels_file.map { meta, tsv -> tsv }
    }

    // ---- WGS arm: genome-wide homology-free discovery + PoN ----
    if (n_wgs > 0) {
        // build the sample<TAB>tumor_cram map the scripts consume (single source of truth)
        ch_cram_map = ch_rows.wgs
            .map { row -> "${row.sample}\t${row.edited_cram}" }
            .collectFile(name: 'wgs_cram_map.tsv', newLine: true, sort: true)

        WGS_WORKLIST(ch_cram_map, file(params.offtarget_shape_model), params.fasta)
        PON_OFFTARGET_FILTER(WGS_WORKLIST.out.worklist, ch_cram_map, params.fasta)
    }

    // ---- paired arm: ECS-truth ⋈ WGS-features → training table + recall-vs-VAF ----
    // Only meaningful when BOTH arms are present.
    ch_truth = Channel.empty()
    if (n_ecs > 0 && n_wgs > 0) {
        HOTSPOT_TO_TABLE(ch_ecs_truth_files.collect(), ch_ss)
        SCORE_HOTSPOTS(HOTSPOT_TO_TABLE.out.table, ch_cram_map,
                       file(params.offtarget_shape_model), params.fasta)
        BUILD_TRAINING_TABLE(SCORE_HOTSPOTS.out.scores, HOTSPOT_TO_TABLE.out.truth, ch_ss)
        RECALL_VS_VAF(BUILD_TRAINING_TABLE.out.training)
        ch_truth = HOTSPOT_TO_TABLE.out.truth
    }

    // ---- reconciled report: whenever there's a WGS worklist to annotate (truth optional) ----
    if (n_wgs > 0) {
        ch_truth_opt = ch_truth.ifEmpty(file("${projectDir}/assets/NO_FILE"))
        RECONCILE_OFFTARGET_REPORT(PON_OFFTARGET_FILTER.out.worklist, ch_truth_opt)
    }
}
