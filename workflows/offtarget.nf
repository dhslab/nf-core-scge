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
include { GENERATE_HOTSPOTS          } from '../subworkflows/local/generate_hotspots.nf'

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

    // ---- preflight: every ecs row needs a target_file OR a spacer (to auto-generate one via
    // GENERATE_HOTSPOTS). Count the auto-generate rows synchronously here so the hotspot
    // subworkflow — whose PREP_CASOFFINDER_REF runs once regardless of guide count — is only
    // invoked when a row actually needs it (never on a fully target_file-provided sheet).
    def n_ecs_auto = 0
    if (n_ecs > 0) {
        def sfi = hdr.indexOf('spacer')
        def tfi = hdr.indexOf('target_file')
        def si2 = hdr.indexOf('sample')
        def eprob = []
        rows.drop(1).each { line ->
            def cols = line.split(',', -1)
            if ((cols[dti] ?: '').trim().toLowerCase() != 'ecs') return
            def sample = si2 >= 0 ? (cols[si2] ?: '').trim() : '?'
            def tf = tfi >= 0 ? (cols[tfi] ?: '').trim() : ''
            def sp = sfi >= 0 ? (cols[sfi] ?: '').trim() : ''
            if (!tf && !sp) {
                eprob << "row '${sample}': ecs row needs a 'target_file' or a 'spacer' (to auto-generate one)"
            } else if (!tf && sp) {
                n_ecs_auto++
            }
        }
        if (eprob) { error "OFFTARGET: ecs preflight failed:\n  " + eprob.join('\n  ') }
        if (n_ecs_auto > 0) { log.info "OFFTARGET: auto-generating hotspots for ${n_ecs_auto} ecs row(s) from their spacer" }
    }

    // ---- preflight: WGS rows depend on DRAGEN sidecar files derived from the tumor CRAM name
    // by convention (bin/worklist_from_vcf.py, bin/score.py): the matched normal '<base>.cram'
    // and the somatic VCF '<base>.hard-filtered.vcf.gz' must sit beside '<base>_tumor.cram'.
    // Validate them up front with a named-file error, so a misnamed or missing sidecar fails
    // clearly here instead of surfacing as a confusing empty worklist or a mid-run task crash.
    // Skipped under -stub-run and -preview, where inputs are placeholder paths that need not exist.
    if (n_wgs > 0 && !workflow.stubRun && !workflow.preview) {
        def eci = hdr.indexOf('edited_cram')
        def si  = hdr.indexOf('sample')
        if (eci < 0) { error "OFFTARGET: --input has no 'edited_cram' column" }
        def problems = []
        rows.drop(1).each { line ->
            def cols = line.split(',', -1)
            if ((cols[dti] ?: '').trim().toLowerCase() != 'wgs') return
            def sample = si >= 0 ? (cols[si] ?: '').trim() : '?'
            def tumor  = (cols[eci] ?: '').trim()
            if (!tumor.endsWith('_tumor.cram')) {
                problems << "row '${sample}': WGS edited_cram must be a DRAGEN '<base>_tumor.cram' (got '${tumor}')"
                return
            }
            def base = tumor - ~/_tumor\.cram$/
            [ (tumor)                                    : 'tumor CRAM',
              ("${base}.cram".toString())                : 'matched-normal CRAM',
              ("${base}.hard-filtered.vcf.gz".toString()): 'DRAGEN somatic VCF' ].each { p, what ->
                if (!file(p).exists()) problems << "row '${sample}': missing ${what}: ${p}"
            }
        }
        if (problems) {
            error "OFFTARGET: DRAGEN sidecar preflight failed. The WGS arm derives the matched\n" +
                  "normal and somatic VCF from the tumor CRAM name; each must exist beside it:\n  " +
                  problems.join('\n  ')
        }
    }

    ch_rows = Channel.fromPath(ch_ss)
        | splitCsv(header: true)
        | branch { row ->
            ecs: (row.datatype ?: '').toLowerCase() == 'ecs'
            wgs: (row.datatype ?: '').toLowerCase() == 'wgs'
        }

    // ---- ECS arm: truth at hotspots ----
    ch_ecs_truth_files = Channel.empty()
    if (n_ecs > 0) {
        // Rows that already point at a target_file VCF use it directly.
        ch_ecs_provided = ch_rows.ecs
            .filter { row -> row.target_file?.trim() }
            .map { row -> tuple([id: row.sample, guide: row.guide],
                                row.edited_cram, row.control_cram, row.target_file) }

        if (n_ecs_auto > 0) {
            // Rows with a spacer but no target_file: auto-generate the hotspot VCF once per
            // guide, then broadcast it to every ecs replicate of that guide.
            ch_ecs_auto = ch_rows.ecs.filter { row -> !row.target_file?.trim() && row.spacer?.trim() }

            ch_guides = ch_ecs_auto
                .map { row -> tuple(row.guide, row.spacer.trim().toUpperCase(),
                                    (row.pam?.trim() ?: params.offtarget_pam)) }
                .unique()
                .map { guide, spacer, pam -> tuple([id: guide], spacer, pam,
                                                    file("${projectDir}/assets/NO_IDT")) }

            GENERATE_HOTSPOTS(ch_guides)

            ch_ecs_auto_in = ch_ecs_auto
                .map { row -> tuple(row.guide, [id: row.sample, guide: row.guide],
                                    row.edited_cram, row.control_cram) }
                .combine(GENERATE_HOTSPOTS.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }, by: 0)
                .map { guide, meta, ed, ctl, vcf -> tuple(meta, ed, ctl, vcf) }

            ch_ecs_in = ch_ecs_provided.mix(ch_ecs_auto_in)
        } else {
            ch_ecs_in = ch_ecs_provided
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
