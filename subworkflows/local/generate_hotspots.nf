/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE_HOTSPOTS  —  gRNA -> predicted off-target sites (target_file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per guide: enumerate off-target sites with Cas-OFFinder (and, when enabled, CRISPRme),
    combine them into <guide>.targets.csv, then expand to <guide>.targets.vcf — the VCF the
    OFFTARGET arm consumes as a per-row target_file.

        (meta, spacer, pam, idt)
           ├─ CASOFFINDER (2bit) ──► <guide>.casoffinder.txt
           └─ CRISPRME  (opt)    ──► <guide>.crisprme.tsv
                 └─ COMBINE_OFFTARGET_SITES ──► <guide>.targets.csv
                       └─ TARGETS_CSV_TO_VCF ──► <guide>.targets.vcf
*/

include { PREP_CASOFFINDER_REF    } from '../../modules/local/prep_casoffinder_ref.nf'
include { CASOFFINDER             } from '../../modules/local/casoffinder.nf'
include { CRISPRME                } from '../../modules/local/crisprme.nf'
include { COMBINE_OFFTARGET_SITES } from '../../modules/local/combine_offtarget_sites.nf'
include { TARGETS_CSV_TO_VCF      } from '../../modules/local/targets_csv_to_vcf.nf'

workflow GENERATE_HOTSPOTS {
    take:
    ch_guides            // tuple(meta[id], spacer, pam, idt_path_or_NO_FILE)

    main:
    ch_versions = Channel.empty()
    // Distinct sentinels for the two optional COMBINE slots — same filename in both slots
    // would collide when staged into one task.
    def no_crisprme = file("${projectDir}/assets/NO_FILE")

    // Reference the search runs against: a prebuilt .2bit, or built once with faToTwoBit.
    if (params.casoffinder_2bit) {
        ch_2bit = Channel.value(file(params.casoffinder_2bit, checkIfExists: true))
    } else {
        PREP_CASOFFINDER_REF(file(params.fasta, checkIfExists: true))
        ch_2bit     = PREP_CASOFFINDER_REF.out.twobit.first()
        ch_versions = ch_versions.mix(PREP_CASOFFINDER_REF.out.versions)
    }

    ch_guide_only = ch_guides.map { meta, spacer, pam, idt -> tuple(meta, spacer, pam) }

    // ---- Cas-OFFinder (always) ----
    CASOFFINDER(ch_guide_only, ch_2bit)
    ch_versions = ch_versions.mix(CASOFFINDER.out.versions)

    // ---- CRISPRme (optional; Cas-OFFinder-only by default) ----
    if (params.run_crisprme) {
        ch_index = Channel.value(file(params.crisprme_index_dir, checkIfExists: true))
        CRISPRME(ch_guide_only, ch_index)
        ch_versions = ch_versions.mix(CRISPRME.out.versions)
        ch_cme = CRISPRME.out.hits
    } else {
        ch_cme = ch_guide_only.map { meta, s, p -> tuple(meta, no_crisprme) }
    }

    // ---- combine per guide, then expand to the target VCF ----
    ch_idt = ch_guides.map { meta, spacer, pam, idt -> tuple(meta, idt) }
    ch_combine = CASOFFINDER.out.hits.join(ch_cme).join(ch_idt)

    COMBINE_OFFTARGET_SITES(ch_combine)
    ch_versions = ch_versions.mix(COMBINE_OFFTARGET_SITES.out.versions)

    TARGETS_CSV_TO_VCF(COMBINE_OFFTARGET_SITES.out.sites,
                       file(params.fasta), file("${params.fasta}.fai"))
    ch_versions = ch_versions.mix(TARGETS_CSV_TO_VCF.out.versions)

    emit:
    sites    = COMBINE_OFFTARGET_SITES.out.sites   // (meta, <guide>.targets.csv)
    vcf      = TARGETS_CSV_TO_VCF.out.vcf           // (meta, <guide>.targets.vcf)
    versions = ch_versions
}
