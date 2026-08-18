#!/usr/bin/env python3
"""panel_overlap.py — measure the DRAGEN systematic-noise panels as standalone filters.

Why this exists
---------------
Two questions were asked of the panels and answered once, in prose, in
docs/CALLER_INTEGRATION_LOG.md — but the script that produced those numbers was never kept, so
nothing could be re-measured when the caller changed underneath them. This restores that
capability and answers both questions on whatever run is handed to it:

  * indels — how many calls overlap a noisy SNV locus in
    `IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz`, and does that overlap track truth?
  * breakends — how many junctions overlap a record in each of the three
    `*_systematic_noise.sv.bedpe.gz` panels?

What it is NOT
--------------
This does not filter anything and is not called by the pipeline. It reuses the shipped matching
code (`review_filter.load_snv_noise` / `snv_noise_mask`, `review_filter_bnd.load_sv_noise` /
`in_noise`) so that what it reports is what the pipeline would do, not a second implementation
that could drift.

The one deliberate departure, and why
-------------------------------------
The shipped SNV loader keeps ONE record per position — the one with the most donors — because the
shipped rule only ever asks a single yes/no question. That is lossy for the question asked here:
an SNV-only record with 20 donors shadows an indel record with 2 at the same coordinate, so the
counts for "any indel-capable record" would come out too low. This script therefore loads every
record per position via `load_snv_noise_all`, reports the three overlap definitions honestly, and
then re-runs the shipped loader + mask and reports whether the two agree. A disagreement is a real
defect in the shipped rule, not a bookkeeping detail, so it is printed rather than smoothed over.

BREAKENDS: the on-target exemption is DISABLED here, on purpose. The shipped rule exempts
`is_target == 1` junctions from the panel check, and on this cohort every gated junction is
on-target — so with the exemption on, the answer to "how many junctions does the panel flag?" is
trivially zero and tells you nothing about the panel. The question is about the panel's
discrimination, so the exemption is turned off and that is stated in the output.

usage:
  panel_overlap.py --tables 'results_cart_bnd/*/*.offtarget_analysis.tsv' \
      --queue-all results_cart_bnd/review/review_queue_all.tsv \
      --snv-noise .../IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz \
      --truth-wgs '.../cart_wgs_merged.xlsx' \
      --sv-noise .../WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz \
      --sv-noise .../IDPF_WGS_hg38_v3.0.0_systematic_noise.sv.bedpe.gz \
      --sv-noise .../WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz
"""
import argparse
import collections
import glob
import bisect
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import noise_model as nm                                              # noqa: E402
import review_filter as rf                                            # noqa: E402
import review_filter_bnd as rfb                                       # noqa: E402

SLOP = rf.SNV_NOISE_SLOP
MIN_DONORS = 3


# ---------------------------------------------------------------------------------------------
# SNV panel
# ---------------------------------------------------------------------------------------------
def load_snv_noise_all(path, positions, slop=SLOP):
    """Same file and same columns as review_filter.load_snv_noise, but keeps EVERY record.

    -> {(chrom, pos): [(max_noise, n_donors, alleles), ...]}

    The shipped loader keeps only the record with the most donors; see the module docstring for
    why that is not good enough to count overlap definitions against each other.
    """
    want = collections.defaultdict(set)
    for c, e in positions:
        for d in range(-slop, slop + 1):
            want[c].add(int(e) + d)
    hits = collections.defaultdict(list)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            s = want.get(p[0])
            if s is None:
                continue
            try:
                pos = int(p[2])
                if pos not in s:
                    continue
                hits[(p[0], pos)].append((float(p[4]), int(p[6]), p[5]))
            except ValueError:
                continue
    return hits


def is_indel_record(rec):
    """The panel's allele column is a comma list; D and I are its indel codes."""
    return "D" in rec[2] or "I" in rec[2]


def overlap_flags(hits_all, chroms, ends, min_donors=MIN_DONORS, slop=SLOP):
    """Three nested overlap definitions per site, as boolean arrays.

      any     -- the panel says anything at all here
      indel   -- the panel says something indel-capable here
      shipped -- indel-capable AND recurrent in >= min_donors donors (what review_filter uses)
    """
    a, i, s = [], [], []
    for c, e in zip(chroms, ends):
        recs = []
        for d in range(-slop, slop + 1):
            recs.extend(hits_all.get((c, int(e) + d), ()))
        ind = [r for r in recs if is_indel_record(r)]
        a.append(bool(recs))
        i.append(bool(ind))
        s.append(any(r[1] >= min_donors for r in ind))
    return np.array(a), np.array(i), np.array(s)


def confusion(flag, label):
    """2x2 for the panel used as a standalone filter.

    The panel flags NOISE, so a flagged row is a predicted negative. Reported from the point of
    view of the filter's job: does flagging remove rejected rows and spare confirmed ones?
    """
    flag, label = np.asarray(flag, bool), np.asarray(label, int)
    tp = int((flag & (label == 0)).sum())        # correctly flagged an artifact
    fp = int((flag & (label == 1)).sum())        # flagged a real edit -- the expensive error
    fn = int((~flag & (label == 0)).sum())       # missed an artifact
    tn = int((~flag & (label == 1)).sum())       # correctly spared a real edit
    div = lambda a, b: a / b if b else float("nan")   # noqa: E731
    return dict(flagged_rejected=tp, flagged_confirmed=fp, missed_rejected=fn,
                spared_confirmed=tn,
                sensitivity=div(tp, tp + fn),     # of artifacts, how many caught
                specificity=div(tn, tn + fp),     # of real edits, how many spared
                precision=div(tp, tp + fp))       # of flagged rows, how many really artifacts


def snv_report(df, queue_all, snv_noise, truth_wgs, min_donors, out, cache=None):
    print("\n" + "=" * 94)
    print("ASK #2 — the SNV panel as an indel filter")
    print("=" * 94)
    print(f"panel : {snv_noise}")
    print(f"slop  : +/-{SLOP} bp     shipped rule: indel-capable AND >= {min_donors} donors")

    positions = set(zip(df.chrom.astype(str), df.end.astype(int)))
    if cache and os.path.exists(cache):
        hits_all = pd.read_pickle(cache)
        print(f"\npanel hits loaded from cache {cache}")
    else:
        print(f"\nstreaming the panel for {len(positions)} distinct positions (~1 GB, 2-3 min)...")
        hits_all = load_snv_noise_all(snv_noise, positions)
        if cache:
            pd.to_pickle(dict(hits_all), cache)
    print(f"panel records retained at queried positions: {sum(len(v) for v in hits_all.values())} "
          f"at {len(hits_all)} positions")

    # The shipped path, run exactly as review_filter runs it, for the agreement check below.
    # Derived from the same full record set so the file is streamed once, not twice: keeping only
    # the highest-donor record per position is precisely what load_snv_noise does.
    hits_shipped = {k: max(v, key=lambda r: r[1]) for k, v in hits_all.items()}

    levels = [("all site-rows", df)]
    if queue_all is not None:
        levels.append(("gated rows", queue_all))
        kept = queue_all[queue_all.why_dropped.fillna("") == ""]
        levels.append(("review queue", kept))

    rows = []
    for name, d in levels:
        any_, ind, ship = overlap_flags(hits_all, d.chrom.astype(str), d.end.astype(int),
                                        min_donors)
        rows.append(dict(population=name, n=len(d),
                         any_record=int(any_.sum()), indel_record=int(ind.sum()),
                         shipped_rule=int(ship.sum())))
    tab = pd.DataFrame(rows)
    for c in ("any_record", "indel_record", "shipped_rule"):
        tab[c + "_pct"] = (tab[c] / tab.n * 100).round(2)
    print("\n-- overlap counts at three population levels --")
    print(tab.to_string(index=False))

    # Agreement between the two loaders. A mismatch means the shipped rule missed an indel record
    # because a higher-donor SNV record at the same position shadowed it.
    _, _, ship_full = overlap_flags(hits_all, df.chrom.astype(str), df.end.astype(int), min_donors)
    ship_pipe = np.array(rf.snv_noise_mask(hits_shipped, df.chrom.astype(str),
                                           df.end.astype(int), min_donors))
    shadowed = int((ship_full & ~ship_pipe).sum())
    print(f"\n-- shipped-loader agreement --")
    print(f"shipped rule, all records kept : {int(ship_full.sum())}")
    print(f"shipped rule, as the pipeline runs it : {int(ship_pipe.sum())}")
    if shadowed:
        print(f"*** {shadowed} site(s) carry an indel record that the one-record-per-position "
              f"loader discards ***\n    (an SNV-only record with more donors shadows it; the "
              f"pipeline under-flags by this much)")
    else:
        print("identical — no indel record is shadowed at any queried position on this run")

    # ---- truth ----
    if truth_wgs:
        lab = nm.wgs_curated_label(truth_wgs)
        d = df.copy()
        d["guide"] = d.sample_name.map(nm.guide_of)
        d["chrom"] = d.chrom.astype(str)
        m = d.merge(lab, on=["guide", "chrom", "start"], how="inner")
        n_neg = int((m.curated_label == 0).sum())
        print(f"\n-- the panel against the curated WGS label --")
        print(f"joined {len(m)} labelled rows: {int((m.curated_label==1).sum())} confirmed / "
              f"{n_neg} human-rejected")
        # Precision on a 75%-negative set is flattered by prevalence, so print what a filter that
        # flagged EVERY row would score. Anything at or below that line has learned nothing.
        print(f"base rate: flagging every row would score precision {n_neg/len(m):.3f} "
              f"at sensitivity 1.000")
        any_, ind, ship = overlap_flags(hits_all, m.chrom, m.end.astype(int), min_donors)
        for defn, flag in (("any panel record", any_), ("indel-capable record", ind),
                           (f"shipped rule (>={min_donors} donors & D/I)", ship)):
            c = confusion(flag, m.curated_label.values)
            print(f"  {defn:38s} flagged {int(flag.sum()):4d}   "
                  f"sens {c['sensitivity']:.3f}  spec {c['specificity']:.3f}  "
                  f"prec {c['precision']:.3f}   (real edits flagged: {c['flagged_confirmed']})")

    # ---- the decomposition that settles the substitution question ----
    if queue_all is not None:
        print("\n-- can the panel stand in for rule 1 (germline in the matched control)? --")
        wd = queue_all.why_dropped.fillna("")
        germ = queue_all[wd == "germline (present in control)"]
        any_, ind, ship = overlap_flags(hits_all, germ.chrom.astype(str), germ.end.astype(int),
                                        min_donors)
        print(f"rows rule 1 dropped as germline: {len(germ)}")
        print(f"  of those, panel flags: any {int(any_.sum())}   indel-capable {int(ind.sum())}"
              f"   shipped rule {int(ship.sum())}")
        print("  Germline is a property of ONE donor. A 46-donor panel reports a population")
        print("  average, so it cannot represent one donor's genotype at any panel size.")
        # Which germline sites the panel DOES see is the informative half: a population panel can
        # only recognise a population-common allele, so the overlap should concentrate on the
        # sites where many donors carry it.
        dons = []
        for c, e in zip(germ.chrom.astype(str), germ.end.astype(int)):
            best = 0
            for d in range(-SLOP, SLOP + 1):
                for r in hits_all.get((c, int(e) + d), ()):
                    best = max(best, r[1])
            dons.append(best)
        dons = np.array(dons)
        print(f"  donor support at the germline sites the panel does see: "
              f"median {int(np.median(dons[dons>0])) if (dons>0).any() else 0}, "
              f"max {int(dons.max())}  (of 46 panel donors)")
        print(f"  germline rows invisible to the panel entirely: {int((dons==0).sum())}"
              f"/{len(germ)} = {(dons==0).mean()*100:.0f}%")

    if out:
        tab.to_csv(out, sep="\t", index=False)
        print(f"\nwrote {out}")


# ---------------------------------------------------------------------------------------------
# SV panels
# ---------------------------------------------------------------------------------------------
def sv_donor_index(path):
    """Breakpoint intervals with the set of DONORS supporting each, plus total bp covered.

    The BEDPE carries no donor-count column and no `##PON SAMPLES` header, so there is no direct
    analogue of the SNV panel's `n_donors` -- which is the very statistic that makes the SNV panel
    usable, since it is what separates a recurrent artifact from a one-donor coincidence.

    It is recoverable anyway: field 7 is a candidate NAME that embeds the donor it came from,
    e.g. `ImpreciseNoiseCandidate_LP7108672-DNA_A06_42_DRAGEN:BND:...`. One record is one candidate
    from one donor, so pooling records by interval and counting distinct donors reconstructs the
    missing column. This measures whether doing so would rescue the large panels.

    -> ({chrom: (starts, ends, donorsets)}, covered_bp, n_records)
    """
    iv = collections.defaultdict(list)
    n_rec = 0
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            name = p[6]
            # strip the NoiseCandidate prefix and the trailing _<idx>_DRAGEN:<type>... suffix
            d = name.split("_", 1)[-1]
            d = d.split("_DRAGEN:")[0]
            d = d.rsplit("_", 1)[0] if d.rsplit("_", 1)[-1].isdigit() else d
            try:
                iv[p[0]].append((int(p[1]), int(p[2]), d))
                iv[p[3]].append((int(p[4]), int(p[5]), d))
            except ValueError:
                continue
            n_rec += 1
    merged, covered = {}, 0
    for c, v in iv.items():
        v.sort()
        st, en, ds = [], [], []
        for s, e, d in v:
            if st and s <= en[-1]:
                en[-1] = max(en[-1], e)
                ds[-1].add(d)
            else:
                st.append(s); en.append(e); ds.append({d})
        merged[c] = (st, en, ds)
        covered += sum(b - a for a, b in zip(st, en))
    return merged, covered, n_rec


def donors_at(merged, chrom, pos, slop):
    if chrom not in merged:
        return set()
    st, en, ds = merged[chrom]
    i = bisect.bisect_right(st, pos + slop) - 1
    if i >= 0 and st[i] - slop <= pos <= en[i] + slop:
        return ds[i]
    return set()


def sv_report(junctions, panels, slops, extra_artifacts=None):
    print("\n" + "=" * 94)
    print("ASK #3 — the SV panels as a breakend filter")
    print("=" * 94)
    print("ON-TARGET EXEMPTION DISABLED. Every gated junction on this cohort is is_target==1, so")
    print("the shipped rule flags nothing by construction; that would answer the wrong question.")

    gated = junctions[junctions.reads >= rfb.MIN_READS]
    sets = [("all junctions", junctions), (f"gated (reads>={rfb.MIN_READS})", gated)]
    if extra_artifacts is not None and len(extra_artifacts):
        sets.append(("dropped junctions (artifacts)", extra_artifacts))

    rows, donor_rows = [], []
    gate_slop = rfb.SV_NOISE_SLOP
    for path in panels:
        name = os.path.basename(path).replace("_systematic_noise.sv.bedpe.gz", "")
        print(f"\nloading {name} ...")
        merged = rfb.load_sv_noise(path)
        n_iv = sum(len(v[0]) for v in merged.values())
        dmerged, covered, n_rec = sv_donor_index(path)
        print(f"  {n_rec} BEDPE records -> {n_iv} merged breakpoint intervals over "
              f"{len(merged)} contigs")
        print(f"  genome covered by those intervals: {covered/1e6:.1f} Mb "
              f"({covered/3.1e9*100:.2f}% of hg38)")
        for slop in slops:
            for label, s in sets:
                flag = [rfb.in_noise(merged, c, p, slop) or rfb.in_noise(merged, c2, p2, slop)
                        for c, p, c2, p2 in zip(s.chrom, s.pos, s.chrom2, s.pos2)]
                rows.append(dict(panel=name, slop=slop, population=label, n=len(s),
                                 flagged=int(np.sum(flag)), covered_mb=round(covered / 1e6, 1)))

        # Donor support behind the hits on the gated junctions, at the shipped slop.
        nd = []
        for c, p, c2, p2 in zip(gated.chrom, gated.pos, gated.chrom2, gated.pos2):
            d = donors_at(dmerged, c, p, gate_slop) | donors_at(dmerged, c2, p2, gate_slop)
            nd.append(len(d))
        nd = np.array(nd)
        donor_rows.append(dict(panel=name, hit=int((nd > 0).sum()),
                               single_donor=int((nd == 1).sum()),
                               ge2=int((nd >= 2).sum()), ge3=int((nd >= 3).sum()),
                               max_donors=int(nd.max()) if len(nd) else 0))

    tab = pd.DataFrame(rows)
    tab["pct"] = (tab.flagged / tab.n * 100).round(1)
    print("\n-- flagged counts (either end in a panel interval) --")
    for slop in slops:
        print(f"\nslop {slop} bp:")
        print(tab[tab.slop == slop].pivot(index="panel", columns="population",
                                          values="flagged").to_string())

    print(f"\n-- donor support behind the hits on the {len(gated)} gated junctions "
          f"(slop {gate_slop}) --")
    print("reconstructed from the BEDPE name field; see sv_donor_index()")
    print(pd.DataFrame(donor_rows).to_string(index=False))
    print("\nIf a panel's hits are mostly SINGLE-donor, a min_donors rule -- the same thing that")
    print("makes the SNV panel usable -- would restore its specificity. If they are multi-donor,")
    print("the panel genuinely covers these loci and no threshold rescues it.")
    return tab


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True,
                    help="glob for *.offtarget_analysis.tsv (quote it)")
    ap.add_argument("--queue-all", help="review_queue_all.tsv (gated rows with why_dropped)")
    ap.add_argument("--snv-noise", help="systematic_noise.snv.bed.gz")
    ap.add_argument("--min-donors", type=int, default=MIN_DONORS)
    ap.add_argument("--truth-wgs", help="cart_wgs_merged.xlsx")
    ap.add_argument("--sv-noise", action="append", default=[],
                    help="systematic_noise.sv.bedpe.gz (repeatable)")
    ap.add_argument("--sv-slop", type=int, nargs="+", default=[0, 50, 200])
    ap.add_argument("--extra-bnd-queue", action="append", default=[],
                    help="bnd_review_queue_all.tsv from another run, to recover an artifact set "
                         "(rows with a non-empty why_dropped)")
    ap.add_argument("--snv-cache", help="pickle the streamed panel hits here and reuse them")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    paths = sorted(glob.glob(a.tables))
    if not paths:
        sys.exit(f"no tables matched {a.tables}")
    df = nm.load_tables(paths)
    print(f"tables {len(paths)}   site-rows {len(df)}")

    queue_all = pd.read_csv(a.queue_all, sep="\t") if a.queue_all else None
    if queue_all is not None:
        print(f"gated rows {len(queue_all)}   queue "
              f"{int((queue_all.why_dropped.fillna('') == '').sum())}")

    os.makedirs(a.outdir, exist_ok=True)

    if a.snv_noise:
        snv_report(df, queue_all, a.snv_noise, a.truth_wgs, a.min_donors,
                   os.path.join(a.outdir, "panel_overlap_snv.tsv"), a.snv_cache)

    if a.sv_noise:
        events = []
        for p in paths:
            d = pd.read_csv(p, sep="\t")
            events += rfb.parse_bnds(d, os.path.basename(p).split(".offtarget")[0])
        junctions = pd.DataFrame(events)
        print(f"\njunction entries {len(junctions)}")

        arts = []
        for p in a.extra_bnd_queue:
            q = pd.read_csv(p, sep="\t")
            arts.append(q[q.why_dropped.fillna("") != ""])
        extra = pd.concat(arts, ignore_index=True) if arts else None

        tab = sv_report(junctions, a.sv_noise, a.sv_slop, extra)
        o = os.path.join(a.outdir, "panel_overlap_sv.tsv")
        tab.to_csv(o, sep="\t", index=False)
        print(f"\nwrote {o}")


if __name__ == "__main__":
    main()
