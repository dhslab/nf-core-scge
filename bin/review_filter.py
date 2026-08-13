#!/usr/bin/env python3
"""Automate the manual review of CRISPR off-target hotspot calls.

Takes find_edited_reads.py output (*.offtarget_analysis.tsv) and returns the subset a human
still needs to look at. Four rules, each one a statement you can defend:

  1. the matched control is clean          -- if the variant is in the unedited sample from the
                                              same donor, Cas9 did not make it
  2. within 10 bp of a PAM position        -- Cas9 cuts ~3 bp from the PAM; real edits sit at
                                              the cut, background indels do not
  3. at least 3 distinct indel lengths     -- NHEJ produces a spectrum of deletion sizes; one
                                              length repeated across every read is an artifact
  4. not a known-bad site                  -- a site carrying indels in samples that were never
                                              edited is germline or a repeat, not an edit
  5. not in a repeat region                -- the off-target panel is built from homology, so it
                                              is enriched for the sequence aligners misplace
                                              indels in

Measured on the 25-sample CAR-T WGS cohort against the corrected truth set (61 real edits,
177 artifacts):

    review queue 238 -> 62 rows, 61/61 real edits retained
    precision 0.256 -> 0.984, recall 1.000

RULE 1 REQUIRES THE FIXED CALLER. Before the control-reporting fix, find_edited_reads.py summed
control support *after* the -x/--max-in-control filter had removed the events carrying it, so
control_indel_reads was 0 at every site and this rule was a no-op. The script warns if it sees
that signature.

RULE 4 HAS TWO SOURCES, and which one you get matters:
  --pon    a panel of normals built from unedited samples (bin/build_offtarget_pon.py).
           Uses no guide information, so it works for a SINGLE GUIDE. Preferred.
  (auto)   cross-guide recurrence -- a site appearing under >=2 distinct guides is a bad
           region. Requires several differently-guided samples in ONE invocation. Used only
           as a fallback when no PoN is given.
On the CAR-T cohort the PoN is strictly better (63 rows / 2 FP vs 65 / 4) and subsumes the
cross-guide rule entirely.

WITH NEITHER A PoN NOR CONTROLS, use --repeats. Repeat annotation needs no cohort, no guide
context and no unedited sample, and on this cohort it stands in for rule 4 at 64 rows / 3 FP
(precision 0.953) -- worse than a PoN but far better than nothing, and it works on day one for
a guide never run before. The GATK 1000g PoN was tested for this role and rejected: it flags 15
artifacts but also 4 of the 61 real edits, because it is a Mutect2 SNV panel from blood normals
rather than an indel-artifact map.

usage:
  review_filter.py IN.tsv [IN2.tsv ...] -o queue.tsv [--pon assets/offtarget_pon.tsv]
                   [--repeats rmsk.bed trf.bed]
"""
import argparse
import bisect
import collections
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

GATE_READS, GATE_VAF = 10, 0.05
SNV_NOISE_MIN_DONORS = 3
SNV_NOISE_SLOP = 2
MAX_CUT_DIST, MIN_DISTINCT_LEN, MAX_GUIDES = 10, 3, 2
MAX_CONTROL_VAF = 0.05
# Default AQ cut for --noise-model. Measured on the 32-sample cohort with the depth floor on:
# 3, 5 and 8 all reproduce the panel-of-normals result exactly (91 rows, 64 confirmed, 9 rejected,
# precision 0.877); 10 drops a confirmed edit. 5 sits in the middle of that plateau.
AQ_MIN = 5.0


def _indel_len(ref, alt):
    m = re.fullmatch(r'DEL(\d+)', str(alt))
    if m:
        return -int(m.group(1))
    m = re.fullmatch(r'INS(\d+)', str(alt))
    if m:
        return int(m.group(1))
    if str(alt).startswith('<') or str(ref) == '.' or str(alt) == '.':
        return None
    return len(str(alt)) - len(str(ref))


def parse_events(s):
    """indel_info -> [(distance_to_pam, indel_length, supporting_reads), ...]

    Field layout per ';'-separated event:
      chrom|pos|chrom2|pos2|strands|ref|alt|distance|distance2|counts|control_alt_counts
    'distance' (index 7) is min(|pos - PAM_position|) over the site's PAM positions.
    """
    out = []
    if not isinstance(s, str):
        return out
    for ev in s.split(';'):
        f = ev.split('|')
        if len(f) < 10:
            continue
        try:
            dist, cnt = abs(int(f[7])), int(f[9])
        except ValueError:
            continue
        L = _indel_len(f[5], f[6])
        if L:
            out.append((dist, L, cnt))
    return out


def add_features(df):
    cut, ndl = [], []
    for s in df.indel_info:
        ev = parse_events(s)
        if not ev:
            cut.append(np.nan)
            ndl.append(0)
            continue
        cut.append(min(e[0] for e in ev))
        ndl.append(len({e[1] for e in ev}))
    df = df.copy()
    df["cut_dist_min"] = cut
    df["n_distinct_len"] = ndl
    return df


def load_snv_noise(path, positions, slop=SNV_NOISE_SLOP):
    """DRAGEN systematic-noise BED -> {(chrom, pos): (max_noise, n_donors, alleles)}.

    Columns: chrom, start, end, mean_noise, max_noise, alleles, n_donors. `end` is the 1-based
    position, matching the `end` column of the caller's table.

    Despite the "snv" in the filename this panel is NOT SNV-only: the allele column carries D
    (1,299,670 records) and I (892,925) codes, so it covers indel noise too. Measured against the
    1,498 gated rows of the 32-sample cohort, a bare interval hit flags 33.1% of artifacts but also
    1.2% of real on-target edits -- 27x enrichment, but not clean enough to use as-is.

    The panel's own fields separate the two cases, which is why this returns them rather than a
    boolean:

        IKZF2  chr2:213,147,790  on-target,  n_donors=1,  alleles 'G',   max 0.027  <- coincidence
        B2M    chr1:28,580,333   off-target, n_donors=17, alleles 'C,D', max 0.078  <- real warning

    The file is ~1 GB and is not shipped with a tabix index, so this streams it once (~2-3 min).
    Only positions in `positions` are retained, so memory stays small.
    """
    want = collections.defaultdict(set)
    for c, e in positions:
        for d in range(-slop, slop + 1):
            want[c].add(int(e) + d)
    hits = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            c = p[0]
            if c not in want:
                continue
            try:
                pos = int(p[2])
            except ValueError:
                continue
            if pos not in want[c]:
                continue
            try:
                prev = hits.get((c, pos))
                rec = (float(p[4]), int(p[6]), p[5])
                # keep the strongest record if several land on one position
                if prev is None or rec[1] > prev[1]:
                    hits[(c, pos)] = rec
            except ValueError:
                continue
    return hits


def snv_noise_mask(hits, chroms, ends, min_donors, slop=SNV_NOISE_SLOP):
    """True where a site sits on a recurrent indel-noise locus.

    Two conditions beyond a plain interval hit, both derived from the measurement above:
      * the locus recurs in at least `min_donors` panel donors -- drops the n=1 coincidences
      * the noise involves an insertion or deletion -- an SNV-only locus says nothing about an
        indel call at the same coordinate
    """
    out = []
    for c, e in zip(chroms, ends):
        flag = False
        for d in range(-slop, slop + 1):
            rec = hits.get((c, int(e) + d))
            if rec and rec[1] >= min_donors and ("D" in rec[2] or "I" in rec[2]):
                flag = True
                break
        out.append(flag)
    return out


def guide_of(sample):
    """NS0011-ABTB1 / CART_NS0027-B2M_1 -> ABTB1 / B2M"""
    return re.sub(r'_\d+$', '', re.sub(r'^(CART_)?NS\d+-', '', str(sample)))


class RepeatIndex:
    """Point-in-interval lookup over one or more BED files.

    Intervals are held as per-chromosome sorted (start, end) arrays and queried with bisect,
    which keeps a 120 MB RepeatMasker BED usable without pulling in pybedtools.
    """

    def __init__(self, paths):
        self.iv = {}
        self.n = 0
        for p in paths:
            op = gzip.open if str(p).endswith(".gz") else open
            with op(p, "rt") as fh:
                for line in fh:
                    if line.startswith(("#", "track", "browser")):
                        continue
                    f = line.split("\t")
                    if len(f) < 3:
                        continue
                    try:
                        self.iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
                    except ValueError:
                        continue
                    self.n += 1
        for c in self.iv:
            self.iv[c].sort()

    def __contains__(self, key):
        chrom, pos = key
        a = self.iv.get(chrom)
        if not a:
            return False
        # rightmost interval whose start <= pos
        i = bisect.bisect_right(a, (pos, float("inf"))) - 1
        return i >= 0 and a[i][0] <= pos <= a[i][1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="*.offtarget_analysis.tsv")
    ap.add_argument("-o", "--out", required=True, help="filtered review queue (TSV)")
    ap.add_argument("--pon", help="panel-of-normals TSV from build_offtarget_pon.py "
                                  "(single-guide safe; preferred over cross-guide recurrence)")
    ap.add_argument("--min-reads", type=int, default=GATE_READS)
    ap.add_argument("--min-vaf", type=float, default=GATE_VAF)
    ap.add_argument("--max-cut-dist", type=int, default=MAX_CUT_DIST)
    ap.add_argument("--min-distinct-len", type=int, default=MIN_DISTINCT_LEN)
    ap.add_argument("--max-guides", type=int, default=MAX_GUIDES,
                    help="off-target sites appearing in >= this many distinct guides are dropped "
                         "(fallback only, ignored when --pon is given)")
    ap.add_argument("--max-control-vaf", type=float, default=MAX_CONTROL_VAF,
                    help="drop sites whose matched control carries indels at >= this fraction")
    ap.add_argument("--snv-noise", metavar="BED",
                    help="DRAGEN systematic-noise BED (rule 6). Supplements the run's own PoN with "
                         "an external panel; on-target sites are exempt. Streams a ~1 GB file.")
    ap.add_argument("--snv-noise-min-donors", type=int, default=SNV_NOISE_MIN_DONORS,
                    help="a noise locus must recur in at least this many panel donors to count "
                         "(default 3; at 1 it starts flagging real on-target edits)")
    ap.add_argument("--repeats", nargs="*", default=None, metavar="BED",
                    help="repeat-annotation BED(s); off-target sites inside a repeat are dropped. "
                         "Needs no controls, cohort or guide context, so this is the rule-4 "
                         "substitute for a brand-new guide with nothing else available.")
    ap.add_argument("--noise-model", default="off", choices=["off", "matched", "loo", "both"],
                    help="REPLACE rule 4 with a beta-binomial test against a control-derived "
                         "background (bin/noise_model.py). 'matched' uses the sample's own "
                         "unedited control and needs no cohort, which is the point: it matches "
                         "the PoN's measured performance (64/64 recall, 9 rejects, precision "
                         "0.877) without one. Default off, so existing runs are unchanged.")
    ap.add_argument("--aq-min", type=float, default=AQ_MIN,
                    help=f"drop sites scoring below this AQ under --noise-model (default "
                         f"{AQ_MIN}). With the depth floor on, 3-8 all reproduce the PoN exactly; "
                         f"10 starts costing confirmed edits")
    ap.add_argument("--no-depth-floor", dest="depth_floor", action="store_false",
                    help="disable the 1/control_depth floor on the posterior background "
                         "(see apply_depth_floor in noise_model.py -- leaving it off lets a clean "
                         "control claim a background it has no power to support)")
    ap.add_argument("--strict-fallback", action="store_true",
                    help="exit non-zero instead of warning when rule 4 has no usable source "
                         "(no --pon, no --noise-model, and only one guide in this invocation)")
    ap.add_argument("--keep-all", action="store_true",
                    help="emit every gated row with a why_dropped column instead of filtering")
    a = ap.parse_args()

    frames = []
    for p in a.inputs:
        d = pd.read_csv(p, sep="\t")
        d["sample_name"] = os.path.basename(p).split(".offtarget")[0]
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    gate = (df.indel_reads >= a.min_reads) & (df.indel_fraction >= a.min_vaf)
    q = add_features(df[gate])
    if q.empty:
        pd.DataFrame(columns=list(df.columns)).to_csv(a.out, sep="\t", index=False)
        print("no rows cleared the gate; wrote an empty queue")
        return
    q["guide"] = q.sample_name.map(guide_of)
    q = q.join(q.groupby(["chrom", "start"]).guide.nunique().rename("n_guides"),
               on=["chrom", "start"])
    n_guides_total = q.guide.nunique()

    # Rule 1. Control VAF, not a bare count: a couple of stray reads in a 60x normal is noise,
    # a third of them is germline.
    ctrl_vaf = q.control_indel_reads.fillna(0) / q.control_reads.replace(0, np.nan)
    q["control_vaf"] = ctrl_vaf.fillna(0).round(4)
    ctrl = q.control_vaf >= a.max_control_vaf

    # Rules 2 and 3, shape at the cut site.
    far = ~(q.cut_dist_min <= a.max_cut_dist)
    mono = q.n_distinct_len < a.min_distinct_len

    # Rule 4. On-target sites are exempt from both forms: they are shared by design.
    is_off = q.get("is_target", pd.Series(0, index=q.index)) == 0
    pon_coverage = None
    if a.noise_model != "off":
        # Rule 4 as a statistical test rather than a blacklist. Imported lazily so that a run
        # without --noise-model never needs scipy.
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from noise_model import (fit_global_prior, control_posterior, apply_depth_floor,
                                 aq_from_sf)
        from scipy.stats import betabinom

        a0, b0 = fit_global_prior(df.control_indel_reads.fillna(0).values,
                                  df.control_reads.fillna(0).values)
        alpha, beta, _nsup, dep = control_posterior(df, q, a.noise_model, a0, b0)
        n_floored = 0
        if a.depth_floor:
            alpha, beta, n_floored = apply_depth_floor(alpha, beta, dep)
        q["bg_rate"] = (alpha / (alpha + beta)).round(6)
        q["AQ"] = aq_from_sf(betabinom.sf(q.indel_reads.values.astype(int) - 1,
                                          q.total_reads.values.astype(int), alpha, beta)).round(2)
        known_bad = q.AQ < a.aq_min
        rule4 = f"beta-binomial vs {a.noise_model} control (AQ<{a.aq_min:g})"
        print(f"noise model           : prior Beta({a0:.4g},{b0:.4g}), "
              f"depth floor {'on' if a.depth_floor else 'OFF'} ({n_floored} rows raised)")
    elif a.pon:
        pon = pd.read_csv(a.pon, sep="\t")
        # A PoN is specific to the guide panel it was built from. Applied to a different guide's
        # sites it silently matches nothing, which is indistinguishable from a clean result --
        # so measure how much of this queue the PoN actually covers before trusting it.
        universe = set(zip(pon.chrom, pon.pos))
        bad = set(zip(pon.loc[pon.get("blacklisted", 1) == 1, "chrom"],
                      pon.loc[pon.get("blacklisted", 1) == 1, "pos"]))
        in_universe = [(c, int(p)) in universe for c, p in zip(q.chrom, q.end)]
        pon_coverage = float(np.mean(in_universe)) if len(q) else 0.0
        known_bad = is_off & pd.Series([(c, int(p)) in bad for c, p in zip(q.chrom, q.end)],
                                       index=q.index)
        rule4 = "panel of normals"
    else:
        known_bad = is_off & (q.n_guides >= a.max_guides)
        rule4 = "cross-guide recurrence"

    # Rule 5. Repeat context. Measured on the CAR-T cohort: RepeatMasker flags 100 of 177
    # artifacts and 0 of 61 real edits, tandem repeats 54 and 0 -- the homology-built off-target
    # panel is enriched for exactly the sequence aligners misplace indels in. On-target exempt.
    if a.repeats:
        rep_idx = RepeatIndex(a.repeats)
        in_repeat = is_off & pd.Series([(c, int(p)) in rep_idx for c, p in zip(q.chrom, q.end)],
                                       index=q.index)
    else:
        rep_idx = None
        in_repeat = pd.Series(False, index=q.index)

    # Rule 6. External systematic-noise panel. Placed LAST in the np.select order deliberately:
    # first match wins, so appending it here cannot change how any previously-dropped row is
    # attributed -- it can only claim rows the five existing rules kept.
    if a.snv_noise:
        _hits = load_snv_noise(a.snv_noise, set(zip(q.chrom, q.end)))
        snv_noise = is_off & pd.Series(
            snv_noise_mask(_hits, q.chrom, q.end, a.snv_noise_min_donors), index=q.index)
    else:
        snv_noise = pd.Series(False, index=q.index)

    rule4_label = (f"indistinguishable from control noise (AQ<{a.aq_min:g})"
                   if a.noise_model != "off" else "known-bad site (%s)" % rule4)
    q["why_dropped"] = np.select(
        [ctrl, far, mono, known_bad, in_repeat, snv_noise],
        ["germline (present in control)", "far from PAM", "single indel length",
         rule4_label, "repeat region",
         "systematic noise (external panel)"], default="")
    keep = q.why_dropped == ""

    (q if a.keep_all else q[keep]).to_csv(a.out, sep="\t", index=False)

    print(f"input rows            : {len(df)}")
    print(f"cleared the gate      : {len(q)}   (reads>={a.min_reads}, VAF>={a.min_vaf})")
    print(f"samples / guides      : {q.sample_name.nunique()} / {n_guides_total}")
    print(f"rule 4 source         : {rule4}"
          + (f"  (covers {pon_coverage:.0%} of this queue)" if pon_coverage is not None else ""))
    for lab, m in [("germline (in control)", ctrl), ("far from PAM", far & ~ctrl),
                   ("single indel length", mono & ~ctrl & ~far),
                   ("known-bad site", known_bad & ~ctrl & ~far & ~mono),
                   ("repeat region", in_repeat & ~ctrl & ~far & ~mono & ~known_bad),
                   ("systematic noise (panel)",
                    snv_noise & ~ctrl & ~far & ~mono & ~known_bad & ~in_repeat)]:
        print(f"  dropped, {lab:24s}: {int(m.sum())}")
    if rep_idx is not None:
        print(f"repeat annotation     : {rep_idx.n} intervals from {len(a.repeats)} file(s)")
    if a.keep_all:
        print(f"REVIEW QUEUE          : {int(keep.sum())}  "
              f"(all {len(q)} gated rows written with why_dropped -> {a.out})")
    else:
        print(f"REVIEW QUEUE          : {int(keep.sum())}  -> {a.out}")

    # --- guards against the two ways this silently degrades ---------------------------------
    if (q.control_indel_reads.fillna(0) == 0).all():
        print("\nWARNING: control_indel_reads is 0 at every site, so the germline rule did "
              "nothing.\n         These inputs predate the find_edited_reads.py control-reporting "
              "fix, which\n         summed control support *after* the -x filter removed it. "
              "Re-run the caller;\n         the germline rule is the single biggest filter "
              "(238 -> 84 on its own).", file=sys.stderr)

    if pon_coverage is not None and pon_coverage < 0.5:
        print(f"\nWARNING: the supplied PoN covers only {pon_coverage:.0%} of the sites in this "
              f"queue.\n         A panel of normals is specific to the guide panel it was built "
              f"from -- this one\n         was almost certainly built for a different guide, so "
              f"rule 4 is barely firing.\n         Score your unedited sample(s) against THIS "
              f"guide's target file and rebuild:\n"
              f"             build_offtarget_pon.py <unedited>.offtarget_analysis.tsv -o pon.tsv",
              file=sys.stderr)

    # Rule 4 has three possible sources and one dangerous hole. The hole is a SILENT one: with no
    # PoN, no noise model and a single guide in the invocation, cross-guide recurrence cannot fire
    # and the output looks identical to a clean result -- so the user reads an unfiltered queue as
    # a precise one. Say so loudly, and let a caller make it fatal.
    if a.noise_model != "off":
        pass                                            # rule 4 is covered by the statistical test
    elif not a.pon and n_guides_total < 2:
        msg = (f"\n{'=' * 78}\n"
               f"WARNING: RULE 4 IS NOT ACTIVE. No --pon, no --noise-model, and only "
               f"{n_guides_total} guide in this\n"
               f"         invocation -- cross-guide recurrence needs >=2 guides passed TOGETHER,\n"
               f"         so it cannot fire. Known-bad sites are NOT being removed and this queue\n"
               f"         is less precise than it looks.\n"
               f"         Fix, in order of preference:\n"
               f"           --noise-model matched   (needs only this sample's own control)\n"
               f"           --pon pon.tsv           (bin/build_offtarget_pon.py)\n"
               f"           pass all guides in ONE invocation\n"
               f"{'=' * 78}")
        print(msg, file=sys.stderr)
        if a.strict_fallback:
            sys.exit("ERROR: --strict-fallback set and rule 4 has no usable source")
    elif not a.pon:
        print(f"\nNOTE: using cross-guide recurrence for rule 4 across {n_guides_total} guides. "
              f"A panel of\n      normals (--pon) or --noise-model matched is strictly better and "
              f"is not affected\n      by how many guides you pass in one invocation.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
