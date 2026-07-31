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

Measured on the 25-sample CAR-T WGS cohort against the corrected truth set (61 real edits,
177 artifacts):

    review queue 238 -> 63 rows, 61/61 real edits retained
    precision 0.256 -> 0.968, recall 1.000

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

usage:
  review_filter.py IN.tsv [IN2.tsv ...] -o queue.tsv [--pon assets/offtarget_pon.tsv]
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

GATE_READS, GATE_VAF = 10, 0.05
MAX_CUT_DIST, MIN_DISTINCT_LEN, MAX_GUIDES = 10, 3, 2
MAX_CONTROL_VAF = 0.05


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


def guide_of(sample):
    """NS0011-ABTB1 / CART_NS0027-B2M_1 -> ABTB1 / B2M"""
    return re.sub(r'_\d+$', '', re.sub(r'^(CART_)?NS\d+-', '', str(sample)))


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
    if a.pon:
        pon = pd.read_csv(a.pon, sep="\t")
        # A PoN is specific to the guide panel it was built from. Applied to a different guide's
        # sites it silently matches nothing, which is indistinguishable from a clean result --
        # so measure how much of this queue the PoN actually covers before trusting it.
        universe = set(zip(pon.chrom, pon.pos))
        bad = set(zip(pon.loc[pon.get("blacklisted", 1) == 1, "chrom"],
                      pon.loc[pon.get("blacklisted", 1) == 1, "pos"]))
        in_universe = [(c, int(p)) in universe for c, p in zip(q.chrom, q.end)]
        pon_coverage = float(np.mean(in_universe)) if len(q) else 0.0
        known_bad = is_off & [(c, int(p)) in bad for c, p in zip(q.chrom, q.end)]
        rule4 = "panel of normals"
    else:
        known_bad = is_off & (q.n_guides >= a.max_guides)
        rule4 = "cross-guide recurrence"

    q["why_dropped"] = np.select(
        [ctrl, far, mono, known_bad],
        ["germline (present in control)", "far from PAM", "single indel length",
         "known-bad site (%s)" % rule4], default="")
    keep = q.why_dropped == ""

    (q if a.keep_all else q[keep]).to_csv(a.out, sep="\t", index=False)

    print(f"input rows            : {len(df)}")
    print(f"cleared the gate      : {len(q)}   (reads>={a.min_reads}, VAF>={a.min_vaf})")
    print(f"samples / guides      : {q.sample_name.nunique()} / {n_guides_total}")
    print(f"rule 4 source         : {rule4}"
          + (f"  (covers {pon_coverage:.0%} of this queue)" if pon_coverage is not None else ""))
    for lab, m in [("germline (in control)", ctrl), ("far from PAM", far & ~ctrl),
                   ("single indel length", mono & ~ctrl & ~far),
                   ("known-bad site", known_bad & ~ctrl & ~far & ~mono)]:
        print(f"  dropped, {lab:24s}: {int(m.sum())}")
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

    if not a.pon and n_guides_total < 2:
        print(f"\nWARNING: no --pon given and only {n_guides_total} guide present, so rule 4 "
              f"cannot fire.\n         The filter has degraded to three rules. Build a panel of "
              f"normals from any\n         unedited samples you have "
              f"(bin/build_offtarget_pon.py) and pass --pon; it\n         uses no guide "
              f"information and works for a single guide.", file=sys.stderr)
    elif not a.pon:
        print(f"\nNOTE: using cross-guide recurrence for rule 4 across {n_guides_total} guides. "
              f"A panel of\n      normals (--pon) is strictly better and is not affected by how "
              f"many guides\n      you pass in one invocation.", file=sys.stderr)


if __name__ == "__main__":
    main()
