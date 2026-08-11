#!/usr/bin/env python3
"""Shortlist the breakend (BND) calls a human still needs to look at, and name what they are.

Takes find_edited_reads.py output (*.offtarget_analysis.tsv), expands the `bnd_info` column into one
row per junction, and applies the breakend analogue of the indel review filter. The interesting
output is not only the shortlist: on this cohort most surviving junctions are **multi-cut
deletions** -- two cuts from the same guide's target set with the intervening segment excised --
which the pipeline has always emitted and never labelled.

  1. enough support            -- >= 3 reads. Breakend support is thin (cohort median 1 read), so
                                  this is the single most discriminating cut available.
  2. near the cut              -- within 10 bp of a PAM position, same rule as indels. The caller
                                  has already capped this at 25 bp via -d/--max-mutation-distance.
  3. breakpoint not promiscuous-- a breakpoint partnering with many unrelated loci is an alignment
                                  hub, not a junction. On-target sites are EXEMPT (see below).
  4. matched control clean     -- NOT applied here. It is already applied per event upstream: the
                                  caller's -x/--max-in-control (default 0) drops any junction with
                                  control support before it is written to bnd_info.

Rule 4 only started working for breakends once the position handed to add_normal_counts was
corrected for BND rows (see docs/CALLER_INTEGRATION_LOG.md). Before that fix, control support at a
breakend was structurally unmeasurable and every junction survived it.

ON-TARGET SITES MUST BE EXEMPT FROM RULE 3. A real Cas9 cut generates junctions to many places, so
the true cut sites are among the most promiscuous breakpoints in the cohort -- measured on the
32-sample CAR-T set, 6 of the 10 breakpoints with >=5 distinct partners are the intended TRAC,
TRBC1, TRBC2 and B2M cut sites. Applying rule 3 without the exemption deletes the real edits.

CALIBRATION NOTE, worth knowing before tuning rule 3. On the 32-sample cohort, once rule 1 is
applied, the ONLY promiscuous breakpoint left is chr1:246,009,98x -- and the caller's control filter
now removes that one on its own. Rule 3 is therefore **currently non-binding**: it changes nothing
on this data. It is kept because it costs nothing, it is the only defence against an alignment hub
that happens to be absent from the matched control, and its column is worth reporting either way.
Do not read "0 dropped by rule 3" as the rule being broken.

usage:
  review_filter_bnd.py IN.tsv [IN2.tsv ...] -o bnd_queue.tsv [--keep-all]
"""
import argparse
import bisect
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd

MIN_READS = 3
MAX_CUT_DIST = 10
MAX_PARTNERS = 5
BIN_SIZE = 1000
ONTARGET_SLOP = 500
SV_NOISE_SLOP = 50


def parse_bnds(df, sample):
    """bnd_info -> one row per junction.

    Field layout per ';'-separated event, identical to indel_info:
      chrom|pos|chrom2|pos2|strands|ref|alt|distance|distance2|counts|control_alt_counts
    'distance' (index 7) is min(|pos - PAM_position|) over the site's PAM positions.
    """
    out = []
    for _, r in df[df.bnd_count > 0].iterrows():
        for ev in str(r.bnd_info).split(';'):
            f = ev.split('|')
            if len(f) < 11:
                continue
            try:
                out.append({
                    'sample_name': sample,
                    'chrom': f[0], 'pos': int(f[1]),
                    'chrom2': f[2], 'pos2': int(f[3]),
                    'strands': f[4],
                    'cut_dist': abs(int(f[7])),
                    'reads': int(f[9]),
                    'control_reads_at_event': int(f[10]),
                    'site_start': r.start, 'site_end': r.end,
                    'is_target': int(r.is_target),
                    'site_total_reads': r.total_reads,
                })
            except ValueError:
                continue
    return out


def load_sv_noise(path):
    """BEDPE -> {chrom: (starts, ends)} of merged breakpoint intervals, both ends pooled.

    Measured on the 32-sample CAR-T cohort against the 36 gated junctions:

        panel                                     real queue flagged   artifacts flagged
        WGS_hg38_v3.1.0    (311k records)              0 / 25              11 / 11
        IDPF_WGS_v3.0.0    (2.6M records)             25 / 25              11 / 11
        WGS_FF_Heme_v3.1.0 (2.2M records)             24 / 25              11 / 11

    Only the small WGS panel discriminates. IDPF and FF_Heme flag the real junctions at a HIGHER
    rate than the noise -- as a blacklist they preferentially delete findings, so they are not
    merely useless here but actively harmful. Slop matters too: at 0 or 50 bp the WGS panel is
    perfect, at 200 bp it starts flagging real junctions.
    """
    iv = collections.defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                iv[p[0]].append((int(p[1]), int(p[2])))
                iv[p[3]].append((int(p[4]), int(p[5])))
            except ValueError:
                continue
    merged = {}
    for c, v in iv.items():
        v.sort()
        out = []
        for s, e in v:
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[c] = ([x[0] for x in out], [x[1] for x in out])
    return merged


def in_noise(merged, chrom, pos, slop):
    if chrom not in merged:
        return False
    st, en = merged[chrom]
    i = bisect.bisect_right(st, pos + slop) - 1
    return i >= 0 and st[i] - slop <= pos <= en[i] + slop


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="*.offtarget_analysis.tsv")
    ap.add_argument("-o", "--out", required=True, help="breakend review queue (TSV)")
    ap.add_argument("--min-reads", type=int, default=MIN_READS)
    ap.add_argument("--max-cut-dist", type=int, default=MAX_CUT_DIST)
    ap.add_argument("--max-partners", type=int, default=MAX_PARTNERS,
                    help="breakpoints with >= this many distinct partner loci are alignment hubs")
    ap.add_argument("--bin-size", type=int, default=BIN_SIZE,
                    help="breakpoints are binned at this resolution before counting partners")
    ap.add_argument("--sv-noise", metavar="BEDPE",
                    help="DRAGEN systematic-noise SV panel (BEDPE, optionally gzipped). Use "
                         "WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz -- see the note below on "
                         "which panels are safe.")
    ap.add_argument("--sv-noise-slop", type=int, default=SV_NOISE_SLOP,
                    help="bp of tolerance when matching a breakpoint to a noise interval "
                         "(default 50; do NOT raise to 200, it starts hitting real junctions)")
    ap.add_argument("--keep-all", action="store_true",
                    help="write every gated junction with a why_dropped column, for auditing")
    a = ap.parse_args()

    events, ontarget = [], {}
    for p in a.inputs:
        d = pd.read_csv(p, sep="\t")
        s = os.path.basename(p).split(".offtarget")[0]
        events += parse_bnds(d, s)
        # The intended cut sites for this sample, used both to exempt them from the promiscuity
        # rule and to recognise a junction whose FAR end also lands on a cut site.
        ontarget[s] = d.loc[d.is_target == 1, ['chrom', 'start', 'end']].values.tolist()

    if not events:
        pd.DataFrame(columns=['sample_name', 'chrom', 'pos', 'chrom2', 'pos2']).to_csv(
            a.out, sep="\t", index=False)
        print("no breakend events in the input; wrote an empty queue")
        return

    b = pd.DataFrame(events)
    n_in = len(b)

    # Promiscuity: distinct partner loci per binned breakpoint, counted across everything supplied
    # in this invocation. Cohort-wide is stronger than per-sample -- an alignment hub recurs.
    b['bin'] = b.chrom + ':' + (b.pos // a.bin_size).astype(str)
    b['partner_bin'] = b.chrom2 + ':' + (b.pos2 // a.bin_size).astype(str)
    b = b.join(b.groupby('bin').partner_bin.nunique().rename('n_partners'), on='bin')

    b['span'] = np.where(b.chrom == b.chrom2, (b.pos2 - b.pos).abs(), pd.NA)
    b['interchromosomal'] = (b.chrom != b.chrom2).astype(int)

    def far_end_on_target(row):
        for c, s, e in ontarget.get(row.sample_name, []):
            if c == row.chrom2 and (s - ONTARGET_SLOP) <= row.pos2 <= (e + ONTARGET_SLOP):
                return 1
        return 0

    b['far_end_on_target'] = b.apply(far_end_on_target, axis=1)

    gate = b.reads >= a.min_reads
    q = b[gate].copy()
    if q.empty:
        b.head(0).to_csv(a.out, sep="\t", index=False)
        print(f"input junctions       : {n_in}\nno junction cleared >= {a.min_reads} reads")
        return

    noise = pd.Series(False, index=q.index)
    if a.sv_noise:
        merged = load_sv_noise(a.sv_noise)
        # Either end landing in a known-noise interval condemns the junction. On-target sites are
        # exempt for the same reason they are exempt from the promiscuity rule.
        noise = pd.Series(
            [(in_noise(merged, c, p, a.sv_noise_slop) or in_noise(merged, c2, p2, a.sv_noise_slop))
             for c, p, c2, p2 in zip(q.chrom, q.pos, q.chrom2, q.pos2)],
            index=q.index) & (q.is_target == 0)

    far = ~(q.cut_dist <= a.max_cut_dist)
    # On-target breakpoints are the most promiscuous in the cohort by construction -- a real cut
    # throws junctions everywhere. Exempt them, exactly as rules 4 and 5 do for indels.
    hub = (q.n_partners >= a.max_partners) & (q.is_target == 0)

    q["why_dropped"] = np.select([far, hub, noise],
                                 ["far from PAM", "promiscuous breakpoint",
                                  "DRAGEN systematic noise"], default="")
    keep = q.why_dropped == ""

    # Classify what survives. `is_target` describes the NEAR end -- the site the junction was found
    # at -- so a junction anchored at a cut site is an on-target editing outcome even when its far
    # end is nowhere in particular. Only a junction whose near end is not a cut site is off-target.
    #
    #   both ends at cut sites, same chrom -> the segment between two cuts was excised
    #   near end at a cut site, same chrom -> one cut, resected and joined further out
    #   near end at a cut site, other chrom -> one cut joined to another chromosome
    #   near end not a cut site            -> genuinely off-target
    at_cut = q.is_target == 1
    q["call"] = np.select(
        [keep & at_cut & (q.far_end_on_target == 1) & (q.interchromosomal == 0),
         keep & at_cut & (q.interchromosomal == 0),
         keep & at_cut,
         keep],
        ["multi-cut deletion",
         "deletion at cut site",
         "translocation at cut site",
         "off-target junction"],
        default="")

    (q if a.keep_all else q[keep]).to_csv(a.out, sep="\t", index=False)

    print(f"input junctions       : {n_in}")
    print(f"cleared the gate      : {len(q)}   (reads>={a.min_reads})")
    print(f"samples               : {b.sample_name.nunique()}")
    for lab, m in [("far from PAM", far), ("promiscuous breakpoint", hub & ~far),
                   ("DRAGEN systematic noise", noise & ~far & ~hub)]:
        print(f"  dropped, {lab:24s}: {int(m.sum())}")
    if int((hub & ~far).sum()) == 0:
        print("           (rule 3 is currently non-binding on this data -- see the module docstring)")
    print(f"BREAKEND QUEUE        : {int(keep.sum())}  -> {a.out}")
    for lab in ("multi-cut deletion", "deletion at cut site", "translocation at cut site",
                "off-target junction"):
        sub = q[keep & (q.call == lab)]
        if len(sub):
            extra = ""
            if lab.endswith("deletion") and sub.span.notna().any():
                sp = pd.to_numeric(sub.span, errors="coerce").dropna()
                extra = f"   span {int(sp.min()):,}-{int(sp.max()):,} bp"
            print(f"    {lab:22s}: {len(sub):4d} in {sub.sample_name.nunique()} sample(s){extra}")

    if int((q.control_reads_at_event > 0).sum()) == 0:
        print("\nNote: control support is 0 on every junction, as expected -- the caller's "
              "-x/--max-in-control\n      filter removes control-supported junctions before they "
              "reach bnd_info. That is the\n      matched-control rule, applied upstream, not a "
              "missing check.", file=sys.stderr)


if __name__ == "__main__":
    main()
