#!/usr/bin/env python3
"""
bnd_snapshots.py -- render a review packet for BREAKENDS, one figure per junction.

The breakend sibling of review_snapshots.py, and deliberately not the same picture.
Three things make a junction different from an indel site:

1. A junction has TWO loci. review_snapshots.py draws one, and its
   `int(r.get("end", r.get("start")))` idiom does not even parse the breakend schema
   (which is chrom/pos + chrom2/pos2). So this is a 2x2 grid: left and right breakpoint,
   each with the edited sample over its matched unedited control -- keeping the
   self-adjudicating pairing the indel packet is built on.

2. One junction is reported as several rows. bnd_review_queue.tsv carries each event
   from BOTH ends, at +/- a few bp of jitter, under both strand orientations. On the
   reference cohort that is 25 rows for 8 real junctions, and rendering per row gives
   one sample five near-identical pictures. We collapse on (sample, {bin, partner_bin})
   -- columns review_filter_bnd.py already computes -- and sum the read support, because
   the per-row count understates the event by up to 6x (ARID4A: rows say 3-4 reads, the
   junction has 18).

3. The evidence is split reads, which the default pileup hides. A breakend's signature
   is a stack of soft clips terminating on one base whose clipped portion maps to the
   partner locus. Those are usually supplementary alignments, and pileup_snapshot's
   default filter drops them. We pass keep_supplementary=True and colour the reads whose
   SA tag points at the partner green, so the junction is visible rather than inferred.

On SA parsing: find_edited_reads.py has a far more capable SA handler
(get_sa_indel_vcf), but it reconstructs alleles and needs a fasta, target positions and
a target index, and importing that module pulls in edlib/joblib/pyranges/scipy. All we
need here is "does a segment of this read land near the partner locus", which is the
first two fields of the SA tag -- a SAM spec constant. Parsed inline on purpose.
"""
import os
import sys
import argparse
from argparse import RawDescriptionHelpFormatter
from collections import defaultdict

import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# pileup_snapshot lives beside this script in bin/
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pileup_snapshot import snapshot  # noqa: E402

BIN_SIZE = 1000          # must match review_filter_bnd.py's --bin-size
SA_SLOP = 1000           # how near the partner an SA segment must land to count


# --------------------------------------------------------------------------- junctions

def junction_bins(row):
    """The two 1 kb bins a junction connects, in canonical (sorted) order."""
    return tuple(sorted([str(row["bin"]), str(row["partner_bin"])]))


def end_of(group, wanted_bin):
    """Resolve one end of a junction to (chrom, representative_pos, all_positions).

    Rows report the junction from both directions, so a given bin appears as `bin` in
    some rows and `partner_bin` in others. The representative position is the one
    carrying the most read support; the rest are drawn as faint marks so the jitter is
    visible rather than hidden.
    """
    acc = defaultdict(int)
    chrom = None
    for _, r in group.iterrows():
        if str(r["bin"]) == wanted_bin:
            acc[int(r["pos"])] += int(r["reads"])
            chrom = str(r["chrom"])
        if str(r["partner_bin"]) == wanted_bin:
            acc[int(r["pos2"])] += int(r["reads"])
            chrom = str(r["chrom2"])
    if not acc:
        raise ValueError(f"no rows for bin {wanted_bin}")
    best = max(acc.items(), key=lambda kv: (kv[1], -kv[0]))[0]
    return chrom, best, sorted(acc)


def collapse(q):
    """Group the queue into junctions. Returns a list of dicts, one per junction."""
    out = []
    q = q.copy()
    q["_key"] = [junction_bins(r) for _, r in q.iterrows()]
    for (sample, key), g in q.groupby(["sample_name", "_key"], sort=False):
        b1, b2 = key
        try:
            c1, p1, marks1 = end_of(g, b1)
            c2, p2, marks2 = end_of(g, b2)
        except ValueError:
            continue
        # Canonical left/right so the figure reads the same way every time.
        if (c1, p1) > (c2, p2):
            (c1, p1, marks1), (c2, p2, marks2) = (c2, p2, marks2), (c1, p1, marks1)
        top = g.loc[g["reads"].astype(int).idxmax()]
        out.append({
            "sample": sample,
            "left": (c1, p1, marks1),
            "right": (c2, p2, marks2),
            "reads": int(g["reads"].astype(int).sum()),
            "n_rows": len(g),
            "depth": int(pd.to_numeric(g["site_total_reads"], errors="coerce").max()),
            "control": int(pd.to_numeric(g["control_reads_at_event"],
                                         errors="coerce").fillna(0).max()),
            "call": str(top.get("call", "") or "breakend"),
            "strands": ",".join(sorted(set(g["strands"].astype(str)))),
            "span": None if c1 != c2 else abs(p2 - p1),
            "interchrom": c1 != c2,
            "far_on_target": int(pd.to_numeric(g["far_end_on_target"],
                                               errors="coerce").fillna(0).max()),
        })
    return out


# ------------------------------------------------------------------------------ SA tag

def sa_supporters(cram, fasta, chrom, pos, window, pchrom, ppos, slop=SA_SLOP):
    """Read names near (chrom,pos) with an SA segment landing near (pchrom,ppos).

    SA tag format is `rname,pos,strand,CIGAR,mapQ,NM;` repeated -- we need fields 0-1.
    """
    names = set()
    try:
        bam = pysam.AlignmentFile(cram, "rc", reference_filename=fasta)
    except Exception:
        return names
    try:
        for r in bam.fetch(chrom, max(0, pos - window), pos + window):
            if r.is_unmapped or r.is_duplicate or not r.has_tag("SA"):
                continue
            for seg in str(r.get_tag("SA")).rstrip(";").split(";"):
                f = seg.split(",")
                if len(f) < 2:
                    continue
                try:
                    if f[0] == pchrom and abs(int(f[1]) - ppos) <= slop:
                        names.add(r.query_name)
                        break
                except ValueError:
                    continue
    finally:
        bam.close()
    return names


# ---------------------------------------------------------------------------- schematic

def draw_schematic(ax, j):
    """A to-scale cartoon of what happened, above the read panels."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    c1, p1, _ = j["left"]
    c2, p2, _ = j["right"]
    y = 0.55
    grey, cut, gone = "#4d4d4d", "#cc2b2b", "#bbbbbb"

    if j["interchrom"]:
        ax.add_patch(Rectangle((0.04, y - 0.05), 0.38, 0.10, color=grey))
        ax.add_patch(Rectangle((0.58, y - 0.05), 0.38, 0.10, color=grey))
        ax.plot([0.42], [y], marker="v", color=cut, ms=9)
        ax.plot([0.58], [y], marker="v", color=cut, ms=9)
        ax.annotate("", xy=(0.575, y - 0.20), xytext=(0.425, y - 0.20),
                    arrowprops=dict(arrowstyle="<->", color=cut, lw=1.4))
        ax.text(0.5, y - 0.42, "interchromosomal junction", ha="center",
                fontsize=8, color=cut)
        ax.text(0.23, y + 0.22, f"{c1}:{p1:,}", ha="center", fontsize=8)
        ax.text(0.77, y + 0.22, f"{c2}:{p2:,}", ha="center", fontsize=8)
        return

    xl, xr = 0.30, 0.70
    ax.add_patch(Rectangle((0.04, y - 0.05), xl - 0.04, 0.10, color=grey))
    ax.add_patch(Rectangle((xr, y - 0.05), 0.96 - xr, 0.10, color=grey))
    ax.add_patch(Rectangle((xl, y - 0.05), xr - xl, 0.10,
                           facecolor=gone, edgecolor=cut, hatch="///", lw=0.8))
    for x in (xl, xr):
        ax.plot([x], [y + 0.16], marker="v", color=cut, ms=9)
    span = j["span"]
    ax.text(0.5, y + 0.30, f"excised {span:,} bp" if span is not None else "excised",
            ha="center", fontsize=8, color=cut)
    ax.text(xl, y - 0.22, f"{c1}:{p1:,}", ha="center", fontsize=8)
    ax.text(xr, y - 0.22, f"{c2}:{p2:,}", ha="center", fontsize=8)
    ax.annotate("", xy=(xr, y - 0.36), xytext=(xl, y - 0.36),
                arrowprops=dict(arrowstyle="-", color=cut, lw=1.2,
                                connectionstyle="bar,fraction=-0.25"))
    ax.text(0.5, y - 0.52, "joined", ha="center", fontsize=8, color=cut)


# --------------------------------------------------------------------------------- main

def render(j, ed, ct, fasta, window, outpath):
    (c1, p1, m1), (c2, p2, m2) = j["left"], j["right"]
    hlL = sa_supporters(ed, fasta, c1, p1, window, c2, p2)
    hlR = sa_supporters(ed, fasta, c2, p2, window, c1, p1)

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.75, 2.2, 2.2],
                          hspace=0.52, wspace=0.10)
    draw_schematic(fig.add_subplot(gs[0, :]), j)

    axLe = fig.add_subplot(gs[1, 0])
    axRe = fig.add_subplot(gs[1, 1])
    # sharex per COLUMN only -- the two breakpoints have unrelated coordinates, so the
    # global sharex=True that review_snapshots.py uses would be meaningless here.
    axLc = fig.add_subplot(gs[2, 0], sharex=axLe)
    axRc = fig.add_subplot(gs[2, 1], sharex=axRe)

    kw = dict(ref=fasta, window=window, keep_supplementary=True, clip_cap=None)
    # NB the green count and the suptitle's junction-read count are different statistics
    # and will not agree: this one is every read with an SA segment within SA_SLOP of the
    # partner, the other is the caller's filtered support (MAPQ, cut distance, dedup).
    # Labelled explicitly so the figure does not appear to contradict itself.
    # Two-line titles: a single line does not fit a half-width column and the left and
    # right panel titles run into each other.
    snapshot(axLe, ed, c1, p1, highlight=hlL, mark=m1, **kw,
             title=f"EDITED  {j['sample']}   LEFT  {c1}:{p1:,}\n"
                   f"{len(hlL)} reads with SA at partner (green)")
    snapshot(axRe, ed, c2, p2, highlight=hlR, mark=m2, **kw,
             title=f"EDITED  {j['sample']}   RIGHT  {c2}:{p2:,}\n"
                   f"{len(hlR)} reads with SA at partner (green)")
    ctname = os.path.basename(str(ct))
    snapshot(axLc, ct, c1, p1, mark=m1, **kw,
             title=f"MATCHED UNEDITED CONTROL\n({ctname})")
    snapshot(axRc, ct, c2, p2, mark=m2, **kw,
             title=f"MATCHED UNEDITED CONTROL\n({ctname})")

    bits = [j["call"].upper()]
    if j["span"] is not None:
        bits.append(f"span {j['span']:,} bp")
    bits.append(f"{j['reads']} junction reads / {j['depth']}x")
    bits.append(f"control {j['control']}")
    bits.append(f"strands {j['strands']}")
    if j["n_rows"] > 1:
        bits.append(f"{j['n_rows']} queue rows")
    fig.suptitle("BREAKEND REVIEW CANDIDATE\n" + "   |   ".join(bits),
                 fontsize=11, fontweight="bold")
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument("--queue", required=True, help="bnd_review_queue.tsv from review_filter_bnd.py")
    ap.add_argument("--cram-map", required=True,
                    help="TSV: sample<TAB>edited_cram<TAB>control_cram")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="bnd_snapshots")
    ap.add_argument("--window", type=int, default=150,
                    help="bp either side of each breakpoint (default 150; a clip stack "
                         "needs more flank than an indel)")
    ap.add_argument("--max-junctions", type=int, default=200,
                    help="safety cap; a runaway queue should not render thousands of PNGs")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    q = pd.read_csv(a.queue, sep="\t")
    if q.empty:
        print("breakend review queue is empty; nothing to render")
        return

    junctions = collapse(q)
    print(f"queue rows      : {len(q)}")
    print(f"junctions       : {len(junctions)}")
    if len(junctions) > a.max_junctions:
        print(f"WARNING: {len(junctions)} junctions, rendering only the first "
              f"{a.max_junctions}. A set this large usually means a filter is not "
              f"firing -- check the review_filter_bnd.py warnings.", file=sys.stderr)
        junctions = junctions[:a.max_junctions]

    cmap = pd.read_csv(a.cram_map, sep="\t", header=None,
                       names=["sample", "edited", "control"]).set_index("sample")

    made, skipped = 0, []
    for i, j in enumerate(junctions):
        s = j["sample"]
        if s not in cmap.index:
            skipped.append(f"{s} (not in cram map)")
            continue
        c1, p1, _ = j["left"]
        c2, p2, _ = j["right"]
        tag = f"{i + 1:03d}_{s}_{c1}_{p1}_{c2}_{p2}".replace("/", "-")
        try:
            render(j, cmap.loc[s, "edited"], cmap.loc[s, "control"],
                   a.fasta, a.window, os.path.join(a.outdir, tag + ".png"))
            made += 1
        except Exception as e:            # one bad junction must not kill the packet
            plt.close("all")
            skipped.append(f"{s} {c1}:{p1}-{c2}:{p2} ({type(e).__name__}: {e})")

    print(f"snapshots made  : {made}  -> {a.outdir}/")
    if skipped:
        print(f"skipped         : {len(skipped)}", file=sys.stderr)
        for line in skipped[:10]:
            print(f"  {line}", file=sys.stderr)


if __name__ == "__main__":
    main()
