#!/usr/bin/env python3
"""
pileup_snapshot.py — render an IGV-style read-pileup image at a locus straight
from the CRAM, so a reviewer confirms an edit by eye without opening IGV. Designed
to be a PIPELINE STEP: point it at score.py's worklist and it emits one PNG per
candidate (auto-review packet).

Reads are packed into rows (IGV-style), aligned blocks drawn as grey bars,
deletions as red gaps, insertions as purple ticks, soft-clips as faint blue. The
candidate position is a dashed vertical line. Title carries the call metadata.
"""
import os, sys, argparse
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REF = "/storage2/fs1/dspencer/Active/spencerlab/rotating_students/btyler/seqs/hg38_mgi_patch.fa"

# CIGAR: 0=M 1=I 2=D 3=N 4=S 5=H 7== 8=X
REF_CONSUME = {0, 2, 3, 7, 8}
QUERY_CONSUME = {0, 1, 4, 7, 8}


def read_blocks(read):
    """Return (aligned_blocks, deletions, insertions, softclips) in ref coords.

    Soft clips are (anchor, length, is_leading). A LEADING clip hangs off the left
    edge of the alignment and must be drawn leftwards from reference_start; a
    TRAILING clip hangs off the right edge and is drawn rightwards from the end.
    Getting this right matters most at a breakend, where a stack of clips all
    terminating on one base is the whole signal.
    """
    blocks, dels, ins, soft = [], [], [], []
    ref = read.reference_start
    cig = read.cigartuples or []
    for i, (op, ln) in enumerate(cig):
        if op in (0, 7, 8):                       # match/mismatch
            blocks.append((ref, ref + ln)); ref += ln
        elif op == 2 or op == 3:                  # deletion / skip
            dels.append((ref, ref + ln)); ref += ln
        elif op == 1:                             # insertion (0-width in ref)
            ins.append((ref, ln))
        elif op == 4:                             # soft clip
            soft.append((ref, ln, i == 0))
    return blocks, dels, ins, soft


def pack_rows(reads, pad=3):
    """Greedy IGV-style row packing by reference start."""
    rows = []   # list of last-end per row
    placed = []
    for r in sorted(reads, key=lambda x: x.reference_start):
        for ri, end in enumerate(rows):
            if r.reference_start > end + pad:
                rows[ri] = r.reference_end; placed.append((ri, r)); break
        else:
            rows.append(r.reference_end); placed.append((len(rows) - 1, r))
    return placed, len(rows)


def snapshot(ax, cram, chrom, pos, ref=REF, window=60, max_rows=120, title="",
             *, keep_supplementary=False, highlight=None, clip_cap=8, mark=None):
    """Draw an IGV-style pileup at one locus into `ax`.

    The keyword-only arguments exist for the breakend view (bnd_snapshots.py) and
    default to the historical behaviour, so the indel callers are unaffected:

      keep_supplementary  keep split-read segments. A junction's supporting reads are
                          usually supplementary alignments, so the default filter hides
                          exactly the evidence a breakend figure is meant to show.
      highlight           set of read names to draw in green instead of strand-tint.
      clip_cap            max drawn soft-clip length; None draws the true length, which
                          keeps clip length readable as a cue at a breakpoint.
      mark                extra positions to mark with a faint dashed vertical, e.g. the
                          jittered breakpoint calls belonging to the same junction.
    """
    highlight = highlight or set()
    bam = pysam.AlignmentFile(cram, "rc", reference_filename=ref)
    lo, hi = pos - window, pos + window

    def drop(r):
        if r.is_unmapped or r.is_secondary or r.is_duplicate:
            return True
        return r.is_supplementary and not keep_supplementary

    reads = [r for r in bam.fetch(chrom, max(0, lo), hi) if not drop(r)]
    placed, nrows = pack_rows(reads)
    n_indel = 0
    for ri, r in placed:
        if ri > max_rows:
            continue
        y = -ri
        blocks, dels, ins, soft = read_blocks(r)
        fwd = not r.is_reverse
        hot = r.query_name in highlight
        if hot:
            col, clipcol = ("#31a354" if fwd else "#a1d99b"), "#006d2c"
        else:
            col, clipcol = ("#8fb3d9" if fwd else "#d9a68f"), "#bcd4ec"
        for a, b in blocks:
            ax.add_patch(Rectangle((a, y - 0.4), b - a, 0.8, color=col, lw=0))
        for a, b in dels:                           # deletion = red span
            if a <= hi and b >= lo:
                ax.plot([a, b], [y, y], color="#cc2b2b", lw=1.6, solid_capstyle="butt")
                if lo <= a <= hi:
                    n_indel += 1
        for a, ln in ins:                           # insertion = purple tick
            if lo <= a <= hi:
                ax.plot([a, a], [y - 0.45, y + 0.45], color="#7a3fbf", lw=1.4)
                n_indel += 1
        for a, ln, leading in soft:                 # soft clip = stub off the read edge
            w = ln if clip_cap is None else min(ln, clip_cap)
            ax.add_patch(Rectangle((a - w if leading else a, y - 0.4), w, 0.8,
                                   color=clipcol, lw=0, alpha=0.9 if hot else 0.6))
    for m in (mark or ()):
        if m != pos and lo <= m <= hi:
            ax.axvline(m, color="k", ls=":", lw=0.6, alpha=0.35)
    ax.axvline(pos, color="k", ls="--", lw=0.8, alpha=0.7)
    ax.set_xlim(lo, hi); ax.set_ylim(-min(nrows, max_rows) - 1, 1)
    ax.set_yticks([]); ax.set_xlabel(f"{chrom}:{pos:,}")
    ax.set_title(title, fontsize=9)
    bam.close()
    return nrows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cram", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--window", type=int, default=60)
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    fig, ax = plt.subplots(figsize=(9, 6))
    n = snapshot(ax, args.cram, args.chrom, args.pos, args.ref, args.window, title=args.title)
    fig.tight_layout(); fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({n} read rows)")


if __name__ == "__main__":
    main()
