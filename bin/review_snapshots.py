#!/usr/bin/env python3
"""Render a review packet: one IGV-style pileup PNG per site in the review queue.

The queue TSV names sites; this turns each one into a picture a human can adjudicate without
opening IGV. Each figure is TWO panels -- the edited sample on top, its matched unedited control
below, same locus, same scale. That layout is the point: most of the calls that survive the
filters are settled by looking at the control panel, because germline and shared alignment
artifacts appear in both while a real edit appears in only one.

Reads --cram-map, a TSV of `sample<TAB>edited_cram<TAB>control_cram` (absolute paths), matching
the convention score_hotspots.nf already uses for its cram list.

usage:
  review_snapshots.py --queue review_queue.tsv --cram-map map.tsv --fasta ref.fa --outdir snapshots
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# pileup_snapshot lives beside this script in bin/
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pileup_snapshot import snapshot  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queue", required=True, help="review_queue.tsv from review_filter.py")
    ap.add_argument("--cram-map", required=True,
                    help="TSV: sample<TAB>edited_cram<TAB>control_cram")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="snapshots")
    ap.add_argument("--window", type=int, default=80, help="bp either side of the site")
    ap.add_argument("--max-sites", type=int, default=500,
                    help="safety cap; a runaway queue should not render thousands of PNGs")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    q = pd.read_csv(a.queue, sep="\t")
    if q.empty:
        print("review queue is empty; nothing to render")
        return

    cmap = pd.read_csv(a.cram_map, sep="\t", header=None,
                       names=["sample", "edited", "control"]).set_index("sample")

    if len(q) > a.max_sites:
        print(f"WARNING: queue has {len(q)} sites, rendering only the first {a.max_sites}. "
              f"A queue this large usually means a filter is not firing -- check the "
              f"review_filter.py warnings.", file=sys.stderr)
        q = q.head(a.max_sites)

    made, skipped = 0, []
    for i, r in q.reset_index(drop=True).iterrows():
        s = r.get("sample_name")
        if s not in cmap.index:
            skipped.append(f"{s} (not in cram map)")
            continue
        ed, ct = cmap.loc[s, "edited"], cmap.loc[s, "control"]
        pos = int(r.get("end", r.get("start")))
        tag = f"{i + 1:03d}_{s}_{r.chrom}_{pos}".replace("/", "-")
        try:
            fig, ax = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
            snapshot(ax[0], ed, r.chrom, pos, ref=a.fasta, window=a.window,
                     title=f"EDITED  {s}   {int(r.indel_reads)}/{int(r.total_reads)} indel reads"
                           f"   VAF {float(r.indel_fraction):.3f}")
            snapshot(ax[1], ct, r.chrom, pos, ref=a.fasta, window=a.window,
                     title=f"MATCHED UNEDITED CONTROL   ({os.path.basename(str(ct))})")
            bits = [f"{r.chrom}:{pos:,}"]
            if pd.notna(r.get("cut_dist_min")):
                bits.append(f"{float(r.cut_dist_min):.0f} bp from PAM")
            if pd.notna(r.get("n_distinct_len")):
                bits.append(f"{int(r.n_distinct_len)} indel lengths")
            if pd.notna(r.get("control_vaf")):
                bits.append(f"control VAF {float(r.control_vaf):.3f}")
            fig.suptitle("REVIEW CANDIDATE\n" + "   |   ".join(bits),
                         fontsize=11, fontweight="bold")
            plt.tight_layout(rect=[0, 0, 1, 0.93])
            fig.savefig(os.path.join(a.outdir, tag + ".png"), dpi=130)
            plt.close(fig)
            made += 1
        except Exception as e:                      # one bad locus must not kill the packet
            plt.close("all")
            skipped.append(f"{s} {r.chrom}:{pos} ({type(e).__name__}: {e})")

    print(f"review sites   : {len(q)}")
    print(f"snapshots made : {made}  -> {a.outdir}/")
    if skipped:
        print(f"skipped        : {len(skipped)}", file=sys.stderr)
        for s in skipped[:10]:
            print(f"  {s}", file=sys.stderr)


if __name__ == "__main__":
    main()
