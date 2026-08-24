#!/usr/bin/env python3
"""
make_igv_session.py — a navigable index + IGV batch script for live review.

The problem this solves
-----------------------
The review material is spread across three coordinate conventions that do not line up:

  * the hotspot table (`training.tsv`) keys a site by its **predicted cut site**
  * the genome-wide worklist and the rendered snapshots key the same event by the
    **observed VCF indel position** — e.g. `rank005_..._chr7_26760821.png` for the hotspot at
    `chr7:26760750`, 71 bp away
  * the tagged BAMs are just BAMs; you have to know where to look

So there is no way to answer "show me the evidence for hotspot X" without doing the arithmetic by
hand, and three of the five strongest candidates have no snapshot at all. This emits one index
that carries all three, plus an IGV batch script so sites can be jumped to live.

Ordering is by `cut_dist` ascending, deliberately: an indel far from the predicted cut is probably
an artifact regardless of how good its VAF looks, so the defensible sites come first. The
highest-ECS-VAF candidate on AAVS1 (`chr7:26760750`, ecs_if 0.302) sits 71 bp out and sorts near
the bottom, which is the correct outcome.

Usage
-----
  make_igv_session.py --outdir results_offtarget_aavs1_igv/offtarget \
                      --out-prefix results_demo/igv_review
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_recall_aavs1 import load_truth, DEFAULT_TP, DEFAULT_TN, SLACK, _near

WINDOW = 150          # bp each side, matches the tagged-BAM window
SNAP_RE = re.compile(r"^rank(\d+)_(.+?)_(chr[^_]+)_(\d+)\.png$")


def index_snapshots(snap_dir):
    """(chrom, pos) -> filename, from the rank###_<sample>_<chrom>_<pos>.png naming."""
    out = {}
    for path in sorted(glob.glob(os.path.join(snap_dir, "*.png"))):
        m = SNAP_RE.match(os.path.basename(path))
        if m:
            out.setdefault((m.group(3), int(m.group(4))), os.path.basename(path))
    return out


def nearest_snapshot(snaps, chrom, pos, tol=120):
    """Snapshots are keyed by VCF position, sites by cut site; allow real slack."""
    best, best_d = None, None
    for (c, p), name in snaps.items():
        if c != chrom:
            continue
        d = abs(p - pos)
        if d <= tol and (best_d is None or d < best_d):
            best, best_d = name, d
    return best, best_d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True,
                    help="the run's offtarget/ directory (training.tsv, snapshots/, *_tagged.bam)")
    ap.add_argument("--out-prefix", required=True, help="prefix for the .csv and .bat outputs")
    ap.add_argument("--tp", default=DEFAULT_TP)
    ap.add_argument("--tn", default=DEFAULT_TN)
    ap.add_argument("--genome", default="hg38",
                    help="IGV genome id or a path to the FASTA used for the run")
    ap.add_argument("--window", type=int, default=WINDOW)
    args = ap.parse_args()

    tp_set, tn_set = load_truth(args.tp, args.tn)
    t = pd.read_csv(os.path.join(args.outdir, "training.tsv"), sep="\t", low_memory=False)
    t["called"] = t["verdict"].astype(str).str.contains("LIKELY EDIT", na=False)
    if "cut_dist" not in t.columns:
        t["cut_dist"] = (t["modal_pos"] - t["start"]).abs()

    t["curated"] = ["TP" if _near(c, p, tp_set, SLACK) else
                    ("TN" if _near(c, p, tn_set, SLACK) else "uncurated")
                    for c, p in zip(t["chrom"], t["start"])]

    # everything worth opening: the confirmed edits plus every call
    sel = t[t["called"] | (t["curated"] == "TP")].copy()

    snaps = index_snapshots(os.path.join(args.outdir, "snapshots"))
    ecs_bams = sorted(glob.glob(os.path.join(args.outdir, "*.tagged.bam")))
    wgs_bams = sorted(glob.glob(os.path.join(args.outdir, "**", "*wgs_tagged.bam"),
                                recursive=True))

    rows = []
    for _, r in sel.iterrows():
        obs = int(r["modal_pos"]) if pd.notna(r.get("modal_pos")) else int(r["start"])
        snap, snap_off = nearest_snapshot(snaps, str(r["chrom"]), int(r["start"]))
        rows.append({
            "curated": r["curated"],
            "sample": r.get("sample"),
            "hotspot": f"{r['chrom']}:{int(r['start'])}",
            "observed_indel": f"{r['chrom']}:{obs}",
            "cut_dist": r["cut_dist"],
            "modal_len": r.get("modal_len"),
            "indel_frac": r.get("indel_frac"),
            "ecs_if": r.get("ecs_if"),
            "score": r.get("score"),
            "verdict": r.get("verdict"),
            "igv_locus": f"{r['chrom']}:{max(1, obs - args.window)}-{obs + args.window}",
            "snapshot": snap or "MISSING",
            "snapshot_offset_bp": snap_off if snap else "",
        })
    idx = pd.DataFrame(rows)
    # curated edits first, then by how defensible the site is
    idx["_o"] = (idx["curated"] != "TP").astype(int)
    idx = idx.sort_values(["_o", "cut_dist"], na_position="last").drop(columns="_o")

    csv_path = args.out_prefix + ".csv"
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    idx.to_csv(csv_path, index=False)

    bat_path = args.out_prefix + ".bat"
    with open(bat_path, "w") as fh:
        fh.write("# IGV batch script — off-target review\n")
        fh.write("# Run: IGV > Tools > Run Batch Script, or `igv.sh -b this_file`\n")
        fh.write("# Reads are coloured by the XC tag: Edited_* vs Unedited_WT vs Skipped_*\n")
        fh.write("new\n")
        fh.write(f"genome {args.genome}\n")
        for b in ecs_bams + wgs_bams:
            fh.write(f"load {os.path.abspath(b)}\n")
        if not (ecs_bams or wgs_bams):
            fh.write("# WARNING: no tagged BAMs found in this run directory\n")
        fh.write("colorBy TAG XC\n")
        fh.write("maxPanelHeight 800\n")
        for _, r in idx.iterrows():
            fh.write(f"\n# {r['curated']:9s} {r['hotspot']:24s} cut_dist={r['cut_dist']} "
                     f"ecs_if={r['ecs_if']} verdict={r['verdict']}\n")
            fh.write(f"goto {r['igv_locus']}\n")

    n_tp = int((idx["curated"] == "TP").sum())
    print(f"\nwrote {csv_path}  ({len(idx)} sites: {n_tp} curated edits, "
          f"{int((idx['curated'] == 'TN').sum())} curated negatives, "
          f"{int((idx['curated'] == 'uncurated').sum())} uncurated)")
    print(f"wrote {bat_path}  (loads {len(ecs_bams)} ECS + {len(wgs_bams)} WGS tagged BAM(s))")
    print(f"snapshots matched: {int((idx['snapshot'] != 'MISSING').sum())} of {len(idx)}")
    print()
    show = ["curated", "hotspot", "cut_dist", "indel_frac", "ecs_if", "verdict", "snapshot"]
    with pd.option_context("display.width", 200, "display.max_rows", 40):
        print(idx[show].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
