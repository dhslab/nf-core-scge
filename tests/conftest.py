"""Shared fixtures for the off-target glue-script tests.

These exercise the pure-pandas glue scripts (no pysam / CRAM / model needed):
  hotspot_to_table.py, join_training_table.py, recall_vs_vaf.py, reconcile_offtarget_report.py

Run with an interpreter that has pandas (the docker-scge container, or any env with
pandas). Scripts are invoked as subprocesses with cwd=tmp so the repo-root vendored
packages never shadow the interpreter's own.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parent.parent / "bin"

# 3 hotspots for one guide (PLCB2_KO):
#   chr12:32679408  real off-target, strong ECS edit (VAF 0.15)  -> WGS should detect
#   chr1:1000000    ECS-negative (VAF 0)                          -> true negative
#   chr7:5000000    ECS low VAF (0.008), below the WGS depth floor -> WGS misses (by design)
HOTSPOTS = [
    ("chr12", 32679408, 32679409, 0.15, 1, "MISMATCHES=1;PAM=AGG"),
    ("chr1",  1000000,  1000001,  0.00, 0, "MISMATCHES=3;PAM=TGG"),
    ("chr7",  5000000,  5000001,  0.008, 0, "MISMATCHES=2;PAM=GGG"),
]
ECS_COLS = ["chrom", "start", "end", "pam_positions", "total_reads", "indel_reads",
            "indel_fraction", "control_reads", "control_indel_reads",
            "control_indel_fraction", "indel_count", "indel_info", "bnd_count",
            "bnd_info", "target_info", "is_target"]


def run(script, *args, cwd):
    """Run a bin/ script as a subprocess; return the CompletedProcess (asserts rc==0)."""
    proc = subprocess.run([sys.executable, str(BIN / script), *map(str, args)],
                          cwd=cwd, capture_output=True, text=True)
    assert proc.returncode == 0, f"{script} failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


def _ecs_table(scale):
    rows = []
    for chrom, start, end, base_vaf, is_tgt, tinfo in HOTSPOTS:
        vaf = round(base_vaf * scale, 4)
        total = 3000
        rows.append([chrom, start, end, str(end), total, int(total * vaf), vaf,
                     2500, 0, 0.0, int(total * vaf), ".", 0, ".", tinfo, is_tgt])
    return pd.DataFrame(rows, columns=ECS_COLS)


@pytest.fixture
def workspace(tmp_path):
    """A tmp dir populated with realistic ECS tables, samplesheet, simulated WGS score
    output and a PoN worklist. Returns the Path."""
    d = tmp_path
    _ecs_table(1.0).to_csv(d / "PLCB2_ecs_1.offtarget_analysis.tsv", sep="\t", index=False)
    _ecs_table(0.8).to_csv(d / "PLCB2_ecs_2.offtarget_analysis.tsv", sep="\t", index=False)

    pd.DataFrame([
        ["PLCB2_ecs_1", "ecs", "PLCB2_KO", "/abs/ed1.cram", "/abs/ctl.cram", "/abs/t.vcf", ""],
        ["PLCB2_ecs_2", "ecs", "PLCB2_KO", "/abs/ed2.cram", "/abs/ctl.cram", "/abs/t.vcf", ""],
        ["PLCB2_wgs_1", "wgs", "PLCB2_KO", "/abs/w1_tumor.cram", "/abs/w1.cram", "", "/abs/w1.hard-filtered.vcf.gz"],
    ], columns=["sample", "datatype", "guide", "edited_cram", "control_cram", "target_file", "vcf"]
    ).to_csv(d / "samplesheet.csv", index=False)

    # simulated SCORE_HOTSPOTS output: chr12 detected, chr7 below HI, chr1 nothing
    pd.DataFrame([
        ["PLCB2_wgs_1", "chr12", 32679408, 1, 0.16, 0.0, 1, 0.85, 0.90, "LIKELY EDIT", 42, 32679408],
        ["PLCB2_wgs_1", "chr1",  1000000,  3, 0.00, 0.0, 0, 0.00, 0.00, "ARTIFACT (shape)", 40, 1000000],
        ["PLCB2_wgs_1", "chr7",  5000000,  2, 0.03, 0.0, 1, 0.30, 0.12, "POSSIBLE — review", 38, 5000000],
    ], columns=["sample", "chrom", "start", "min_mm", "indel_frac", "ctrl_if", "modal_len",
                "conc_ratio", "score", "verdict", "spanning", "modal_pos"]
    ).to_csv(d / "wgs_hotspot_scores.csv", index=False)

    # simulated PoN worklist: chr12 hit 1bp off a known hotspot, chr3 novel
    pd.DataFrame([
        [1, "PLCB2_wgs_1", "chr12", 32679409, "A", 0.16, 0.90, "LIKELY EDIT"],
        [2, "PLCB2_wgs_1", "chr3",  9999999,  "AT", 0.20, 0.88, "LIKELY EDIT"],
    ], columns=["rank", "sample", "chrom", "start", "alt", "dragen_af", "score", "verdict_pon"]
    ).to_csv(d / "worklist_pon.csv", index=False)
    return d
