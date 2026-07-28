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


# ---------------------------------------------------------------------------
# Synthetic aligned-read workspace for find_edited_reads.py
#
# The glue fixtures above are pure tables; the read-level tagging path needs
# real alignments. This builds the smallest workspace that reproduces the two
# things that make tagging non-trivial:
#
#   * two targets 200 bp apart -- far enough that pyranges clusters them
#     separately (slack = --target-window = 150), close enough that their
#     +/-150 fetch windows overlap, so some reads are visited twice;
#   * a read pair sitting in that overlap which one target calls Unedited_WT
#     and the other cannot evaluate, i.e. the tag-conflict case.
# ---------------------------------------------------------------------------

ECS_CHROM = "chr1"
ECS_CONTIG_LEN = 3000
ECS_TARGET_A = 1001        # 1-based; the edited site
ECS_TARGET_B = 1201        # 1-based; a quiet site
ECS_DEL_LEN = 5
ECS_N_EDIT = 12            # read pairs carrying the deletion
ECS_N_WT = 8               # read pairs spanning target A cleanly
ECS_N_DUP = 2              # duplicate-flagged pairs
ECS_N_LOWMAPQ = 2          # MAPQ below the default floor of 20
ECS_N_MISMATCH = 2         # NM above the default ceiling of 4
ECS_N_OVERLAP = 6          # pairs inside the window overlap, spanning neither target
ECS_N_SPAN_B = 3           # pairs in the overlap that DO span target B (conflict case)


def _write_ecs_fasta(path):
    import random
    random.seed(7)
    seq = "".join(random.choice("ACGT") for _ in range(ECS_CONTIG_LEN))
    with open(path, "w") as fh:
        fh.write(f">{ECS_CHROM}\n")
        for i in range(0, ECS_CONTIG_LEN, 60):
            fh.write(seq[i:i + 60] + "\n")
    import pysam
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def _ecs_pair(fasta, name, r1_start, r1_cigar, r2_start, mapq=60, dup=False,
              r1_nm=0, r2_nm=0):
    """One properly-paired FR pair (R1 forward, R2 reverse) with real sequence."""
    import pysam

    def qseq(start, cigar):
        out, ref = [], start
        for op, ln in cigar:
            if op == 0:                       # M
                out.append(fasta.fetch(ECS_CHROM, ref, ref + ln)); ref += ln
            elif op == 2:                     # D
                ref += ln
            elif op == 1:                     # I
                out.append("A" * ln)
            elif op == 4:                     # S
                out.append("T" * ln)
        return "".join(out)

    span = r2_start + 100 - r1_start
    reads = []
    for is_read1, start, cigar, mate_start, nm in (
            (True, r1_start, r1_cigar, r2_start, r1_nm),
            (False, r2_start, [(0, 100)], r1_start, r2_nm)):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = qseq(start, cigar)
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = mapq
        a.cigar = cigar
        a.next_reference_id = 0
        a.next_reference_start = mate_start
        a.template_length = span if is_read1 else -span
        a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
        a.is_paired = True
        a.is_proper_pair = True
        a.is_read1 = is_read1
        a.is_read2 = not is_read1
        a.is_reverse = not is_read1
        a.mate_is_reverse = is_read1
        a.is_duplicate = dup
        a.set_tag("NM", nm, value_type="i")
        a.set_tag("MC", "100M", value_type="Z")
        reads.append(a)
    return reads


def _build_ecs_cram(fasta, fasta_path, out_bam, reads):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": ECS_CHROM, "LN": ECS_CONTIG_LEN}]}
    reads.sort(key=lambda r: r.reference_start)
    tmp = str(out_bam) + ".tmp.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", str(out_bam), tmp)
    Path(tmp).unlink()
    pysam.index(str(out_bam))
    # production opens its inputs as CRAM, so hand the caller a CRAM
    cram = str(out_bam)[:-4] + ".cram"
    with pysam.AlignmentFile(str(out_bam)) as src, \
            pysam.AlignmentFile(cram, "wc", template=src,
                                reference_filename=str(fasta_path)) as out:
        for r in src:
            out.write(r)
    pysam.index(cram)
    return cram


@pytest.fixture
def ecs_reads_workspace(tmp_path):
    """Synthetic reference + targets VCF + edited/control CRAMs for the ECS caller.

    Returns a dict of paths plus the read counts the tags should reconcile with.
    """
    pysam = pytest.importorskip("pysam")
    d = tmp_path
    fasta_path = d / "ref.fa"
    fasta = _write_ecs_fasta(fasta_path)

    header = pysam.VariantHeader()
    header.contigs.add(ECS_CHROM, length=ECS_CONTIG_LEN)
    targets = d / "targets.vcf"
    with pysam.VariantFile(str(targets), "w", header=header) as vout:
        for pos in (ECS_TARGET_A, ECS_TARGET_B):
            rec = vout.new_record()
            rec.chrom = ECS_CHROM
            rec.pos = pos
            rec.id = "."
            rec.ref = fasta.fetch(ECS_CHROM, pos - 1, pos)
            rec.alts = ("N",)
            rec.filter.add("PASS")
            vout.write(rec)

    del_cigar = [(0, 60), (2, ECS_DEL_LEN), (0, 40)]
    edited = []
    for i in range(ECS_N_EDIT):
        edited += _ecs_pair(fasta, f"edit{i}", 940, del_cigar, 1180, r1_nm=ECS_DEL_LEN)
    for i in range(ECS_N_WT):
        edited += _ecs_pair(fasta, f"wt{i}", 940, [(0, 100)], 1180)
    for i in range(ECS_N_DUP):
        edited += _ecs_pair(fasta, f"dup{i}", 940, [(0, 100)], 1180, dup=True)
    for i in range(ECS_N_LOWMAPQ):
        edited += _ecs_pair(fasta, f"lowq{i}", 940, [(0, 100)], 1180, mapq=3)
    for i in range(ECS_N_MISMATCH):
        edited += _ecs_pair(fasta, f"mm{i}", 940, [(0, 100)], 1180, r1_nm=6, r2_nm=6)
    for i in range(ECS_N_OVERLAP):
        edited += _ecs_pair(fasta, f"both{i}", 1080, [(0, 100)], 1220)
    # R1 lands in both padded windows: target A cannot evaluate it, target B
    # calls it reference. Unedited_WT must win, and it must be written once.
    for i in range(ECS_N_SPAN_B):
        edited += _ecs_pair(fasta, f"spanb{i}", 1150, [(0, 100)], 1220)

    control = []
    for i in range(ECS_N_EDIT + ECS_N_WT):
        control += _ecs_pair(fasta, f"ctl{i}", 940, [(0, 100)], 1180)

    return {
        "dir": d,
        "fasta": fasta_path,
        "targets": targets,
        "edited": _build_ecs_cram(fasta, fasta_path, d / "edited.bam", edited),
        "control": _build_ecs_cram(fasta, fasta_path, d / "control.bam", control),
        "n_edit_reads": ECS_N_EDIT,
        "del_len": ECS_DEL_LEN,
        "total_records": len(edited),
    }
