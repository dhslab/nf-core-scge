"""add_normal_counts must count control alt support in FRAGMENTS, not reads.

Depth has always been fragment-level -- `len(total_reads)` where `total_reads` is a set of
query_name -- so R1 and R2 of one physical molecule contribute one unit of depth. Alt support was
counted per read, so an overlapping mate pair that both showed the variant contributed TWO to the
numerator and ONE to the denominator. That biases the background rate upward wherever mates
overlap, which on real CAR-T control CRAMs is ~49% of fragments, and it is one route to the
`control_alt_counts > control_total_counts` rows that make the (k, n) pair invalid for betabinom.

These tests pin the invariant: one fragment, one vote, however many of its reads are visible.

Requires the caller's own imports (pysam, edlib, pyranges) -- run inside
ghcr.io/dhslab/docker-scge-offtarget if the host interpreter lacks them.
"""
import importlib.util
import random
from pathlib import Path

import pandas as pd
import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("edlib")

# find_edited_reads imports pyranges at module scope for target-window merging, which
# add_normal_counts never touches. Stub it if absent so this test runs on any interpreter that
# has pysam and edlib, rather than silently skipping and leaving the fix unverified.
try:                                                    # pragma: no cover
    import pyranges                                     # noqa: F401
except ImportError:                                     # pragma: no cover
    import sys, types
    sys.modules["pyranges"] = types.ModuleType("pyranges")

BIN = Path(__file__).resolve().parents[1] / "bin"
CHROM, CONTIG_LEN = "chr1", 2000
DEL_POS, DEL_LEN = 1001, 4          # 1-based deletion start, length


def _load_caller():
    spec = importlib.util.spec_from_file_location("find_edited_reads",
                                                  BIN / "find_edited_reads.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def workspace(tmp_path):
    """A reference plus a control BAM whose alt-supporting pairs OVERLAP each other."""
    random.seed(11)
    seq = "".join(random.choice("ACGT") for _ in range(CONTIG_LEN))
    fa = tmp_path / "ref.fa"
    with open(fa, "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, CONTIG_LEN, 60):
            fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(fa))
    fasta = pysam.FastaFile(str(fa))

    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}]}
    bam = tmp_path / "control.bam"

    def deleted_seq(start, length):
        """Read sequence carrying the deletion, as ref-minus-the-deleted-bases."""
        left = fasta.fetch(CHROM, start, DEL_POS - 1)
        right = fasta.fetch(CHROM, DEL_POS - 1 + DEL_LEN,
                            start + length + DEL_LEN)
        return left + right

    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        # 3 pairs where BOTH mates span the deletion -- the double-counting case.
        for i in range(3):
            start = DEL_POS - 41
            left_len = DEL_POS - 1 - start
            right_len = 60 - left_len
            cig = [(0, left_len), (2, DEL_LEN), (0, right_len)]
            s = deleted_seq(start, 60)
            for is_r1 in (True, False):
                a = pysam.AlignedSegment()
                a.query_name = f"overlapping_pair_{i}"     # SAME name for both mates
                a.query_sequence = s
                a.flag = 99 if is_r1 else 147
                a.reference_id = 0
                a.reference_start = start
                a.mapping_quality = 60
                a.cigartuples = cig
                a.next_reference_id = 0
                a.next_reference_start = start
                a.template_length = 60 if is_r1 else -60
                a.query_qualities = pysam.qualitystring_to_array("I" * len(s))
                out.write(a)
        # 2 pairs where only R1 spans it -- one read, one fragment, unambiguous.
        for i in range(2):
            start = DEL_POS - 41
            left_len = DEL_POS - 1 - start
            right_len = 60 - left_len
            a = pysam.AlignedSegment()
            a.query_name = f"single_end_pair_{i}"
            a.query_sequence = deleted_seq(start, 60)
            a.flag = 99
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, left_len), (2, DEL_LEN), (0, right_len)]
            a.next_reference_id = 0
            a.next_reference_start = 1500
            a.template_length = 600
            a.query_qualities = pysam.qualitystring_to_array("I" * 60)
            out.write(a)
    pysam.index(str(bam))
    return {"fasta": fasta, "bam": bam}


def _run(workspace):
    mod = _load_caller()
    ref = workspace["fasta"].fetch(CHROM, DEL_POS - 2, DEL_POS - 1 + DEL_LEN)
    alt = ref[0]
    # alttype drives generate_contig's BND-vs-linear branch; a plain deletion is linear.
    df = pd.DataFrame([{"chrom": CHROM, "pos": DEL_POS - 1, "ref": ref, "alt": alt,
                        "alttype": "DEL"}])
    with pysam.AlignmentFile(str(workspace["bam"]), "rb") as bam:
        reads = list(bam.fetch(CHROM, DEL_POS - 400, DEL_POS + 400))
    return mod.add_normal_counts(df, reads, workspace["fasta"])


def test_overlapping_mates_count_once(workspace):
    """5 supporting FRAGMENTS (3 overlapping pairs + 2 single) -- not 8 supporting reads."""
    out = _run(workspace)
    assert int(out["control_alt_counts"].iloc[0]) == 5


def test_alt_never_exceeds_depth(workspace):
    """The (k, n) pair must be valid: betabinom is undefined for k > n."""
    out = _run(workspace)
    assert int(out["control_alt_counts"].iloc[0]) <= int(out["control_total_counts"].iloc[0])


def test_depth_is_fragment_level(workspace):
    """5 fragments were written; both mates of three of them are present."""
    out = _run(workspace)
    assert int(out["control_total_counts"].iloc[0]) == 5
