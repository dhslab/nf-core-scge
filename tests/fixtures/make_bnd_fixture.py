#!/usr/bin/env python
"""Regenerate the committed breakend fixture in tests/fixtures/bnd/.

modules/local/review_bnd_snapshots.nf.test runs the REAL junction renderer in the real
container, so it needs real CRAMs -- but a cohort CRAM is tens of GB and cannot live in
git. This builds the smallest synthetic pair that still exercises the thing the figure
exists to show.

    chr1, 8 kb
    left  breakpoint @ 2000   (1-based)
    right breakpoint @ 6000   -- a 4,000 bp "deletion", the multi-cut case
                                 seen 8 times in the CAR-T cohort

    edited.cram
        20 clean pairs spanning each breakpoint          -> the background
        14 SPLIT reads: primary clipped at one breakpoint, plus a
           supplementary segment at the partner, carrying an SA tag       <- the signal
    control.cram
        20 clean pairs at each breakpoint, no split reads, no SA tags     <- the contrast

The split reads are the whole point. bin/pileup_snapshot.py used to drop
`is_supplementary` unconditionally, which hid exactly the alignments that constitute a
junction; `snapshot(..., keep_supplementary=True, highlight=...)` is what this fixture
proves is wired up. A fixture without SA tags would render two empty pileups and pass.

The queue is real, not synthetic: tests/fixtures/bnd_review_queue_cart.tsv is the actual
25-row cohort queue. This script writes a 5-row queue in the same schema pointing at the
synthetic contig, because the nf-test needs coordinates that exist in these CRAMs. The
row-to-junction collapse is covered against the real queue in tests/test_bnd_snapshots.py.

Run (from the repo root, in an env with pysam):

    python tests/fixtures/make_bnd_fixture.py
"""
import random
from pathlib import Path

import pysam

OUT = Path(__file__).resolve().parent / "bnd"
CHROM = "chr1"
CONTIG_LEN = 8000
LEFT = 2000          # 1-based breakpoint
RIGHT = 6000
READ_LEN = 100
N_CLEAN = 20         # clean pairs per breakpoint
N_SPLIT = 14         # split reads carrying the junction
SAMPLE = "BNDDEMO"

# 1-based -> 0-based
L0, R0 = LEFT - 1, RIGHT - 1


def write_fasta(path):
    random.seed(11)
    seq = "".join(random.choice("ACGT") for _ in range(CONTIG_LEN))
    with open(path, "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, CONTIG_LEN, 60):
            fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def seg(fasta, name, start, cigar, *, supplementary=False, reverse=False, sa=None,
        mapq=60):
    """One alignment record. `cigar` is a list of (op, len); S is filled with the
    reference base so the clip is not a run of a single letter."""
    out, ref = [], start
    for op, ln in cigar:
        if op in (0, 7, 8):
            out.append(fasta.fetch(CHROM, ref, ref + ln)); ref += ln
        elif op in (2, 3):
            ref += ln
        elif op == 1:
            out.append("A" * ln)
        elif op == 4:
            # soft clip: sequence that belongs at the OTHER side of the junction
            out.append(fasta.fetch(CHROM, max(0, ref), max(0, ref) + ln))
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "".join(out)
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigar = cigar
    a.is_paired = False
    a.is_reverse = reverse
    a.is_supplementary = supplementary
    a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
    a.set_tag("NM", 0, value_type="i")
    if sa:
        a.set_tag("SA", sa, value_type="Z")
    return a


def clean_pile(fasta, prefix, centre, n):
    """Reads tiled across a breakpoint with no junction evidence."""
    reads = []
    for i in range(n):
        start = centre - READ_LEN // 2 - 25 + (i % 10) * 5
        reads.append(seg(fasta, f"{prefix}{i}", start, [(0, READ_LEN)],
                         reverse=bool(i % 2)))
    return reads


def split_pair(fasta, name, i):
    """One read spanning the junction: primary clipped at LEFT, supplementary at RIGHT.

    Both segments carry an SA tag naming the other, which is what
    bnd_snapshots.sa_supporters() keys on. The clip lengths sum to the read length, as a
    real split alignment's do.
    """
    keep = 55 + (i % 8)                 # aligned on the left side of the cut
    clip = READ_LEN - keep
    lstart = L0 - keep + 1              # so the alignment ends exactly at the breakpoint
    rstart = R0                         # partner picks up at the far breakpoint

    sa_right = f"{CHROM},{rstart + 1},+,{keep}S{clip}M,60,0;"
    sa_left = f"{CHROM},{lstart + 1},+,{keep}M{clip}S,60,0;"
    return [
        seg(fasta, name, lstart, [(0, keep), (4, clip)], sa=sa_right),
        seg(fasta, name, rstart, [(4, keep), (0, clip)], supplementary=True, sa=sa_left),
    ]


def build(fasta, fasta_path, out_cram, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}]}
    reads.sort(key=lambda r: r.reference_start)
    tmp = str(out_cram) + ".tmp.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    sorted_bam = str(out_cram) + ".sorted.bam"
    pysam.sort("-o", sorted_bam, tmp)
    with pysam.AlignmentFile(sorted_bam) as src, \
            pysam.AlignmentFile(str(out_cram), "wc", template=src,
                                reference_filename=str(fasta_path)) as out:
        for r in src:
            out.write(r)
    pysam.index(str(out_cram))
    Path(tmp).unlink()
    Path(sorted_bam).unlink()


QUEUE_COLS = ["sample_name", "chrom", "pos", "chrom2", "pos2", "strands", "cut_dist",
              "reads", "control_reads_at_event", "site_start", "site_end", "is_target",
              "site_total_reads", "bin", "partner_bin", "n_partners", "span",
              "interchromosomal", "far_end_on_target", "why_dropped", "call"]


def write_queue(path):
    """5 rows describing ONE junction, in the cohort's own reporting pattern:
    both directions, a few bp of jitter, both strand orientations."""
    lbin, rbin = f"{CHROM}:{LEFT // 1000}", f"{CHROM}:{RIGHT // 1000}"
    rows = [
        (LEFT,     RIGHT,     "-+", 3, lbin, rbin),
        (LEFT + 2, RIGHT - 1, "+-", 4, lbin, rbin),
        (RIGHT - 2, LEFT + 1, "-+", 3, rbin, lbin),
        (RIGHT,    LEFT,      "-+", 4, rbin, lbin),
        (RIGHT,    LEFT + 1,  "+-", 4, rbin, lbin),
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(QUEUE_COLS) + "\n")
        for pos, pos2, strands, reads, b, pb in rows:
            fh.write("\t".join(str(x) for x in [
                SAMPLE, CHROM, pos, CHROM, pos2, strands, 2, reads, 0,
                pos - 20, pos + 20, 1, 40, b, pb, 1, abs(pos2 - pos),
                "False", 1, "", "multi-cut deletion"]) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fasta_path = OUT / "ref.fa"
    fasta = write_fasta(fasta_path)

    edited = clean_pile(fasta, "clnL", L0, N_CLEAN) + clean_pile(fasta, "clnR", R0, N_CLEAN)
    for i in range(N_SPLIT):
        edited += split_pair(fasta, f"split{i}", i)

    control = clean_pile(fasta, "ctlL", L0, N_CLEAN) + clean_pile(fasta, "ctlR", R0, N_CLEAN)

    build(fasta, fasta_path, OUT / "edited.cram", edited)
    build(fasta, fasta_path, OUT / "control.cram", control)
    write_queue(OUT / "bnd_review_queue.tsv")

    # The renderer takes CRAM paths from a map rather than staged files -- same reason
    # production does: the queue names arbitrary samples and staging every cohort CRAM
    # would copy TBs. Paths are filled in by the test at run time; this is the template.
    (OUT / "cram_map.tsv").write_text(f"{SAMPLE}\tedited.cram\tcontrol.cram\n")

    print(f"wrote {OUT}: {len(edited)} edited records ({N_SPLIT} split), "
          f"{len(control)} control records")


if __name__ == "__main__":
    main()
