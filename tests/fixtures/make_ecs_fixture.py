#!/usr/bin/env python
"""Regenerate the committed ECS fixture in tests/fixtures/ecs/.

The nf-test suite runs the REAL ECS caller in the real container, so it needs real
CRAMs — but a cohort CRAM is tens of GB and cannot live in git. This builds a
deliberately tiny synthetic one instead: a 3 kb contig, two targets, 35 read pairs.
The whole fixture is a few KB, so CI can check it out and run the caller on it.

The read construction is NOT duplicated here. It is imported from tests/conftest.py,
which is the single definition of what a synthetic ECS pileup looks like, so the
pytest suite and the nf-test suite can never drift apart on the fixture's meaning.

What the fixture contains (see conftest for the constants):

    chr1, 3000 bp
    target A @ 1001  -- the edited site: 12 read pairs carry a 5 bp deletion
    target B @ 1201  -- a quiet site 200 bp away, close enough that the +/-150 bp
                        fetch windows overlap, so some reads are visited twice

    edited.cram   12 edited + 8 clean + 2 duplicate + 2 low-MAPQ + 2 high-mismatch
                  + 6 in-the-overlap + 3 spanning-B-only  = 35 pairs / 70 records
    control.cram  20 clean pairs at the same site (no edit)

The awkward cases are the point: a duplicate, a MAPQ-3 pair and an NM-6 pair must be
classed as skipped rather than counted, and the pairs in the window overlap must be
written to the tagged BAM exactly once even though two targets both visit them.

Run (from the repo root, inside the off-target container):

    python tests/fixtures/make_ecs_fixture.py

Regenerate only when the read model changes; the outputs are committed.
"""
import importlib
import importlib.util
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUT = HERE / "ecs"


class _PytestShim(types.ModuleType):
    """Just enough `pytest` to import conftest.py outside a test run.

    The off-target container ships pandas and pysam but not pytest, and this generator
    has to run there because that is where the pinned pysam lives. conftest only touches
    pytest at import time for the @pytest.fixture decorators, so an identity decorator is
    a faithful stand-in. We only ever call the plain builder functions, never a fixture.
    """

    def __init__(self):
        super().__init__("pytest")

    @staticmethod
    def fixture(func=None, **_kw):
        return func if func is not None else (lambda f: f)

    @staticmethod
    def importorskip(name, *_a, **_kw):
        return importlib.import_module(name)

    class approx:                                     # noqa: N801 - mirrors pytest's name
        def __init__(self, *_a, **_kw):
            pass


def _load_conftest():
    """Import tests/conftest.py as a module without needing pytest to collect it."""
    sys.modules.setdefault("pytest", _PytestShim())
    spec = importlib.util.spec_from_file_location(
        "ecs_conftest", REPO / "tests" / "conftest.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    import pysam

    c = _load_conftest()
    OUT.mkdir(parents=True, exist_ok=True)

    fasta_path = OUT / "ref.fa"
    fasta = c._write_ecs_fasta(fasta_path)

    # targets VCF — the --target-file the caller takes
    header = pysam.VariantHeader()
    header.contigs.add(c.ECS_CHROM, length=c.ECS_CONTIG_LEN)
    targets = OUT / "targets.vcf"
    with pysam.VariantFile(str(targets), "w", header=header) as vout:
        for pos in (c.ECS_TARGET_A, c.ECS_TARGET_B):
            rec = vout.new_record()
            rec.chrom = c.ECS_CHROM
            rec.pos = pos
            rec.id = "."
            rec.ref = fasta.fetch(c.ECS_CHROM, pos - 1, pos)
            rec.alts = ("N",)
            rec.filter.add("PASS")
            vout.write(rec)

    del_cigar = [(0, 60), (2, c.ECS_DEL_LEN), (0, 40)]
    edited = []
    for i in range(c.ECS_N_EDIT):
        edited += c._ecs_pair(fasta, f"edit{i}", 940, del_cigar, 1180, r1_nm=c.ECS_DEL_LEN)
    for i in range(c.ECS_N_WT):
        edited += c._ecs_pair(fasta, f"wt{i}", 940, [(0, 100)], 1180)
    for i in range(c.ECS_N_DUP):
        edited += c._ecs_pair(fasta, f"dup{i}", 940, [(0, 100)], 1180, dup=True)
    for i in range(c.ECS_N_LOWMAPQ):
        edited += c._ecs_pair(fasta, f"lowq{i}", 940, [(0, 100)], 1180, mapq=3)
    for i in range(c.ECS_N_MISMATCH):
        edited += c._ecs_pair(fasta, f"mm{i}", 940, [(0, 100)], 1180, r1_nm=6, r2_nm=6)
    for i in range(c.ECS_N_OVERLAP):
        edited += c._ecs_pair(fasta, f"both{i}", 1080, [(0, 100)], 1220)
    for i in range(c.ECS_N_SPAN_B):
        edited += c._ecs_pair(fasta, f"spanb{i}", 1150, [(0, 100)], 1220)

    control = []
    for i in range(c.ECS_N_EDIT + c.ECS_N_WT):
        control += c._ecs_pair(fasta, f"ctl{i}", 940, [(0, 100)], 1180)

    c._build_ecs_cram(fasta, fasta_path, OUT / "edited.bam", edited)
    c._build_ecs_cram(fasta, fasta_path, OUT / "control.bam", control)

    # the intermediate BAMs are a build artifact of the CRAM writer; only ship CRAMs
    for stem in ("edited", "control"):
        for ext in (".bam", ".bam.bai"):
            p = OUT / f"{stem}{ext}"
            if p.exists():
                p.unlink()

    print(f"wrote fixture to {OUT}")
    for p in sorted(OUT.iterdir()):
        print(f"  {p.name:20s} {p.stat().st_size:>8,d} bytes")
    print(f"\ntotal {sum(p.stat().st_size for p in OUT.iterdir()):,d} bytes")
    print(f"edited records: {len(edited)}  control records: {len(control)}")


if __name__ == "__main__":
    sys.exit(main())
