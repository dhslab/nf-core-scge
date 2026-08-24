"""bin/bnd_snapshots.py -- the row-to-junction collapse, which is the whole design.

The breakend queue reports one physical event several times: once from each end, at a
few bp of position jitter, and under both strand orientations. On the 32-sample CAR-T
cohort that is 25 rows describing 8 junctions. Rendering per row would produce 25
figures, three of which are the same inversion drawn from opposite directions, and the
per-row read counts (3, 4, 3, 4, 4) would understate an event carrying 18 reads.

So `collapse()` is not a convenience -- it is the difference between a figure that says
what happened and a stack of figures that say it five times, quieter each time. This
file pins that behaviour against the REAL queue, checked in at
tests/fixtures/bnd_review_queue_cart.tsv (3.5 KB), rather than a synthetic one: the
jitter pattern and the both-directions reporting are precisely the properties a
hand-written fixture would smooth away.

No CRAM, no pysam, no plotting -- this is the arithmetic only. The rendering path is
covered by modules/local/review_bnd_snapshots.nf.test, which runs the real script.
"""
import importlib.util
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "bnd_review_queue_cart.tsv"
BIN = Path(__file__).resolve().parent.parent / "bin"


def _load():
    """Import bin/bnd_snapshots.py by path -- bin/ is not a package."""
    spec = importlib.util.spec_from_file_location("bnd_snapshots", BIN / "bnd_snapshots.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bs = _load()


@pytest.fixture(scope="module")
def queue():
    return pd.read_csv(FIXTURE, sep="\t", dtype=str)


@pytest.fixture(scope="module")
def junctions(queue):
    return bs.collapse(queue)


def by_sample(junctions):
    return {j["sample"]: j for j in junctions}


# --------------------------------------------------------------------------------
# The collapse itself
# --------------------------------------------------------------------------------

def test_the_real_queue_is_25_rows(queue):
    """Guards the fixture: if this changes, every count below is measuring something else."""
    assert len(queue) == 25
    assert queue["sample_name"].nunique() == 8


def test_25_rows_collapse_to_8_junctions(junctions):
    assert len(junctions) == 8
    assert {j["sample"] for j in junctions} == {
        "ARID4A-KO-DNA", "BRAF-KO-CART-DNA", "PDCD4-KO-CART", "KLF12-KO-DNA",
        "PDE7A-KO-DNA", "IKZF2-KO-DNA", "RXRB-KO-DNA", "ZEB2-KO-DNA",
    }


def test_rows_reported_from_opposite_ends_land_in_one_junction(junctions):
    """ARID4A is the clean case: 2 rows call it chr14:58301->58330, 3 call it the reverse.

    A key that did not sort its two bins would produce two junctions here, and the figure
    would draw the same inversion twice, mirrored.
    """
    j = by_sample(junctions)["ARID4A-KO-DNA"]
    assert j["n_rows"] == 5
    assert j["left"][0] == j["right"][0] == "chr14"
    assert j["left"][1] < j["right"][1]


def test_every_row_is_accounted_for(queue, junctions):
    assert sum(j["n_rows"] for j in junctions) == len(queue)
    assert sum(j["reads"] for j in junctions) == queue["reads"].astype(int).sum() == 150


# --------------------------------------------------------------------------------
# Aggregation -- the numbers the figure prints
# --------------------------------------------------------------------------------

def test_aggregate_support_exceeds_any_single_row(queue, junctions):
    """The reason per-row counts must not be quoted: ARID4A's largest row says 4 reads.

    18 is the sum over the rows describing that junction. It is support counted from both
    ends of the same event, so it is an aggregate of the caller's row-level evidence --
    not an independent read count, and deliberately more conservative than the SA-tag
    evidence the figure highlights (48 of 60 SA reads at this locus point at the partner).
    """
    j = by_sample(junctions)["ARID4A-KO-DNA"]
    rows = queue[queue["sample_name"] == "ARID4A-KO-DNA"]["reads"].astype(int)

    assert j["reads"] == 18 == rows.sum()
    assert j["reads"] > rows.max() * 4


@pytest.mark.parametrize("sample,left,right,span,reads", [
    # positions are the highest-support end, not the first seen
    ("ARID4A-KO-DNA", 58301619, 58330082, 28463, 18),
    ("PDCD4-KO-CART", 110887716, 110890582, 2866, 35),
    ("IKZF2-KO-DNA", 213022046, 213147793, 125747, 12),
    ("BRAF-KO-CART-DNA", 140801460, 140834645, 33185, 10),
])
def test_representative_ends_and_span(junctions, sample, left, right, span, reads):
    j = by_sample(junctions)[sample]
    assert j["left"][1] == left
    assert j["right"][1] == right
    assert j["span"] == span
    assert j["reads"] == reads


def test_depth_and_control_are_carried(junctions):
    """site_total_reads is the denominator a reviewer needs; control must stay visible."""
    j = by_sample(junctions)["ARID4A-KO-DNA"]
    assert j["depth"] == 129
    assert j["control"] == 0
    assert all(x["control"] == 0 for x in junctions), "cohort has no control support"


def test_jittered_positions_are_kept_as_marks(junctions):
    """Every distinct end position survives, so the figure can show the jitter."""
    j = by_sample(junctions)["ARID4A-KO-DNA"]
    assert sorted(j["left"][2]) == [58301619, 58301620, 58301622, 58301623]
    assert j["left"][1] in j["left"][2]
    assert j["right"][1] in j["right"][2]


# --------------------------------------------------------------------------------
# Canonical ordering -- so a figure is reproducible
# --------------------------------------------------------------------------------

def test_left_is_always_before_right(junctions):
    for j in junctions:
        assert (j["left"][0], j["left"][1]) < (j["right"][0], j["right"][1])


def test_ordering_is_independent_of_row_order(queue):
    """Reversing the input must not mirror any figure."""
    forward = {j["sample"]: (j["left"], j["right"]) for j in bs.collapse(queue)}
    reverse = {j["sample"]: (j["left"], j["right"])
               for j in bs.collapse(queue.iloc[::-1].reset_index(drop=True))}
    assert forward == reverse


def test_metadata_passthrough(junctions):
    j = by_sample(junctions)["BRAF-KO-CART-DNA"]
    assert j["call"] == "inversion at cut site"
    assert j["interchrom"] is False
    assert set(j["strands"].split(",")) == {"-+", "+-"}

    assert by_sample(junctions)["ARID4A-KO-DNA"]["call"] == "multi-cut inversion"


def test_opposite_strands_are_flagged_inverted(junctions):
    """collapse() must carry orientation, because the schematic branches on it.

    +-/-+ means the two joined segments run in opposite directions: the segment between the
    cuts was flipped, not removed. Drawing the excision cartoon for one of those is worse
    than drawing nothing, since a reader trusts the picture over the caption -- which is
    exactly what shipped until 2026-08-20.
    """
    for j in junctions:
        assert j["inverted"] is True, f"{j['sample']} lost its orientation"


# --------------------------------------------------------------------------------
# Helpers, directly
# --------------------------------------------------------------------------------

def test_junction_bins_is_order_insensitive():
    a = bs.junction_bins({"bin": "chr14:58301", "partner_bin": "chr14:58330"})
    b = bs.junction_bins({"bin": "chr14:58330", "partner_bin": "chr14:58301"})
    assert a == b == ("chr14:58301", "chr14:58330")


def test_end_of_rejects_an_unknown_bin(queue):
    g = queue[queue["sample_name"] == "ARID4A-KO-DNA"]
    with pytest.raises(ValueError):
        bs.end_of(g, "chr1:999999")


def test_empty_queue_collapses_to_nothing(queue):
    assert bs.collapse(queue.iloc[0:0]) == []
