"""bin/bnd_from_indels_to_vcf.py -- the breakend VCF behind the report's SV section.

Why this file exists, in one sentence: for every run ever done before this suite, that
script emitted a VCF containing zero records and nobody noticed.

The cause was `'\\t'` and `'\\n'` written into the source with doubled backslashes -- so
the header was split on the two-character string backslash-t, never matched, and every
row took the `continue` branch. The output was one line of literal backslash-n text,
`wc -l` = 0. It is downstream of nothing that checks it: compile_report_json.nf feeds it
to the per-sample HTML report as --offtarget_svs, and a blank SV section looks exactly
like a sample with no SVs.

So the first test here is deliberately dumb -- it asserts the file contains real tab and
newline bytes -- because that is the assertion whose absence cost 32 samples' worth of
reports. The rest cover what the rewrite added: mate linking and the SR/CTRL evidence
fields that the old `parts[:5]` slice threw away.

Fixtures are REAL bnd_info payloads lifted from
results_cart_nopon/ARID4A-KO-DNA/ARID4A-KO-DNA.offtarget_analysis.tsv, including the one
genuinely reciprocal pair in that sample (the 28 kb chr14 deletion, reported from both
ends). Synthetic strings would not have caught the bracket-orientation question, which
is why they are copied rather than invented.
"""
import csv
from pathlib import Path

import pytest

from conftest import run

BIN = Path(__file__).resolve().parent.parent / "bin"

# The 19 columns find_edited_reads.py writes. Only bnd_info is read by the script under
# test, but the header must be complete or the DictReader lookup is not being exercised
# the way production exercises it.
COLS = ["chrom", "start", "end", "pam_positions", "total_reads", "indel_reads",
        "indel_fraction", "control_reads", "control_indel_reads",
        "control_indel_fraction", "indel_count", "indel_info", "bnd_count",
        "bnd_info", "control_bnd_reads", "n_control_filtered", "min_cut_distance",
        "target_info", "is_target"]

# chrom|pos|chrom2|pos2|strands|ref|alt|distance|distance2|counts|control_alt_counts
#   0    1     2     3     4     5   6      7         8        9         10
MATE_A = "chr14|58301619|chr14|58330081|-+|G|G]chr14:58330081]|1|1|3|0"
MATE_B = "chr14|58330081|chr14|58301619|+-|A|A]chr14:58301619]|1|1|1|0"
# no reciprocal partner in the file -- the common case
LONE   = "chr2|71504476|chr2|32916560|-+|C|CCAC]chr2:32916560]|8|118233|1|0"
# a second lone junction whose ALT uses the opposite bracket form
LONE2  = "chr3|75672394|chr2|32916405|++|C|]chr2:32916405]GCC|2|118078|7|2"


def write_tsv(path, bnd_infos):
    """One analysis row per bnd_info payload; every other column is filler."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        for i, bi in enumerate(bnd_infos):
            row = ["chr14", 1000 + i, 1001 + i, str(1001 + i), 100, 0, 0.0,
                   100, 0, 0.0, 0, ".", bi.count(";") + 1, bi, 0, 0, 1, ".", 0]
            w.writerow(row)


def convert(tmp_path, bnd_infos, sample="TEST"):
    """Run the script the way modules/local/*.nf runs it; return the VCF text."""
    tsv = tmp_path / "in.offtarget_analysis.tsv"
    out = tmp_path / "out.vcf"
    write_tsv(tsv, bnd_infos)
    run("bnd_from_indels_to_vcf.py", "--meta_id", sample,
        "--indels_path", tsv, "--outfile", out, cwd=tmp_path)
    return out.read_text()


def records(text):
    """The data lines, as dicts keyed by VCF column."""
    keys = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    return [dict(zip(keys, ln.split("\t")))
            for ln in text.splitlines() if ln and not ln.startswith("#")]


def info_of(rec):
    return dict(kv.split("=", 1) for kv in rec["INFO"].split(";") if "=" in kv)


# --------------------------------------------------------------------------------
# The regression that was missing
# --------------------------------------------------------------------------------

def test_output_uses_real_tabs_and_newlines(tmp_path):
    """The exact bug: escapes written doubled, so nothing was ever a delimiter."""
    text = convert(tmp_path, [LONE])

    assert "\\t" not in text, "literal backslash-t in output -- the escape bug is back"
    assert "\\n" not in text, "literal backslash-n in output -- the escape bug is back"
    assert text.count("\n") > 1, "single-line output -- newlines are not newlines"

    header = [ln for ln in text.splitlines() if ln.startswith("#CHROM")]
    assert len(header) == 1
    assert header[0].split("\t") == ["#CHROM", "POS", "ID", "REF", "ALT",
                                     "QUAL", "FILTER", "INFO"]


def test_emits_records_at_all(tmp_path):
    """Before the fix this was 0 for every sample in every run."""
    assert len(records(convert(tmp_path, [LONE, LONE2]))) == 2


# --------------------------------------------------------------------------------
# Content
# --------------------------------------------------------------------------------

def test_ref_and_alt_come_from_the_caller_verbatim(tmp_path):
    """Fields 5 and 6 are already a valid VCF REF/ALT -- they must not be re-derived.

    Both bracket forms appear here on purpose: `t]p]` (a piece to the right joins after
    this base) and `]p]t` (a piece to the left joins before it). Re-deriving orientation
    from `strands` is what the old code attempted, and it only handled '++'.
    """
    recs = {r["CHROM"]: r for r in records(convert(tmp_path, [LONE, LONE2]))}

    assert recs["chr2"]["REF"] == "C"
    assert recs["chr2"]["ALT"] == "CCAC]chr2:32916560]"
    assert recs["chr3"]["REF"] == "C"
    assert recs["chr3"]["ALT"] == "]chr2:32916405]GCC"


def test_support_counts_survive_into_info(tmp_path):
    """SR/CTRL are bnd_info fields 9 and 10, which the old `parts[:5]` slice discarded.

    Without them the report's SV table can only show coordinates, so a 1-read junction
    and a 7-read junction look identical to a reviewer.
    """
    recs = {r["CHROM"]: r for r in records(convert(tmp_path, [LONE, LONE2]))}

    assert info_of(recs["chr2"])["SR"] == "1"
    assert info_of(recs["chr2"])["CTRL"] == "0"
    assert info_of(recs["chr3"])["SR"] == "7"      # field 9 of LONE2
    assert info_of(recs["chr3"])["CTRL"] == "2"    # field 10 of LONE2

    text = convert(tmp_path, [LONE])
    assert "##INFO=<ID=SR," in text
    assert "##INFO=<ID=CTRL," in text


def test_every_record_declares_svtype_bnd(tmp_path):
    for r in records(convert(tmp_path, [MATE_A, MATE_B, LONE])):
        assert info_of(r)["SVTYPE"] == "BND"


# --------------------------------------------------------------------------------
# Mate linking
# --------------------------------------------------------------------------------

def test_reciprocal_mates_resolve_to_each_other(tmp_path):
    """The chr14 pair is one 28 kb deletion reported from both ends.

    MATEID must point at a record that actually exists in this file and whose own MATEID
    points back -- a dangling MATEID is worse than none, because a VCF reader will follow
    it.
    """
    recs = records(convert(tmp_path, [MATE_A, MATE_B]))
    by_id = {r["ID"]: r for r in recs}
    assert len(recs) == 2

    for r in recs:
        info = info_of(r)
        assert "MATEID" in info, f"{r['ID']} lost its mate"
        assert info["MATEID"] in by_id, f"{r['ID']} points at a nonexistent mate"
        assert info_of(by_id[info["MATEID"]])["MATEID"] == r["ID"]

    # both halves of one event carry one EVENT id
    assert len({info_of(r)["EVENT"] for r in recs}) == 1


def test_unmated_junction_has_no_mateid(tmp_path):
    """Most junctions are seen from one end only. They must not invent a partner."""
    recs = records(convert(tmp_path, [LONE]))
    assert len(recs) == 1
    assert "MATEID" not in info_of(recs[0])


def test_duplicate_events_collapse(tmp_path):
    """bnd_info repeats a junction across rows; the VCF must not repeat it.

    Both the within-row and across-row cases, since dedup happens over the whole file.
    """
    text = convert(tmp_path, [";".join([LONE, LONE]), LONE])
    assert len(records(text)) == 1


def test_ids_are_unique_and_output_is_coordinate_sorted(tmp_path):
    recs = records(convert(tmp_path, [MATE_A, MATE_B, LONE, LONE2]))
    ids = [r["ID"] for r in recs]
    assert len(ids) == len(set(ids))

    coords = [(r["CHROM"], int(r["POS"])) for r in recs]
    assert coords == sorted(coords)


# --------------------------------------------------------------------------------
# Degenerate input
# --------------------------------------------------------------------------------

@pytest.mark.parametrize("payload", [".", "", "chr1|100|chr2", "not|enough|fields|at|all"])
def test_malformed_payloads_produce_a_valid_empty_vcf(tmp_path, payload):
    """A short or absent payload is skipped, not crashed on, and the header still lands.

    An empty-but-well-formed VCF is what a sample with no breakends should produce; the
    point is that it is reached by skipping bad rows rather than by failing to parse the
    good ones.
    """
    text = convert(tmp_path, [payload])
    assert text.startswith("##fileformat=VCFv4.2\n")
    assert "#CHROM\tPOS\tID" in text
    assert records(text) == []


def test_sample_id_is_recorded(tmp_path):
    assert "##sampleId=ARID4A-KO-DNA" in convert(tmp_path, [LONE], sample="ARID4A-KO-DNA")
