"""Logic tests for the Unified CRISPR Off-Target Workflow glue scripts.

Scenario (see conftest): one guide, two ECS replicates + one WGS sample, three hotspots
(a real edit, a true negative, and a below-floor low-VAF edit). These lock in the
ECS⋈WGS join and the recall-vs-VAF behaviour that the Nextflow layer depends on.
"""
import subprocess
import sys

import pandas as pd
import pytest
from conftest import BIN, run


def test_hotspot_to_table(workspace):
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", "--edit-threshold", "0.0", cwd=workspace)

    truth = pd.read_csv(workspace / "ecs_hotspot_truth.csv").set_index("chrom")
    # ECS truth takes the strongest evidence across replicates (max VAF)
    assert truth.loc["chr12", "ecs_if"] == 0.15
    assert truth.loc["chr12", "ecs_is_edit"] == 1
    assert truth.loc["chr1", "ecs_is_edit"] == 0
    assert truth.loc["chr7", "ecs_is_edit"] == 1          # 0.008 > threshold 0.0
    assert int(truth["ecs_is_edit"].sum()) == 2

    table = pd.read_csv(workspace / "wgs_hotspot_input_table.csv")
    assert len(table) == 3                                 # 1 WGS sample x 3 hotspots
    assert (table["sample_name"] == "PLCB2_wgs_1").all()
    assert (table["indel_fraction"] == 1.0).all()          # gate-passer, not real signal


def test_join_training_table(workspace):
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    run("join_training_table.py", "--wgs-scores", "wgs_hotspot_scores.csv",
        "--truth", "ecs_hotspot_truth.csv", "--samplesheet", "samplesheet.csv", cwd=workspace)

    tr = pd.read_csv(workspace / "training.tsv", sep="\t").set_index("chrom")
    assert len(tr) == 3
    assert set(tr["label"]) == {0, 1}
    # WGS features and ECS truth land on the same row, keyed on (guide, chrom, start)
    assert tr.loc["chr12", "label"] == 1
    assert tr.loc["chr12", "ecs_if"] == 0.15
    assert tr.loc["chr12", "score"] == 0.90
    assert tr.loc["chr7", "label"] == 1 and tr.loc["chr7", "score"] == 0.12


def test_recall_vs_vaf(workspace):
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    run("join_training_table.py", "--wgs-scores", "wgs_hotspot_scores.csv",
        "--truth", "ecs_hotspot_truth.csv", "--samplesheet", "samplesheet.csv", cwd=workspace)
    proc = run("recall_vs_vaf.py", "--training", "training.tsv", "--hi", "0.60",
               "--target-recall", "0.80", cwd=workspace)

    # overall: 1 of 2 credible ECS edits recovered (chr12 yes @0.15, chr7 no @0.008).
    # Both clear the denominator gates (3000x depth -> 450 and 24 indel reads).
    assert "overall WGS recall of credible ECS edits: 0.50" in proc.stdout
    assert "2 EVALUABLE sites" in proc.stdout
    m = pd.read_csv(workspace / "recall_vs_vaf.csv")
    recalls = dict(zip(m["vaf_bin"], m["recall"]))
    assert recalls["[0.005, 0.01)"] == 0.0                 # below the depth floor
    assert recalls["[0.1, 0.2)"] == 1.0                    # above it
    assert (workspace / "recall_vs_vaf.png").exists()

    # the denominator must be stamped into the file so it can't be read out of context
    for col in ("n_unevaluable", "n_evaluable", "recall_incl_unevaluable",
                "denom_min_ecs_vaf", "denom_min_ecs_reads", "denom_excluded_as_noise"):
        assert col in m.columns, f"{col} missing from recall_vs_vaf.csv"
    assert (m["denom_min_ecs_vaf"] == 0.005).all()
    assert (m["denom_min_ecs_reads"] == 5).all()
    # chr1 (VAF 0, no read support) must be excluded as noise, never counted as a miss
    assert (m["denom_excluded_as_noise"] >= 0).all()
    assert "[0.0, 0.005)" not in recalls


def test_recall_denominator_excludes_ecs_noise(workspace):
    """A high-VAF ECS call with almost no read support is assay noise, not a missed edit:
    raising --min-ecs-reads above its support must DROP it from the denominator rather
    than scoring it as a recall failure."""
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    run("join_training_table.py", "--wgs-scores", "wgs_hotspot_scores.csv",
        "--truth", "ecs_hotspot_truth.csv", "--samplesheet", "samplesheet.csv", cwd=workspace)

    # chr7 has 24 indel reads; require 100 and it must leave the denominator entirely
    proc = run("recall_vs_vaf.py", "--training", "training.tsv", "--min-ecs-reads", "100",
               cwd=workspace)
    m = pd.read_csv(workspace / "recall_vs_vaf.csv")
    assert "[0.005, 0.01)" not in set(m["vaf_bin"])        # chr7 gone, not scored as a miss
    assert "overall WGS recall of credible ECS edits: 1.00" in proc.stdout


def test_reconcile_report(workspace):
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    run("reconcile_offtarget_report.py", "--worklist", "worklist_pon.csv",
        "--truth", "ecs_hotspot_truth.csv", "--pad", "25", "--verdict-col", "verdict_pon",
        cwd=workspace)

    rep = pd.read_csv(workspace / "offtarget_report.csv").set_index("chrom")
    # chr12 candidate is 1bp off the known hotspot -> matched within pad, ECS-confirmed
    assert rep.loc["chr12", "is_hotspot"] == 1
    assert rep.loc["chr12", "ecs_confirmed"] == 1
    assert rep.loc["chr12", "ecs_if"] == 0.15
    # chr3 has no predicted hotspot nearby -> novel candidate
    assert rep.loc["chr3", "is_hotspot"] == 0


def test_join_empty_on_coord_mismatch(workspace):
    """Guardrail: when both arms have rows but coords drift, the join is empty and the
    script must FAIL — a silent empty training.tsv would otherwise pass as a clean run."""
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    w = pd.read_csv(workspace / "wgs_hotspot_scores.csv")
    w["start"] += 5                                        # simulate a 5bp drift
    w.to_csv(workspace / "wgs_shifted.csv", index=False)

    # invoke directly (not the rc==0 run() helper) because we expect a non-zero exit
    proc = subprocess.run(
        [sys.executable, str(BIN / "join_training_table.py"),
         "--wgs-scores", "wgs_shifted.csv", "--truth", "ecs_hotspot_truth.csv",
         "--samplesheet", "samplesheet.csv"],
        cwd=workspace, capture_output=True, text=True)

    assert proc.returncode != 0                            # both arms had rows -> hard fail
    assert "ERROR" in proc.stderr and "empty join" in proc.stderr


# ── verdict() call-arity contract ──────────────────────────────────────────────
# score.verdict() returns a 3-tuple (verdict, score, call_basis). worklist_from_vcf.py
# imports it as S.verdict and unpacks it too, so changing the arity in one file silently
# breaks the other — it did, and only surfaced mid-run on the cluster because the failing
# branch needs a real CRAM. This walks the AST of every caller so a future arity change
# fails here instead of an hour into a cohort run.
def test_verdict_callers_unpack_three_values():
    import ast

    offenders = []
    for path in sorted(BIN.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            call = node.value
            if not isinstance(call, ast.Call):
                continue
            fn = call.func
            name = (fn.attr if isinstance(fn, ast.Attribute)
                    else fn.id if isinstance(fn, ast.Name) else None)
            if name != "verdict":
                continue
            for tgt in node.targets:
                if not isinstance(tgt, ast.Tuple):
                    offenders.append(f"{path.name}:{node.lineno} not a tuple unpack")
                elif len(tgt.elts) != 3:
                    offenders.append(f"{path.name}:{node.lineno} unpacks {len(tgt.elts)}, want 3")

    assert not offenders, "verdict() arity mismatch: " + "; ".join(offenders)


# --------------------------------------------------------------------------
# find_edited_reads.py --tagged-bam-out (read-level tags for IGV review)
#
# The caller visits the same alignment record once per overlapping target, so
# the risk here is a BAM with duplicate records, which IGV silently renders as
# doubled depth. These run the real caller on the synthetic CRAM workspace.
# --------------------------------------------------------------------------
import collections

import conftest as C


def _load_find_edited_reads():
    """Import bin/find_edited_reads.py as a module (it is a script, not a package)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("find_edited_reads",
                                                  BIN / "find_edited_reads.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_caller(ws, *extra):
    run("find_edited_reads.py",
        "--fasta", ws["fasta"], "--edited-bam", ws["edited"],
        "--control-bam", ws["control"], "--target-file", ws["targets"],
        "-o", ws["dir"] / "out.tsv", *extra, cwd=ws["dir"])
    return pd.read_csv(ws["dir"] / "out.tsv", sep="\t")


def test_tagged_bam_off_by_default(ecs_reads_workspace):
    ws = ecs_reads_workspace
    tsv = _run_caller(ws)
    assert tsv["indel_reads"].sum() == ws["n_edit_reads"]
    # the feature is opt-in: nothing extra on disk unless asked for
    assert not list(ws["dir"].glob("*.tagged.bam*"))


def test_tagged_bam_tags_every_read_exactly_once(ecs_reads_workspace):
    pysam = pytest.importorskip("pysam")
    ws = ecs_reads_workspace
    out = ws["dir"] / "tagged.bam"
    tsv = _run_caller(ws, "--tagged-bam-out", out)

    assert out.exists() and (ws["dir"] / "tagged.bam.bai").exists()

    keys, tags, positions = [], collections.Counter(), []
    spanb_read1 = []
    with pysam.AlignmentFile(str(out)) as bam:
        assert bam.header.to_dict()["HD"]["SO"] == "coordinate"
        for read in bam:
            assert read.has_tag("XC"), f"{read.query_name} has no XC tag"
            keys.append((read.query_name, read.flag, read.reference_start))
            tags[read.get_tag("XC")] += 1
            positions.append(read.reference_start)
            if read.query_name.startswith("spanb") and read.is_read1:
                spanb_read1.append(read.get_tag("XC"))

    # (a) one record per alignment: the whole point of the second pass
    duplicated = [k for k, n in collections.Counter(keys).items() if n > 1]
    assert not duplicated, f"duplicate records in tagged BAM: {duplicated[:5]}"
    assert len(keys) == ws["total_records"]

    # (b) coordinate-sorted, so IGV will load it
    assert positions == sorted(positions)

    # (c) tag counts reconcile with the TSV. Tags are per alignment record while
    #     indel_reads is per fragment (the caller collapses mates by read name), so
    #     these are only equal because just R1 carries the deletion in this fixture.
    #     On real data the record count runs ~2x the TSV number.
    edit_tag = f"Edited_Deletion_{ws['del_len']}bp"
    assert tags[edit_tag] == ws["n_edit_reads"] == tsv["indel_reads"].sum()

    # (d) skipped reads are visibly skipped, not silently absent
    assert tags["Skipped_Duplicate"] == 2 * C.ECS_N_DUP
    assert tags["Skipped_LowMapQ"] == 2 * C.ECS_N_LOWMAPQ
    assert tags["Skipped_Mismatches"] == 2 * C.ECS_N_MISMATCH
    # reads sitting in the +/-150 bp pad that span neither target: the implicit
    # default, never stored in the tag map (see DEFAULT_TAG)
    assert tags["Skipped_NoSpan"] == 2 * C.ECS_N_OVERLAP + C.ECS_N_SPAN_B

    # (e) conflict resolution: these reads fall in both padded windows, where one
    #     target cannot evaluate them and the other calls them reference
    assert spanb_read1 == ["Unedited_WT"] * C.ECS_N_SPAN_B


def test_tagged_bam_honours_custom_tag_name(ecs_reads_workspace):
    pysam = pytest.importorskip("pysam")
    ws = ecs_reads_workspace
    out = ws["dir"] / "tagged.bam"
    _run_caller(ws, "--tagged-bam-out", out, "--tagged-bam-tag", "YC")
    with pysam.AlignmentFile(str(out)) as bam:
        read = next(iter(bam))
    assert read.has_tag("YC") and not read.has_tag("XC")


def test_read_tag_precedence_is_order_independent():
    """A read seen at two overlapping targets keeps the most specific call."""
    m = _load_find_edited_reads()

    ordered = ["Edited_BND_chr19", "Edited_Deletion_5bp", "Edited_SoftClip",
               "Unedited_WT", "Skipped_NoSpan", "Skipped_LowMapQ"]
    ranks = [m.read_tag_rank(t) for t in ordered]
    assert ranks == sorted(ranks) and len(set(ranks)) == len(ranks)

    class FakeRead:
        query_name, flag, reference_start = "r1", 99, 100

    for first, second in ((("Unedited_WT"), "Edited_Deletion_5bp"),
                          ("Edited_Deletion_5bp", "Unedited_WT"),
                          ("Skipped_NoSpan", "Unedited_WT")):
        tags = {}
        m.record_read_tag(tags, FakeRead(), first)
        m.record_read_tag(tags, FakeRead(), second)
        winner = min([first, second], key=m.read_tag_rank)
        assert tags[("r1", 99, 100)] == winner


def test_classify_read_tag_vocabulary():
    m = _load_find_edited_reads()
    assert m.classify_read_tag(None, "CIGAR") == "Skipped_Unevaluable"
    assert m.classify_read_tag({"alttype": "REF"}, "REF") == "Unedited_WT"
    assert m.classify_read_tag(
        {"alttype": "BND", "chrom2": "PLVM_CD19_CARv4_cd34"}, "SA"
    ) == "Edited_BND_PLVM_CD19_CARv4_cd34"
    assert m.classify_read_tag(
        {"alttype": "DEL", "ref": "ATTTTT", "alt": "A"}, "CIGAR") == "Edited_Deletion_5bp"
    assert m.classify_read_tag(
        {"alttype": "INS", "ref": "A", "alt": "ACCC"}, "CIGAR") == "Edited_Insertion_3bp"
    # symbolic alleles carry no length, so fall back to the breakpoint span
    assert m.classify_read_tag(
        {"alttype": "DUP", "ref": "A", "alt": "<DUP>", "pos": 100, "pos2": 112},
        "SA") == "Edited_Duplication_12bp"
    # soft-clip calls are realignments, reported by mechanism not implied size
    assert m.classify_read_tag(
        {"alttype": "DEL", "ref": "ATTTTT", "alt": "A"}, "SOFTCLIP") == "Edited_SoftClip"


def test_merge_windows_collapses_overlapping_targets():
    m = _load_find_edited_reads()
    merged = m.merge_windows([("chr1", 850, 1151), ("chr1", 1050, 1351),
                              ("chr1", 5000, 5300), ("chr2", 10, 300)])
    assert merged == [("chr1", 850, 1351), ("chr1", 5000, 5300), ("chr2", 10, 300)]
