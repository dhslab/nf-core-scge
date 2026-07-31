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
# ── offtarget_metrics.py: PR-AUC / F-beta on a hand-built frame ────────────────
# Every expected value below is computable with a calculator, so a regression shows up
# as a wrong number rather than a plausible-looking one.
#
#   11 credible positives (ecs_if >= 0.005, ecs_indel_reads >= 5), one of which has NO
#      score (INSUFFICIENT COVERAGE) and is excluded from the ranking metrics
#    5 ECS-negatives (label 0, ecs_is_edit 0)
#    2 germline-demoted rows (label 0, ecs_is_edit 1) -> excluded from the DEFAULT
#      negative set; the shape ranker is not the germline filter
#    3 sub-credibility label==1 rows -> AMBIGUOUS, excluded, never counted as negatives
#
# Of the 10 scored positives, 4 are called LIKELY EDIT; 1 of the 5 negatives is.
#   precision = 4/5  = 0.8      recall = 4/10 = 0.4
#   F1 = 2(.8)(.4)/(.8+.4)            = 0.5333...
#   F2 = 5(.8)(.4)/(4(.8)+.4)         = 1.6/3.6  = 0.4444...
#   F5 = 26(.8)(.4)/(25(.8)+.4)       = 8.32/20.4 = 0.40784...
#   recall_incl_unevaluable = 4/11    = 0.36363...
# Scores separate the two classes perfectly, so PR-AUC and ROC-AUC are exactly 1.0.
METRIC_COLS = ["sample", "guide", "chrom", "start", "score", "verdict",
               "ecs_if", "ecs_is_edit", "ecs_indel_reads", "label"]

_EDIT = "LIKELY EDIT"
_ART = "ARTIFACT (shape)"


def _metrics_frame():
    rows = []
    pos_scores = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45]
    for i, sc in enumerate(pos_scores):
        rows.append(("s1", "G", f"chr{i + 1}", 1000 + i, sc,
                     _EDIT if i < 4 else _ART, 0.20, 1, 50, 1))
    # 11th credible positive: no WGS coverage -> no score, out of the ranking metrics
    rows.append(("s1", "G", "chrU", 9000, float("nan"), "INSUFFICIENT COVERAGE",
                 0.20, 1, 50, 1))
    neg_scores = [0.40, 0.30, 0.20, 0.10, 0.05]
    for i, sc in enumerate(neg_scores):
        rows.append(("s1", "G", f"chrN{i}", 2000 + i, sc,
                     _EDIT if i == 0 else _ART, 0.0, 0, 0, 0))
    for i in range(2):                       # germline-demoted: real indel, in normal
        rows.append(("s1", "G", f"chrG{i}", 3000 + i, 0.95,
                     "GERMLINE/ARTIFACT (in normal)", 0.30, 1, 90, 0))
    for i in range(3):                       # sub-credibility ECS signal -> ambiguous
        rows.append(("s1", "G", f"chrA{i}", 4000 + i, 0.10, _ART, 0.001, 1, 2, 1))
    return pd.DataFrame(rows, columns=METRIC_COLS)


def _run_metrics(tmp_path, df, *extra):
    df.to_csv(tmp_path / "training.tsv", sep="\t", index=False)
    proc = run("offtarget_metrics.py", "--training", "training.tsv",
               "--out-json", "m.json", "--out-txt", "m.txt", *extra, cwd=tmp_path)
    import json
    return json.loads((tmp_path / "m.json").read_text()), (tmp_path / "m.txt").read_text(), proc


def test_offtarget_metrics_exact_values(tmp_path):
    out, txt, _ = _run_metrics(tmp_path, _metrics_frame())
    m, op = out["metrics"], out["metrics"]["operating_point"]

    # denominator composition
    assert (m["n_pos"], m["n_neg"]) == (10, 5)          # NaN-score positive not ranked
    assert m["n_excluded_no_score"] == 1
    assert m["n_excluded_no_score_pos"] == 1
    assert out["counts"]["n_positives_credible"] == 11
    assert out["counts"]["n_ambiguous"] == 3            # excluded, NOT negatives
    assert out["counts"]["n_germline_demoted_excluded"] == 2

    # perfect ranking -> PR-AUC and ROC-AUC are exactly 1
    assert m["pr_auc"] == pytest.approx(1.0)
    assert m["roc_auc"] == pytest.approx(1.0)
    assert m["prevalence"] == pytest.approx(10 / 15)

    # operating point, hand-computed
    assert (op["tp"], op["fp"], op["fn"], op["tn"]) == (4, 1, 6, 4)
    assert op["precision"] == pytest.approx(0.8)
    assert op["recall"] == pytest.approx(0.4)
    assert op["f1"] == pytest.approx(2 * 0.8 * 0.4 / (0.8 + 0.4))
    assert op["f2"] == pytest.approx(1.6 / 3.6)
    assert op["f5"] == pytest.approx(8.32 / 20.4)
    assert op["recall_incl_unevaluable"] == pytest.approx(4 / 11)

    # F_beta is monotone in beta: with precision > recall it DECREASES as beta rises, so
    # F1 is an endpoint. F1 sitting in the middle means the beta wiring is inverted.
    assert op["f1"] > op["f2"] > op["f5"]

    # both denominators must be named in the human-readable report
    assert "manual review" in txt.lower()
    assert "UNREVIEWED" in txt
    assert "NOT vs human review" in txt


def test_offtarget_metrics_fbeta_favours_recall(tmp_path):
    """Mirror image: when recall > precision, F_beta RISES with beta."""
    df = _metrics_frame()
    # call every scored positive plus all 5 negatives -> P = 10/15, R = 10/10
    df["verdict"] = df["verdict"].where(df["verdict"] == "INSUFFICIENT COVERAGE", _EDIT)
    # drop the 3 sub-credibility rows only (the ECS-negatives also have a low ecs_if,
    # so filtering on ecs_if alone would silently take the negatives with them)
    df = df[~((df["label"] == 1) & (df["ecs_indel_reads"] < 5))]
    out, _txt, _ = _run_metrics(tmp_path, df, "--negatives", "all_label0")
    op = out["metrics"]["operating_point"]
    assert op["recall"] == pytest.approx(1.0)
    assert op["precision"] == pytest.approx(10 / 17)    # 5 ECS-neg + 2 germline as FP
    assert op["f1"] < op["f2"] < op["f5"]


def test_offtarget_metrics_negative_set_switch(tmp_path):
    """The germline-demoted rows are excluded by default and included on request; the
    alternative is always reported as a labelled sensitivity block either way."""
    default, _t, _ = _run_metrics(tmp_path, _metrics_frame())
    strict, _t2, _ = _run_metrics(tmp_path, _metrics_frame(), "--negatives", "all_label0")
    assert default["metrics"]["n_neg"] == 5
    assert strict["metrics"]["n_neg"] == 7
    assert any("all label==0" in s["definition"] for s in default["sensitivity"])


def test_offtarget_metrics_ambiguous_never_silently_negative(tmp_path):
    """Sub-credibility label==1 rows are the training-table analogue of the unreviewed
    manual-review rows: excluded by default, counted as negatives only on request."""
    default, _t, _ = _run_metrics(tmp_path, _metrics_frame())
    forced, _t2, _ = _run_metrics(tmp_path, _metrics_frame(), "--ambiguous-as-negative")
    assert default["metrics"]["n_neg"] == 5
    assert forced["metrics"]["n_neg"] == 8              # + the 3 ambiguous rows
    assert forced["counts"]["n_ambiguous"] == 0


@pytest.mark.parametrize("mutate,expect", [
    (lambda d: d[d["label"] == 1], "single-class"),          # no negatives at all
    (lambda d: d.assign(ecs_if=0.0, ecs_indel_reads=0), "single-class"),  # no credible pos
])
def test_offtarget_metrics_degenerate_exits_zero(tmp_path, mutate, expect):
    """A run with no credible positives is a legitimate result, not an error: the file
    still has to be published, with nulls and an explanatory note."""
    out, txt, proc = _run_metrics(tmp_path, mutate(_metrics_frame()))
    assert proc.returncode == 0
    assert out["metrics"] is None
    assert expect in out["note"]
    assert "manual review" in txt.lower()                # denominators still documented


def test_offtarget_metrics_all_nan_scores_exits_zero(tmp_path):
    """All-NaN scores must not be filled with 0 — that would fabricate confident
    negatives. Ranking metrics go null, the exclusion count is reported."""
    df = _metrics_frame()
    df["score"] = float("nan")
    out, _txt, proc = _run_metrics(tmp_path, df)
    assert proc.returncode == 0
    assert out["metrics"]["pr_auc"] is None
    assert out["metrics"]["n_excluded_no_score"] == 16   # 11 pos + 5 neg


# --------------------------------------------------------------------------
# Reference-context + cut-site features (added after auditing the model against
# the manual-review rules). These are pure functions over a sequence / record
# list, so they are unit-testable without CRAMs.
# --------------------------------------------------------------------------
sys.path.insert(0, str(BIN))


class _FakeFasta:
    """Minimal pysam.FastaFile stand-in: fetch(chrom, start, end) -> str."""

    def __init__(self, seq):
        self.seq = seq

    def fetch(self, chrom, start, end):
        return self.seq[max(0, start):end]


def test_homopolymer_run_only_counts_runs_touching_the_site():
    from features import homopolymer_run
    # a run AT the site is reported
    assert homopolymer_run("CGAAAAAGC", 4) == 5
    # the same run far from the site is NOT: a homopolymer 8bp away does not
    # explain an indel here, and counting it would flag every read near any repeat
    assert homopolymer_run("AAAAACGCGCG", 9) == 1
    # adjacency counts (a run ending one base before the site still explains a slip)
    assert homopolymer_run("AAAAACG", 5) == 5
    assert homopolymer_run("ACGTACGT", 4) == 1
    assert homopolymer_run("", 0) == 0


def test_repeat_context_flags_homopolymer_and_tandem_repeat():
    from features import repeat_context
    hp = repeat_context(_FakeFasta("ACGT" * 10 + "A" * 10 + "ACGT" * 10), "c", 45)
    assert hp["homopolymer_len"] >= 10
    # a pure trinucleotide repeat is fully covered but has no homopolymer
    cag = repeat_context(_FakeFasta("CAG" * 30), "c", 45)
    assert cag["repeat_frac"] == 1.0
    assert cag["homopolymer_len"] == 1


def test_repeat_context_returns_zeros_off_contig():
    """A missing/unreadable reference must yield zeros, never fabricated signal."""
    from features import repeat_context

    class Boom:
        def fetch(self, *a):
            raise KeyError("no such contig")

    assert repeat_context(Boom(), "chrZ", 100) == {"homopolymer_len": 0, "repeat_frac": 0.0}


def _rec(spans=True, indel=None, mapq=60, softclip=()):
    return {"spans": spans, "indel": indel, "mapq": mapq, "softclip": list(softclip)}


def test_cut_dist_is_distance_from_observed_indel_to_predicted_cut():
    from features import features_from_records
    # 12 reads carrying a 5bp deletion at ref position 1000, plus 8 clean reads
    recs = [_rec(indel=(1000, 5)) for _ in range(12)] + [_rec() for _ in range(8)]
    f = features_from_records(recs, min_span=8, cut_pos=1002)
    assert f["modal_pos"] == 1000
    assert f["cut_dist"] == 2          # |1000 - 1002|
    f_far = features_from_records(recs, min_span=8, cut_pos=1071)
    assert f_far["cut_dist"] == 71


def test_cut_dist_is_nan_not_zero_when_no_indel_is_observed():
    """0 would assert the indel sits exactly on the cut. There is no indel at all."""
    import math
    from features import features_from_records
    f = features_from_records([_rec() for _ in range(20)], min_span=8, cut_pos=1000)
    assert f["modal_pos"] is None
    assert math.isnan(f["cut_dist"])
    # and NaN when no cut site was supplied, rather than a silent 0
    f2 = features_from_records([_rec(indel=(1000, 5)) for _ in range(20)], min_span=8)
    assert math.isnan(f2["cut_dist"])


def test_new_features_are_declared_model_inputs():
    from features import MODEL_FEATURES
    for f in ("cut_dist", "homopolymer_len", "repeat_frac"):
        assert f in MODEL_FEATURES, f"{f} must be a model input, not just a report column"
