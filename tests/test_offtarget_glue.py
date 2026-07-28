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
