"""Logic tests for the Unified CRISPR Off-Target Workflow glue scripts.

Scenario (see conftest): one guide, two ECS replicates + one WGS sample, three hotspots
(a real edit, a true negative, and a below-floor low-VAF edit). These lock in the
ECS⋈WGS join and the recall-vs-VAF behaviour that the Nextflow layer depends on.
"""
import subprocess
import sys

import pandas as pd
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
