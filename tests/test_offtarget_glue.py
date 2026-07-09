"""Logic tests for the Unified CRISPR Off-Target Workflow glue scripts.

Scenario (see conftest): one guide, two ECS replicates + one WGS sample, three hotspots
(a real edit, a true negative, and a below-floor low-VAF edit). These lock in the
ECS⋈WGS join and the recall-vs-VAF behaviour that the Nextflow layer depends on.
"""
import pandas as pd
from conftest import run


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

    # overall: 1 of 2 ECS-positive sites recovered (chr12 yes @0.15, chr7 no @0.008)
    assert "overall WGS recall of ECS-confirmed edits (score>=0.6): 0.50" in proc.stdout
    m = pd.read_csv(workspace / "recall_vs_vaf.csv")
    recalls = dict(zip(m["vaf_bin"], m["recall"]))
    assert recalls["[0.005, 0.01)"] == 0.0                 # below the depth floor
    assert recalls["[0.1, 0.2)"] == 1.0                    # above it
    assert (workspace / "recall_vs_vaf.png").exists()


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
    """Guardrail: an ECS/WGS coordinate drift yields an empty join, and the script warns."""
    run("hotspot_to_table.py",
        "--ecs-tables", "PLCB2_ecs_1.offtarget_analysis.tsv", "PLCB2_ecs_2.offtarget_analysis.tsv",
        "--samplesheet", "samplesheet.csv", cwd=workspace)
    w = pd.read_csv(workspace / "wgs_hotspot_scores.csv")
    w["start"] += 5                                        # simulate a 5bp drift
    w.to_csv(workspace / "wgs_shifted.csv", index=False)
    proc = run("join_training_table.py", "--wgs-scores", "wgs_shifted.csv",
               "--truth", "ecs_hotspot_truth.csv", "--samplesheet", "samplesheet.csv", cwd=workspace)

    assert "empty join" in proc.stderr
    assert len(pd.read_csv(workspace / "training.tsv", sep="\t")) == 0
