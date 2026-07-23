// BUILD_TRAINING_TABLE — join WGS hotspot features ⋈ ECS truth on (guide, chrom, start).
// Emits training.tsv (WGS features + ECS VAF + label) for the offline model trainer.
// Training is deliberately NOT in this DAG; the deployed model stays a fixed asset.
// Retrain it with the separate `-entry TRAIN` workflow (workflows/train.nf):
//   -entry TRAIN --input <this training.tsv>  ->  a new wgs_shape_model.pkl.
process BUILD_TRAINING_TABLE {
    tag "cohort"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path wgs_scores
    path truth
    path samplesheet

    output:
    path "training.tsv", emit: training
    path "versions.yml", emit: versions

    script:
    // Invoke the repo script by explicit path, NOT bare name: the container bakes an older
    // copy at /opt/scge-offtarget/bin that would otherwise shadow the repo's bin/ (Nextflow
    // only *appends* the pipeline bin/ to PATH). Explicit `python ${projectDir}/bin/...` is
    // immune to both PATH ordering and the script's executable bit.
    """
    python ${projectDir}/bin/join_training_table.py \\
        --wgs-scores ${wgs_scores} \\
        --truth ${truth} \\
        --samplesheet ${samplesheet} \\
        --germline-max-ctrl-if ${params.offtarget_germline_max_ctrl_if} \\
        --out training.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch training.tsv versions.yml"
}
