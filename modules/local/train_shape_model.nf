// TRAIN_SHAPE_MODEL — offline trainer for the off-target WGS shape model.
// Reads a training.tsv (BUILD_TRAINING_TABLE output) and fits the Stage-2 shape
// ranker consumed by score.py / worklist_from_vcf.py. Runs only under -entry TRAIN;
// the deployed OFFTARGET pipeline keeps the model as a fixed asset by design.
process TRAIN_SHAPE_MODEL {
    tag "train"
    label 'process_low'
    container "ghcr.io/dhslab/docker-scge-offtarget:260710"

    input:
    path training

    output:
    path "wgs_shape_model.pkl", emit: model
    path "train_metrics.json",  emit: metrics
    path "versions.yml",        emit: versions

    script:
    // Explicit python ${projectDir}/bin path: immune to PATH ordering and the exec bit,
    // and to any stale copy the container may bake in.
    """
    python ${projectDir}/bin/train_shape_model.py \\
        --training ${training} \\
        --learning-rate ${params.offtarget_train_learning_rate} \\
        --max-iter ${params.offtarget_train_max_iter} \\
        --max-depth ${params.offtarget_train_max_depth} \\
        --holdout-frac ${params.offtarget_train_holdout_frac} \\
        --seed ${params.offtarget_train_seed} \\
        --out wgs_shape_model.pkl \\
        --metrics train_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    "touch wgs_shape_model.pkl train_metrics.json versions.yml"
}
