/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Off-Target Model Trainer  (entry: -entry TRAIN)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Offline trainer for the WGS shape model. Deliberately kept OUT of the OFFTARGET
    DAG (the deployed model stays a fixed asset); this is the explicit retrain loop:

      1. run -entry OFFTARGET (paired ECS+WGS) on a labeled cohort  -> training.tsv
      2. run -entry TRAIN --input training.tsv                       -> wgs_shape_model.pkl
      3. deploy the new model:  -entry OFFTARGET --offtarget_shape_model <new.pkl>

    --input is the training.tsv from step 1 (a BUILD_TRAINING_TABLE output).
*/

include { TRAIN_SHAPE_MODEL } from '../modules/local/train_shape_model.nf'

workflow TRAIN_WORKFLOW {

    if (!params.input) {
        error "TRAIN: provide --input training.tsv (a BUILD_TRAINING_TABLE output from -entry OFFTARGET)"
    }
    ch_training = file(params.input, checkIfExists: true)

    TRAIN_SHAPE_MODEL(ch_training)
}
