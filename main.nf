#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ITS2_CURATION } from './workflows/its2_curation'
include { BUILD_REFERENCE } from './workflows/build_reference'
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_completion'

workflow {
    WorkflowMain.initialise(workflow, params, log)
    WorkflowPipeline.initialise(params, log)

    if (params.custom_db) {
        ITS2_CURATION()
    } else if (params.reference_base) {
        BUILD_REFERENCE()
    } else {
        throw new IllegalArgumentException(
            'Please provide --custom_db for curation or --reference_base to install references.'
        )
    }

    PIPELINE_COMPLETION()
}
