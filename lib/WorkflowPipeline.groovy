class WorkflowPipeline {
    static void initialise(params, log) {
        if (!params.custom_db && !params.reference_base) {
            log.info 'No custom ITS2 database supplied.'
        }
    }
}
