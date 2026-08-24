class WorkflowMain {
    static void initialise(workflow, params, log) {
        log.info header(workflow)
    }

    static String header(workflow) {
        """
================================================================================
${workflow.manifest.description} | version ${workflow.manifest.version}
================================================================================
"""
    }
}
