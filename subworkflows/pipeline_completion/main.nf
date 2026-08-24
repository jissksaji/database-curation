workflow PIPELINE_COMPLETION {
    workflow.onComplete = {
        def outputDir = new File("${params.outdir}/pipeline_info")
        outputDir.mkdirs()
        new File(outputDir, 'pipeline_report.txt').text = """version: ${workflow.manifest.version}
session: ${workflow.sessionId}
success: ${workflow.success}
commandLine: ${workflow.commandLine}
"""
    }
}
