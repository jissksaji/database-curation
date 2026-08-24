process MULTIQC {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.31--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0'}"

    input:
    path qc_files
    path config

    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_report_data', emit: data
    path 'versions.yml', emit: versions

    script:
    """
    multiqc --filename multiqc_report.html --config ${config} .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed -E 's/^multiqc, version //')
    END_VERSIONS
    """
}
