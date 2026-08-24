/*
 * Copy this directory to create a new process module.
 *
 * Set `conda` to the module's environment.yml and/or `container` to an OCI image.
 * Run the pipeline with -profile conda, docker, singularity, apptainer, or podman.
 */
process MODULE_NAME {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    // container 'quay.io/your-org/your-image:tag'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('output.txt'), emit: output

    script:
    """
    cp '${input_file}' output.txt
    """
}
