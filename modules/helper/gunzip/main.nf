process GUNZIP {

    tag "${meta.id}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path('extracted/*'), emit: files

    script:
    archive_name = archive.name.toLowerCase()
    output_name = archive.name.replaceFirst(/(?i)\.gz$/, '')

    if (!archive_name.endsWith('.zip') && !archive_name.endsWith('.gz')) {
        throw new IllegalArgumentException(
            "GUNZIP only supports .zip and .gz files: ${archive.name}"
        )
    }

    """
    mkdir -p extracted

    if [[ "${archive_name}" == *.zip ]]; then
        unzip -q "${archive}" -d extracted
    else
        gzip -dc "${archive}" > "extracted/${output_name}"
    fi
    """
}
