process SCREEN_TAXA {

    tag "${meta.id}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.below_minimum.fasta"), emit: low_count
    tuple val(meta), path("family_groups/*"), emit: families, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    screen_taxa.py \
        --input "${fasta}" \
        --minimum ${params.screen_taxa_min_sequences} \
        --low-output "${prefix}.below_minimum.fasta" \
        --family-dir family_groups

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
