process DEREPLICATE_FULL_LENGTH {

    tag "${meta.id}"
    label 'medium_parallel'

    // Reuse the VSEARCH environment already used by fungal screening.
    conda "${moduleDir}/../fungal_screening/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.derep.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Dereplicate exact full-length sequences without species grouping.
    vsearch \
        --derep_fulllength "${fasta}" \
        --output "${prefix}.derep.fasta" \
        --sizeout \
        --notrunclabels \
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/^vsearch v//')
    END_VERSIONS
    """
}
