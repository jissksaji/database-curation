process CUSTOM_DB_FILTER {

    // Filters the custom database using FASTA headers, sequence length,
    // and the maximum allowed number of consecutive N bases.

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.n_filtered.fasta"), emit: n_filtered
    tuple val(meta), path("*.cleaned.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_pattern = params.custom_db_filter_pattern_exclude
    def include_pattern = params.custom_db_include_pattern
    def min_length = params.custom_db_min_length.toString().toInteger()
    def max_length = params.custom_db_max_length.toString().toInteger()
    def max_consecutive_ns = params.custom_db_max_consecutive_ns.toString().toInteger()
    def n_pattern = "N{${max_consecutive_ns + 1}}N*"
    def seqkit_threads = Math.max(1, task.cpus.intdiv(3))
    def max_length_option = max_length > 0 ? "--max-len ${max_length}" : ""

    """
    # First file:
    # Remove unwanted headers, short sequences, and long runs of N bases.
    # The maximum-length and ITS2 inclusion filters are not applied here.
    seqkit grep --by-name --invert-match --use-regexp --ignore-case \\
        --threads ${seqkit_threads} \\
        --pattern "${exclude_pattern}" \\
        ${db} |
        seqkit seq \\
            --threads ${seqkit_threads} \\
            --min-len ${min_length} \\
            |
        seqkit grep --by-seq --invert-match --use-regexp --ignore-case \\
            --threads ${seqkit_threads} \\
            --pattern "${n_pattern}" \\
            --out-file ${prefix}.n_filtered.fasta

    # Second file:
    # Keep ITS2 headers and apply the optional maximum sequence length.
    seqkit grep --by-name --use-regexp --ignore-case \\
        --threads ${seqkit_threads} \\
        --pattern "${include_pattern}" \\
        ${prefix}.n_filtered.fasta |
        seqkit seq \\
            --threads ${seqkit_threads} \\
            ${max_length_option} \\
            --out-file ${prefix}.cleaned.fasta

    # Stop here with a useful message instead of sending an empty file to ITSx.
    if [[ ! -s "${prefix}.cleaned.fasta" ]]; then
        echo "WARNING: ${prefix}.cleaned.fasta is empty after filtering." >&2
        echo "Check the inclusion pattern and the length/N filters." >&2
        echo "For accession-only FASTA headers, use --custom_db_include_pattern '.*'." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
