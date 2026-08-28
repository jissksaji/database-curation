process ZSCORE_SCREENING {

    tag "${meta.id}:${pairwise_tsv.baseName}"
    label 'short_serial'

    input:
    tuple val(meta), path(pairwise_tsv), path(fasta_directory)

    output:
    tuple val(meta), path("*.zscore.tsv"), emit: scores
    tuple val(meta), path("*.outliers.fasta"), emit: outliers
    tuple val(meta), path("*.non_outliers.fasta"), emit: non_outliers
    path "versions.yml", emit: versions

    script:
    def prefix = pairwise_tsv.baseName.replaceAll(/\.allpairs_global\.blast6$/, '')

    """
    zscore_screening.py \
        --pairwise "${pairwise_tsv}" \
        --fasta-directory "${fasta_directory}" \
        --output-prefix "${prefix}" \
        --z-threshold ${params.zscore_threshold} \
        --min-sequences ${params.zscore_min_sequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
