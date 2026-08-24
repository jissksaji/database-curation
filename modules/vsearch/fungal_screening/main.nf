process FUNGAL_SCREENING {

    tag "${meta.id}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path unite_db

    output:
    tuple val(meta), path("*.fungal_similarity.tsv"), emit: similarity
    tuple val(meta), path("*.fungal_like.fasta"), emit: fungal_like
    // Main output: the incoming lineage FASTA with fungal matches removed.
    tuple val(meta), path("*.no_strong_fungal_hit.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Screen the lineage FASTA against the UNITE fungal database.
    vsearch \
        --usearch_global "${fasta}" \
        --db "${unite_db}" \
        --id ${params.fungal_screening_identity} \
        --query_cov ${params.fungal_screening_query_cov} \
        --mincols ${params.fungal_screening_mincols} \
        --strand both \
        --maxaccepts 10 \
        --maxrejects 0 \
        --maxhits 1 \
        --userout "${prefix}.fungal_similarity.tsv" \
        --userfields query+target+id+qcov+tcov+alnlen+mism+gaps \
        --threads ${task.cpus}

    # Split the original FASTA so its taxid and lineage headers are preserved.
    cut -f 1 "${prefix}.fungal_similarity.tsv" | sort -u > fungal_ids.txt

    if [[ -s fungal_ids.txt ]]; then
        seqkit grep \
            --pattern-file fungal_ids.txt \
            "${fasta}" \
            --out-file "${prefix}.fungal_like.fasta"

        seqkit grep \
            --invert-match \
            --pattern-file fungal_ids.txt \
            "${fasta}" \
            --out-file "${prefix}.no_strong_fungal_hit.fasta"
    else
        touch "${prefix}.fungal_like.fasta"
        cp "${fasta}" "${prefix}.no_strong_fungal_hit.fasta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/^vsearch v//')
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """
}
