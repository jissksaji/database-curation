process MMSEQS_EASY_SEARCH {

    tag "${meta.id}:${query_fasta.baseName}"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(query_fasta)
    tuple path(search_database), path(lineage_tsv)

    output:
    tuple val(meta), path("*.mmseqs_easy_search.tsv"), emit: results
    path "versions.yml", emit: versions

    script:
    def prefix = query_fasta.baseName.replaceAll(/\.outliers$/, '')

    """
    if [[ -s "${query_fasta}" ]]; then
        mmseqs easy-search \
            "${query_fasta}" \
            "${search_database}" \
            mmseqs_hits.raw.tsv \
            mmseqs_tmp \
            --search-type ${params.mmseqs_easy_search_search_type} \
            --min-seq-id ${params.mmseqs_easy_search_min_seq_id} \
            --min-aln-len ${params.mmseqs_easy_search_min_aln_len} \
            -c ${params.mmseqs_easy_search_coverage} \
            --cov-mode ${params.mmseqs_easy_search_cov_mode} \
            --alignment-mode ${params.mmseqs_easy_search_alignment_mode} \
            -e ${params.mmseqs_easy_search_evalue} \
            -s ${params.mmseqs_easy_search_sensitivity} \
            --max-seqs ${params.mmseqs_easy_search_max_seqs} \
            --threads ${task.cpus} \
            --format-mode 4 \
            --format-output 'query,target,pident,alnlen,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,evalue,bits' \
            --remove-tmp-files 1
    else
        printf '%s\n' \
            'query\ttarget\tpident\talnlen\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\tevalue\tbits' \
            > mmseqs_hits.raw.tsv
    fi

    format_mmseqs_results.py \
        --input mmseqs_hits.raw.tsv \
        --query-fasta "${query_fasta}" \
        --lineages "${lineage_tsv}" \
        --output "${prefix}.mmseqs_easy_search.tsv"

    printf '"%s":\n    mmseqs2: %s\n' \
        "${task.process}" \
        "\$(mmseqs version)" \
        > versions.yml
    """
}
