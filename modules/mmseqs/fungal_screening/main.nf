process FUNGAL_SCREENING {

    tag "${meta.id}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(flanked_fasta), path(lineage_fasta)
    path unite_db

    output:
    tuple val(meta), path("*.all_passing_hits.tsv"), emit: hits
    tuple val(meta), path("*.fungal_like.fasta"), emit: fungal_like
    tuple val(meta), path("*.no_strong_fungal_hit.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def searchType = params.fungal_screening_search_type.toString().toInteger()
    def minSeqId = params.fungal_screening_min_seq_id.toString().toDouble()
    def minAlnLen = params.fungal_screening_min_aln_len.toString().toInteger()
    def evalue = params.fungal_screening_evalue.toString().toDouble()
    def sensitivity = params.fungal_screening_sensitivity.toString().toDouble()
    def maxSeqs = params.fungal_screening_max_seqs.toString().toInteger()

    """
    # Search the flanked ITS2 sequences against the UNITE fungal database.
    mmseqs easy-search \
        "${flanked_fasta}" \
        "${unite_db}" \
        "${prefix}.all_passing_hits.tsv" \
        mmseqs_tmp \
        --search-type ${searchType} \
        --min-seq-id ${minSeqId} \
        --min-aln-len ${minAlnLen} \
        -e ${evalue} \
        -s ${sensitivity} \
        --max-seqs ${maxSeqs} \
        --threads ${task.cpus} \
        --format-mode 4 \
        --format-output 'query,target,pident,alnlen,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,evalue,bits' \
        --remove-tmp-files 1

    # TODO: When primer-binding sites are added, keep only fungal hits that
    # overlap those sites. For now, any hit passing the search thresholds is used.
    # MMseqs format mode 4 writes a header. Collect each query with a hit.
    awk -F '\t' 'NR > 1 { print \$1 }' "${prefix}.all_passing_hits.tsv" \
        | sort -u \
        > fungal_ids.txt

    # Apply the hit IDs to the original lineage ITS2 FASTA. Flanking bases are
    # used only for the search and are not passed into dereplication.
    if [[ -s fungal_ids.txt ]]; then
        seqkit grep \
            --pattern-file fungal_ids.txt \
            "${lineage_fasta}" \
            --out-file "${prefix}.fungal_like.fasta"

        seqkit grep \
            --invert-match \
            --pattern-file fungal_ids.txt \
            "${lineage_fasta}" \
            --out-file "${prefix}.no_strong_fungal_hit.fasta"
    else
        touch "${prefix}.fungal_like.fasta"
        cp "${lineage_fasta}" "${prefix}.no_strong_fungal_hit.fasta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version)
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """
}
