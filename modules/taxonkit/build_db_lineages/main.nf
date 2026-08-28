process BUILD_DB_LINEAGES {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(accession_taxid), path(fasta)
    path ncbi_taxonomy

    output:
    tuple val(meta), path("*.lineages.tsv"), emit: lineages
    tuple val(meta), path("*.lineages.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Add the lineage to the accession and taxid table.
    taxonkit reformat2 \
        --data-dir "${ncbi_taxonomy}" \
        --taxid-field 2 \
        --format 'd:{domain|superkingdom};p:{phylum};c:{class};o:{order};f:{family};g:{genus};s:{species}' \
        "${accession_taxid}" \
        > "${prefix}.lineages.tsv"

    # Join each lineage back to its FASTA sequence. Taxonomy is kept inside a
    # whitespace-free identifier so ITSx preserves the complete header.
    awk -f "${moduleDir}/join_lineages.awk" \
        "${prefix}.lineages.tsv" \
        "${fasta}" \
        > "${prefix}.lineages.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version | head -n 1)
    END_VERSIONS
    """
}
