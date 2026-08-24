process BUILD_DB_TAXIDS {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path genbank2taxid

    output:
    tuple val(meta), path("*.db_taxids.tsv"),        emit: taxids
    tuple val(meta), path("*.db_taxids_counts.tsv"), emit: taxids_counts
    tuple val(meta), path("*.accession_taxid.tsv"),  emit: accession_taxid
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    # Get one accession from each FASTA header.
    # This supports normal headers (>ABC123.1 description) and SINTAX headers
    # (>ABC123.1;tax=k:...,s:...;).
    awk '/^>/ {
        header = \$0
        sub(/^>/, "", header)
        split(header, fields, /[[:space:];]+/)
        print fields[1]
    }' "${fasta}" > db_accessions.txt

    # Read only the database accessions into memory, then scan the much larger
    # accession-to-taxid file once.
    if [[ -s db_accessions.txt ]]; then
        awk -F'\\t' -v accession_taxid_file="${prefix}.accession_taxid.tsv" \
            -v missing_file="${prefix}.missing_accessions.tsv" '
            FILENAME == ARGV[1] {
                sub(/(\\.[0-9]+)+\$/, "", \$1)
                database_accessions[\$1] = 1
                next
            }

            {
                accession = \$1
                sub(/(\\.[0-9]+)+\$/, "", accession)

                if (accession in database_accessions) {
                    # NCBI accession2taxid files store the taxid in column 3.
                    # A simple two-column mapping stores it in column 2.
                    taxid = (NF >= 3) ? \$3 : \$2

                    print accession "\\t" taxid > accession_taxid_file
                    print taxid
                    delete database_accessions[accession]
                }
            }

            END {
                for (accession in database_accessions) {
                    print accession > missing_file
                }
            }
        ' db_accessions.txt "${genbank2taxid}" |
            sort -n |
            tee >(sort -u > "${prefix}.db_taxids.tsv") |
            uniq -c |
            awk '{print \$2 "\\t" \$1}' \
            > "${prefix}.db_taxids_counts.tsv"
    fi

    # All outputs must exist, including when ITSx found no plant sequences.
    touch "${prefix}.db_taxids.tsv"
    touch "${prefix}.db_taxids_counts.tsv"
    touch "${prefix}.accession_taxid.tsv"
    touch "${prefix}.missing_accessions.tsv"

    rm db_accessions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
