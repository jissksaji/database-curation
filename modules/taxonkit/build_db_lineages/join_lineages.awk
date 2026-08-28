BEGIN {
    FS = "\t"
}

# Read the lineage table first and build a whitespace-free FASTA identifier.
# ITSx preserves the complete identifier, so taxonomy remains attached.
FILENAME == ARGV[1] {
    accession = $1
    taxid = $2
    lineage = $3
    gsub(/[[:space:]]+/, "_", lineage)
    lineage_header[accession] = accession ";taxid=" taxid ";" lineage
    next
}

# Find the accession at the start of each FASTA header.
/^>/ {
    header = $0
    sub(/^>/, "", header)
    split(header, fields, /[[:space:];]+/)

    accession = fields[1]
    sub(/(\.[0-9]+)+$/, "", accession)

    keep_sequence = accession in lineage_header
    if (keep_sequence) {
        print ">" lineage_header[accession]
    }
    next
}

# Print sequence lines only when their accession has a lineage.
keep_sequence {
    print
}
