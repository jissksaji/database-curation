BEGIN {
    FS = "\t"
}

# Read the lineage table first and save each complete output header.
FILENAME == ARGV[1] {
    lineage_header[$1] = $0
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
