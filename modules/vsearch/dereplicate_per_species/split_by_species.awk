# Read the species name from each FASTA header.
/^>/ {
    species = "unclassified"

    if (match($0, /;s:[^;\t]*/)) {
        species = substr($0, RSTART + 3, RLENGTH - 3)
    }

    # Use simple numbered filenames so spaces in species names are safe.
    if (!(species in species_file)) {
        species_count++
        species_file[species] = sprintf("%s/species_%06d.fasta", output_dir, species_count)
    }

    current_file = species_file[species]
}

# Write the header and sequence to its species file.
current_file != "" {
    print $0 >> current_file
    close(current_file)
}
