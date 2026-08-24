process DEREPLICATE_PER_SPECIES {

    tag "${meta.id}"
    label 'medium_parallel'

    // Reuse the VSEARCH environment already used by fungal screening.
    conda "${moduleDir}/../fungal_screening/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.species_derep.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir species_fastas

    # Make one FASTA file for every s: species in the header.
    awk \
        -v output_dir="species_fastas" \
        -f "${moduleDir}/split_by_species.awk" \
        "${fasta}"

    # Dereplicate each species and add it to the final FASTA.
    touch "${prefix}.species_derep.fasta"

    for species_fasta in species_fastas/*.fasta; do
        [[ -e "\${species_fasta}" ]] || continue

        vsearch \
            --derep_fulllength "\${species_fasta}" \
            --output "\${species_fasta}.derep.fasta" \
            --sizeout \
            --notrunclabels \
            --threads ${task.cpus}

        cat "\${species_fasta}.derep.fasta" \
            >> "${prefix}.species_derep.fasta"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/^vsearch v//')
    END_VERSIONS
    """
}
