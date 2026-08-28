process DEREPLICATE_PER_SPECIES {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.taxon_derep.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir taxon_fastas

    # Put records with the same species in one FASTA. If species is missing,
    # put records with the same genus in one FASTA instead.
    awk \
        -v output_dir="taxon_fastas" \
        -f "${moduleDir}/../split_by_taxon.awk" \
        "${fasta}"

    # Dereplicate every species/genus FASTA separately, then combine them.
    touch "${prefix}.taxon_derep.fasta"

    for taxon_fasta in taxon_fastas/*.fasta; do
        [[ -e "\${taxon_fasta}" ]] || continue

        vsearch \
            --derep_fulllength "\${taxon_fasta}" \
            --output "\${taxon_fasta}.derep.fasta" \
            --sizeout \
            --notrunclabels

        cat "\${taxon_fasta}.derep.fasta" >> "${prefix}.taxon_derep.fasta"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/^vsearch v//')
    END_VERSIONS
    """
}
