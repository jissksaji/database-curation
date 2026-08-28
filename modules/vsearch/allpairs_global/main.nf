process VSEARCH_ALLPAIRS_GLOBAL {

    tag "${meta.id}:${family_dir.baseName}"
    label 'medium_parallel'

    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(family_dir)

    output:
    tuple val(meta), path("*.allpairs_global.blast6.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${family_dir.baseName}"
    def acceptAll = params.vsearch_allpairs_acceptall ? '--acceptall' : ''

    """
    printf '%s\\n' 'taxon\tquery\ttarget\tpident\talnlen\tmismatches\tgaps\tqstart\tqend\ttstart\ttend\tevalue\tbits' \
        > "${prefix}.allpairs_global.blast6.tsv"

    # Each FASTA is already one species or one genus fallback group.
    for taxon_fasta in "${family_dir}"/*.fasta; do
        [[ -e "\${taxon_fasta}" ]] || continue

        taxon_name=\$(basename "\${taxon_fasta}" .fasta)

        vsearch \
            --allpairs_global "\${taxon_fasta}" \
            --blast6out /dev/stdout \
            --id ${params.vsearch_allpairs_id} \
            --iddef ${params.vsearch_allpairs_iddef} \
            ${acceptAll} \
            --maxaccepts ${params.vsearch_allpairs_maxaccepts} \
            --maxrejects ${params.vsearch_allpairs_maxrejects} \
            --minwordmatches ${params.vsearch_allpairs_minwordmatches} \
            --qmask ${params.vsearch_allpairs_qmask} \
            --threads ${task.cpus} \
        | awk -v taxon="\${taxon_name}" \
            'BEGIN { OFS = "\t" } { print taxon, \$0 }' \
            >> "${prefix}.allpairs_global.blast6.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/^vsearch v//')
    END_VERSIONS
    """
}
