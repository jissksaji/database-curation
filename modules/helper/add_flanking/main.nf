process ADD_FLANKING {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(lineage_fasta), path(positions), path(source_fasta)

    output:
    tuple val(meta), path("*.flanked.fasta"), emit: fasta
    tuple val(meta), path("*.flanking_missing.txt"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forward = params.add_flanking_forward.toString().toInteger()
    def reverse = params.add_flanking_reverse.toString().toInteger()

    """
    add_flanking.py \
        --its2-fasta "${lineage_fasta}" \
        --positions "${positions}" \
        --source-fasta "${source_fasta}" \
        --forward ${forward} \
        --reverse ${reverse} \
        --output "${prefix}.flanked.fasta" \
        --missing "${prefix}.flanking_missing.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """
}
