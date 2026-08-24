process AMPLICON_LENGTH {

    tag "${meta.id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.amplicon_lengths.tsv"), emit: tsv
    tuple val(meta), path("*.stats.txt"),            emit: stats
    path "versions.yml",                            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit fx2tab --length --name ${fasta} \
        > ${prefix}.amplicon_lengths.tsv

    seqkit stats --all --tabular ${fasta} \
        > ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
