process ITSX {

    tag "${meta.id}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    val run_fungi

    output:
    tuple val(meta), path("*_plants_itsx.ITS2.fasta"), emit: plants_its2
    tuple val(meta), path("*_plants_itsx.ITS2.full_and_partial.fasta"), emit: plants_its2_partial
    tuple val(meta), path("*_plants_itsx.full.fasta"), emit: plants_full
    tuple val(meta), path("*_plants_itsx.full_and_partial.fasta"), emit: plants_full_and_partial
    tuple val(meta), path("*_plants_itsx.positions.txt"), emit: plants_positions
    tuple val(meta), path("*_plants_itsx.extraction.results"), emit: plants_detailed_results
    tuple val(meta), path("*_plants_itsx.summary.txt"), emit: plants_summary
    tuple val(meta), path("*_plants_itsx_no_detections.fasta"), emit: plants_not_found
    tuple val(meta), path("*_plants_itsx.problematic.txt"), optional: true, emit: plants_problematic

    tuple val(meta), path("*_fungi_itsx.ITS2.fasta"), optional: true, emit: fungi_its2
    tuple val(meta), path("*_fungi_itsx.ITS2.full_and_partial.fasta"), optional: true, emit: fungi_its2_partial
    tuple val(meta), path("*_fungi_itsx.full.fasta"), optional: true, emit: fungi_full
    tuple val(meta), path("*_fungi_itsx.full_and_partial.fasta"), optional: true, emit: fungi_full_and_partial
    tuple val(meta), path("*_fungi_itsx.positions.txt"), optional: true, emit: fungi_positions
    tuple val(meta), path("*_fungi_itsx.extraction.results"), optional: true, emit: fungi_detailed_results
    tuple val(meta), path("*_fungi_itsx.summary.txt"), optional: true, emit: fungi_summary
    tuple val(meta), path("*_fungi_itsx_no_detections.fasta"), optional: true, emit: fungi_not_found
    tuple val(meta), path("*_fungi_itsx.problematic.txt"), optional: true, emit: fungi_problematic

    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fungi_command = run_fungi
        ? """
    ITSx \\
        -i ${fasta} \\
        -o ${prefix}_fungi_itsx \\
        -t F \\
        --cpu ${task.cpus} \\
        --multi_thread T \\
        --save_regions ITS2 \\
        --preserve T \\
        --partial 1 \\
        --positions T \\
        --detailed_results T \\
        --summary T \\
        --graphical F \\
        --not_found T
    """
        : ''

    """
    ITSx \
        -i ${fasta} \
        -o ${prefix}_plants_itsx \
        -t T \
        --cpu ${task.cpus} \
        --multi_thread T \
        --save_regions ITS2 \
        --preserve T \
        --partial 1 \
        --positions T \
        --detailed_results T \
        --summary T \
        --graphical F \
        --not_found T

    ${fungi_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        itsx: \$(ITSx --help 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
