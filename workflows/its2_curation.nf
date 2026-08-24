include { MULTIQC } from '../modules/multiqc/main'
include { CUSTOM_DB_FILTER } from '../modules/seqkit/custom_db_filter/main'
include { ITSX } from '../modules/itsx/main'
include { AMPLICON_LENGTH } from '../modules/seqkit/amplicon_length/main'
include { BUILD_DB_TAXIDS } from '../modules/helper/build_db_taxids/main'
include { BUILD_DB_LINEAGES } from '../modules/taxonkit/build_db_lineages/main'
include { FUNGAL_SCREENING } from '../modules/vsearch/fungal_screening/main'
include { DEREPLICATE_PER_SPECIES } from '../modules/vsearch/dereplicate_per_species/main'
include { DEREPLICATE_FULL_LENGTH } from '../modules/vsearch/dereplicate_full_length/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/custom/dumpsoftwareversions/main'
include { INPUT_CHECK } from '../modules/input_check'

workflow ITS2_CURATION {
    main:
    log.info "ITS2 database-curation template (${workflow.manifest.version})"

    if (!params.custom_db) {
        throw new IllegalArgumentException(
            'Please provide a database with --custom_db.'
        )
    }
    if (!params.reference_base) {
        throw new IllegalArgumentException(
            'Please provide the installed taxonomy directory with --reference_base.'
        )
    }

    itsx_mode = params.itsx.toString().toUpperCase().replaceAll(/\s+/, '')
    if (!(itsx_mode in ['T', 'F', 'T,F', 'F,T'])) {
        throw new IllegalArgumentException(
            '--itsx must be T, F, T,F, or F,T.'
        )
    }
    run_fungi = itsx_mode.contains('F')

    dereplication_mode = params.dereplication.toString().toLowerCase()
    if (!(dereplication_mode in ['species', 'full_length'])) {
        throw new IllegalArgumentException(
            '--dereplication must be species or full_length.'
        )
    }

    custom_db_file = file(params.custom_db, checkIfExists: true)
    accession2taxid_file = file(params.references.accession2taxid, checkIfExists: true)
    taxdump_directory = file(params.references.taxdump, type: 'dir', checkIfExists: true)
    fungal_screening_db_file = file(params.fungal_screening_db, checkIfExists: true)

    custom_db_channel = channel.of(
        tuple([id: custom_db_file.baseName], custom_db_file)
    )
    accession2taxid_channel = channel.value(accession2taxid_file)
    taxdump_channel = channel.value(taxdump_directory)
    fungal_screening_db_channel = channel.value(fungal_screening_db_file)

    INPUT_CHECK(custom_db_channel)
    CUSTOM_DB_FILTER(INPUT_CHECK.out.db)
    ITSX(CUSTOM_DB_FILTER.out.fasta, run_fungi)
    AMPLICON_LENGTH(ITSX.out.plants_its2)
    BUILD_DB_TAXIDS(ITSX.out.plants_its2, accession2taxid_channel)

    // Match each accession/taxid table with the FASTA it came from.
    lineage_input = BUILD_DB_TAXIDS.out.accession_taxid
        .map { meta, accession_taxid -> tuple(meta.id, meta, accession_taxid) }
        .join(
            ITSX.out.plants_its2.map { meta, fasta -> tuple(meta.id, fasta) }
        )
        .map { id, meta, accession_taxid, fasta ->
            tuple(meta, accession_taxid, fasta)
        }
    BUILD_DB_LINEAGES(lineage_input, taxdump_channel)
    FUNGAL_SCREENING(BUILD_DB_LINEAGES.out.fasta, fungal_screening_db_channel)

    // Run only the dereplication method selected by the user.
    if (dereplication_mode == 'species') {
        DEREPLICATE_PER_SPECIES(FUNGAL_SCREENING.out.fasta)
        ch_dereplicated_db = DEREPLICATE_PER_SPECIES.out.fasta
        ch_dereplication_versions = DEREPLICATE_PER_SPECIES.out.versions
    } else {
        DEREPLICATE_FULL_LENGTH(FUNGAL_SCREENING.out.fasta)
        ch_dereplicated_db = DEREPLICATE_FULL_LENGTH.out.fasta
        ch_dereplication_versions = DEREPLICATE_FULL_LENGTH.out.versions
    }

    /*
     * Curation modules should mix their QC artefacts and versions into these
     * channels. Keeping the aggregation here makes MultiQC independent of the
     * order in which curation stages are added.
     */
    // Curation modules append their process-level versions and QC artefacts to
    // these channels. The filter is the first curation stage.
    ch_versions = CUSTOM_DB_FILTER.out.versions
        .mix(ITSX.out.versions)
        .mix(AMPLICON_LENGTH.out.versions)
        .mix(BUILD_DB_TAXIDS.out.versions)
        .mix(BUILD_DB_LINEAGES.out.versions)
        .mix(FUNGAL_SCREENING.out.versions)
        .mix(ch_dereplication_versions)
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.collectFile(name: 'collated_versions.yml'))
    ch_qc_files = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
    ch_multiqc_config = channel.fromPath(params.multiqc_config, checkIfExists: true).collect()

    MULTIQC(
        ch_qc_files.collect(),
        ch_multiqc_config
    )

    emit:
    n_filtered_db = CUSTOM_DB_FILTER.out.n_filtered
    // Final curated database: taxonomy added and fungal matches removed.
    custom_db = ch_dereplicated_db
    fungi_db = ITSX.out.fungi_its2
    plants_db = ITSX.out.plants_its2
    filtered_db = CUSTOM_DB_FILTER.out.fasta
    amplicon_lengths = AMPLICON_LENGTH.out.tsv
    stats = AMPLICON_LENGTH.out.stats
    taxids = BUILD_DB_TAXIDS.out.taxids
    taxid_counts = BUILD_DB_TAXIDS.out.taxids_counts
    accession_taxid = BUILD_DB_TAXIDS.out.accession_taxid
    missing_accessions = BUILD_DB_TAXIDS.out.missing
    lineages = BUILD_DB_LINEAGES.out.lineages
    taxonomic_fasta = BUILD_DB_LINEAGES.out.fasta
    fungal_similarity = FUNGAL_SCREENING.out.similarity
    fungal_like = FUNGAL_SCREENING.out.fungal_like
    no_strong_fungal_hit = FUNGAL_SCREENING.out.fasta
    dereplicated_db = ch_dereplicated_db
    qc = MULTIQC.out.html
    software_versions = CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    versions = ch_versions.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.versions)
}
