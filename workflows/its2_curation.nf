include { MULTIQC } from '../modules/multiqc/main'
include { CUSTOM_DB_FILTER } from '../modules/seqkit/custom_db_filter/main'
include { ITSX } from '../modules/itsx/main'
include { AMPLICON_LENGTH } from '../modules/seqkit/amplicon_length/main'
include { BUILD_DB_TAXIDS } from '../modules/helper/build_db_taxids/main'
include { BUILD_DB_LINEAGES } from '../modules/taxonkit/build_db_lineages/main'
include { ADD_FLANKING } from '../modules/helper/add_flanking/main'
include { FUNGAL_SCREENING } from '../modules/mmseqs/fungal_screening/main'
include { DEREPLICATE_PER_SPECIES } from '../modules/vsearch/dereplicate_per_species/main'
include { DEREPLICATE_FULL_LENGTH } from '../modules/vsearch/dereplicate_full_length/main'
include { SCREEN_TAXA } from '../modules/helper/screen_taxa/main'
include { VSEARCH_ALLPAIRS_GLOBAL } from '../modules/vsearch/allpairs_global/main'
include { ZSCORE_SCREENING } from '../modules/custom/zscore_screening/main'
include { MMSEQS_EASY_SEARCH } from '../modules/mmseqs/easy_search/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/custom/dumpsoftwareversions/main'
include { INPUT_CHECK } from '../modules/input_check'

workflow ITS2_CURATION {
    main:
    log.info "ITS2 database-curation template (${workflow.manifest.version})"

    /* Validate user options. */
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

    /* Resolve input files and create the initial channels. */
    custom_db_file = file(params.custom_db, checkIfExists: true)
    custom_db_meta = [id: custom_db_file.baseName]
    accession2taxid_file = file(
        params.references.accession2taxid,
        checkIfExists: true
    )
    taxdump_directory = file(
        params.references.taxdump,
        type: 'dir',
        checkIfExists: true
    )

    fungal_screening_db_path = params.fungal_screening_db.toString()
    if (
        fungal_screening_db_path ==
        'assets/sh_general_release_dynamic_19.02.2025.fasta'
    ) {
        fungal_screening_db_path = "${projectDir}/${fungal_screening_db_path}"
    }
    fungal_screening_db_file = file(
        fungal_screening_db_path,
        checkIfExists: true
    )

    ch_custom_db = channel.of(tuple(custom_db_meta, custom_db_file))
    ch_accession2taxid = channel.value(accession2taxid_file)
    ch_taxdump = channel.value(taxdump_directory)
    ch_fungal_screening_db = channel.value(fungal_screening_db_file)

    /* Filter the source database and attach taxonomy once. */
    INPUT_CHECK(ch_custom_db)
    ch_input_db = INPUT_CHECK.out.db
    CUSTOM_DB_FILTER(ch_input_db)

    // The full accession-to-taxid table also covers all later MMseqs2 hits.
    BUILD_DB_TAXIDS(ch_input_db, ch_accession2taxid)

    ch_lineage_input = BUILD_DB_TAXIDS.out.accession_taxid
        .join(CUSTOM_DB_FILTER.out.fasta)
    BUILD_DB_LINEAGES(ch_lineage_input, ch_taxdump)

    /* Extract ITS2 and remove strong fungal matches. */
    ITSX(BUILD_DB_LINEAGES.out.fasta, run_fungi)
    AMPLICON_LENGTH(ITSX.out.plants_its2)

    ch_flanking_input = ITSX.out.plants_its2
        .join(ITSX.out.plants_positions)
        .join(BUILD_DB_LINEAGES.out.fasta)
    ADD_FLANKING(ch_flanking_input)

    ch_fungal_screening_input = ADD_FLANKING.out.fasta
        .join(ITSX.out.plants_its2)
    FUNGAL_SCREENING(ch_fungal_screening_input, ch_fungal_screening_db)

    /* Dereplicate and prepare family-level comparison groups. */
    DEREPLICATE_PER_SPECIES(FUNGAL_SCREENING.out.fasta)
    SCREEN_TAXA(DEREPLICATE_PER_SPECIES.out.fasta)

    if (dereplication_mode == 'species') {
        ch_dereplicated_db = DEREPLICATE_PER_SPECIES.out.fasta
        ch_dereplication_versions = DEREPLICATE_PER_SPECIES.out.versions
    } else {
        DEREPLICATE_FULL_LENGTH(DEREPLICATE_PER_SPECIES.out.fasta)
        ch_dereplicated_db = DEREPLICATE_FULL_LENGTH.out.fasta
        ch_dereplication_versions = DEREPLICATE_PER_SPECIES.out.versions
            .mix(DEREPLICATE_FULL_LENGTH.out.versions)
    }

    // SCREEN_TAXA can emit one directory or a list; normalize to one per item.
    ch_family_groups = SCREEN_TAXA.out.families
        .flatMap { meta, family_dirs ->
            def directory_list = family_dirs instanceof List
                ? family_dirs
                : [family_dirs]
            directory_list.collect { family_dir -> tuple(meta, family_dir) }
        }

    /* Detect sequence outliers within each taxonomic group. */
    VSEARCH_ALLPAIRS_GLOBAL(ch_family_groups)

    // Key both streams by database and family before joining them.
    ch_pairwise_by_family = VSEARCH_ALLPAIRS_GLOBAL.out.tsv
        .map { meta, pairwise_tsv ->
            def family_name = pairwise_tsv.name.replaceFirst(
                /\.allpairs_global\.blast6\.tsv$/,
                ''
            )
            tuple("${meta.id}:${family_name}", meta, pairwise_tsv)
        }
    ch_fastas_by_family = ch_family_groups
        .map { meta, family_dir ->
            tuple("${meta.id}:${family_dir.baseName}", family_dir)
        }
    ch_zscore_input = ch_pairwise_by_family
        .join(ch_fastas_by_family)
        .map { _key, meta, pairwise_tsv, family_dir ->
            tuple(meta, pairwise_tsv, family_dir)
        }
    ZSCORE_SCREENING(ch_zscore_input)

    /* Search all outliers once against the original database. */
    ch_mmseqs_queries = ZSCORE_SCREENING.out.outliers
        .map { _meta, outliers -> outliers }
        .collectFile(
            name: "${custom_db_file.baseName}.outliers.fasta",
            newLine: true
        )
        .map { outliers -> tuple(custom_db_meta, outliers) }
    ch_mmseqs_reference = ch_input_db
        .join(BUILD_DB_LINEAGES.out.lineages)
        .map { _meta, database, lineages -> tuple(database, lineages) }
    MMSEQS_EASY_SEARCH(
        ch_mmseqs_queries,
        ch_mmseqs_reference
    )

    /* Collect versions and build the final report. */
    ch_versions = CUSTOM_DB_FILTER.out.versions
        .mix(ITSX.out.versions)
        .mix(AMPLICON_LENGTH.out.versions)
        .mix(BUILD_DB_TAXIDS.out.versions)
        .mix(BUILD_DB_LINEAGES.out.versions)
        .mix(ADD_FLANKING.out.versions)
        .mix(FUNGAL_SCREENING.out.versions)
        .mix(ch_dereplication_versions)
        .mix(SCREEN_TAXA.out.versions)
        .mix(VSEARCH_ALLPAIRS_GLOBAL.out.versions)
        .mix(ZSCORE_SCREENING.out.versions)
        .mix(MMSEQS_EASY_SEARCH.out.versions)

    ch_collated_versions = ch_versions.collectFile(
        name: 'collated_versions.yml'
    )
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_collated_versions)

    ch_qc_files = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
    ch_multiqc_config = channel
        .fromPath(params.multiqc_config, checkIfExists: true)
        .collect()

    MULTIQC(
        ch_qc_files.collect(),
        ch_multiqc_config
    )

    emit:
    // Main databases.
    custom_db = ch_dereplicated_db
    fungi_db = ITSX.out.fungi_its2
    plants_db = ITSX.out.plants_its2
    filtered_db = CUSTOM_DB_FILTER.out.fasta
    n_filtered_db = CUSTOM_DB_FILTER.out.n_filtered

    // Sequence and taxonomy QC.
    amplicon_lengths = AMPLICON_LENGTH.out.tsv
    stats = AMPLICON_LENGTH.out.stats
    taxids = BUILD_DB_TAXIDS.out.taxids
    taxid_counts = BUILD_DB_TAXIDS.out.taxids_counts
    accession_taxid = BUILD_DB_TAXIDS.out.accession_taxid
    missing_accessions = BUILD_DB_TAXIDS.out.missing
    lineages = BUILD_DB_LINEAGES.out.lineages
    taxonomic_fasta = BUILD_DB_LINEAGES.out.fasta
    flanked_db = ADD_FLANKING.out.fasta
    flanking_missing = ADD_FLANKING.out.missing

    // Fungal screening and dereplication.
    fungal_hits = FUNGAL_SCREENING.out.hits
    fungal_like = FUNGAL_SCREENING.out.fungal_like
    no_strong_fungal_hit = FUNGAL_SCREENING.out.fasta
    dereplicated_db = ch_dereplicated_db

    // Taxon comparison and outlier screening.
    grouped_dereplicated_db = ch_family_groups
    below_minimum_taxa = SCREEN_TAXA.out.low_count
    allpairs_global = VSEARCH_ALLPAIRS_GLOBAL.out.tsv
    zscore_results = ZSCORE_SCREENING.out.scores
    zscore_outliers = ZSCORE_SCREENING.out.outliers
    zscore_non_outliers = ZSCORE_SCREENING.out.non_outliers
    mmseqs_easy_search = MMSEQS_EASY_SEARCH.out.results

    // Reports and versions.
    qc = MULTIQC.out.html
    software_versions = CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    versions = ch_versions.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.versions)
}
