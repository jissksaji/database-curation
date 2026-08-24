workflow INPUT_CHECK {
    take:
    custom_db

    main:
    custom_db.map { meta, database_file ->
        def file_name = database_file.getName().toLowerCase()
        def valid_file = file_name.endsWith('.fasta') ||
                         file_name.endsWith('.fna') ||
                         file_name.endsWith('.txt')

        if (!valid_file) {
            throw new IllegalArgumentException(
                "--custom_db must be a .fasta, .fna, or .txt file: ${database_file}"
            )
        }

        tuple(meta, database_file)
    }.set { validated_db }

    emit:
    db = validated_db
}
