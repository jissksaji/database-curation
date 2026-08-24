include { GUNZIP as GUNZIP_TAXDUMP } from '../modules/helper/gunzip/main'
include { GUNZIP as GUNZIP_ACCESSION2TAXID } from '../modules/helper/gunzip/main'

workflow BUILD_REFERENCE {
    main:
    if (!params.reference_base) {
        throw new IllegalArgumentException(
            'Please provide an installation directory with --reference_base.'
        )
    }

    taxdump = file(params.references.taxdump_url, checkIfExists: true)
    accession2taxid = file(
        params.references.accession2taxid_url,
        checkIfExists: true
    )

    ch_taxdump = channel.of(
        tuple([id: 'taxdump'], taxdump)
    )
    ch_accession2taxid = channel.of(
        tuple([id: 'accession2taxid'], accession2taxid)
    )

    GUNZIP_TAXDUMP(ch_taxdump)
    GUNZIP_ACCESSION2TAXID(ch_accession2taxid)
}
