# ITS2 database curation

Starter repository for building and maintaining a curated ITS2 reference database.

The BarBeQuE primer-benchmarking implementation has been removed. This repository
now provides a small, runnable Nextflow foundation for ITS2 curation work.

## Input

Provide one custom ITS2 reference database with a `.fasta`, `.fna`, or `.txt`
extension.

## Getting started

Run with the custom database and the installed taxonomy reference directory:

```bash
nextflow run . -profile conda \
    --custom_db custom_its2.fasta \
    --reference_base /path/to/references
```

ITSx runs the plant profiles (`T`) by default. Add fungal screening with
`--itsx T,F` (or `--itsx F`); the plant run finishes first and the fungal run
then writes separate FASTA, position, summary, and detailed-result files:

```bash
nextflow run . -profile conda \
    --custom_db custom_its2.fasta \
    --itsx T,F
```

The input is validated before curation begins and is emitted on the reference
channel with the stable ID `custom`.

To install the NCBI taxonomy references without running ITS2 curation, provide a
reference base. The reference version defaults to `1.1`:

```bash
nextflow run main.nf -profile conda \
    --reference_base /path/to/references
```

Use `--reference_version <version>` to select a different version directory.

The `BUILD_REFERENCE` workflow downloads and extracts only the NCBI new taxonomy
dump and GenBank nucleotide accession-to-taxid mapping. Both are installed under
`/path/to/references/database_curation/1.1/new_taxdump`, including
`new_taxdump/nucl_gb.accession2taxid`. Supplying only `--reference_base` installs
the references. Supplying it together with `--custom_db` runs ITS2 curation and
taxonomy annotation.

## Execution environments

Nextflow profiles are available for `local`, `conda`, `docker`, `singularity`,
`apptainer`, `podman`, `pixi`, and `slurm`. For example:

```bash
nextflow run . -profile docker
nextflow run . -profile conda
```

Pixi manages the development and launch environment. Install Pixi, then use its
repository task to run Nextflow in the declared environment:

```bash
pixi run pipeline -- --custom_db custom_its2.fasta
```

Add a module-specific `environment.yml` whenever a process needs Conda software;
the `conda` profile will activate that environment for the process. Container
profiles expect the corresponding module process to declare its container image.

The curation stages are defined in
[`workflows/its2_curation.nf`](workflows/its2_curation.nf), and the versioned
installation step is defined in
[`workflows/build_reference.nf`](workflows/build_reference.nf).

## Create the new repository

When you are ready to detach this template from BarBeQuE, remove its Git metadata,
initialize a new repository, and set your new remote:

```bash
rm -rf .git
git init
git add .
git commit -m "Initialize ITS2 database-curation repository"
git remote add origin <new-repository-url>
```
