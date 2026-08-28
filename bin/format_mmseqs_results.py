#!/usr/bin/env python3

import argparse
import csv
import re


VERSION_PATTERN = re.compile(r"\.\d+$")
RANK_CODES = {
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
RANKS = tuple(RANK_CODES.values())
ALIGNMENT_COLUMNS = (
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "qstart",
    "qend",
    "qlen",
    "tstart",
    "tend",
    "tlen",
    "evalue",
    "bits",
)


def accession(sequence_id):
    """Return the accession before taxonomy, size, or version annotations."""
    value = sequence_id.strip().lstrip(">").split(";", 1)[0].split()[0]
    return VERSION_PATTERN.sub("", value)


def normalise_taxon(value):
    return "_".join(value.strip().split())


def parse_taxonomy(value):
    taxonomy = {rank: "" for rank in RANKS}
    for field in value.rstrip("\r\n").split(";"):
        if ":" not in field:
            continue
        code, name = field.split(":", 1)
        rank = RANK_CODES.get(code)
        if rank:
            taxonomy[rank] = normalise_taxon(name)
    return taxonomy


def collect_hit_accessions(input_path):
    wanted = set()
    with open(input_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query", "target", *ALIGNMENT_COLUMNS}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "MMseqs table is missing columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            wanted.add(accession(row["query"]))
            wanted.add(accession(row["target"]))
    return wanted


def read_query_taxonomy(fasta_path, wanted):
    taxonomy_by_accession = {}
    with open(fasta_path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line.rstrip("\r\n")
            identifier = accession(header)
            if identifier in wanted:
                taxonomy_by_accession[identifier] = parse_taxonomy(header)
    return taxonomy_by_accession


def read_lineages(lineage_path, wanted):
    """Scan the lineage table once and retain only accessions in MMseqs hits."""
    taxonomy_by_accession = {}
    remaining = set(wanted)
    with open(lineage_path) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t", 2)
            if len(fields) < 3:
                continue
            identifier = accession(fields[0])
            if identifier not in remaining:
                continue
            taxonomy_by_accession[identifier] = parse_taxonomy(fields[2])
            remaining.remove(identifier)
            if not remaining:
                break
    return taxonomy_by_accession


def name_for(identifier, taxonomy):
    for rank in ("species", "genus", "family", "order", "class"):
        if taxonomy[rank]:
            return taxonomy[rank]
    return identifier


def add_sequence_columns(output_row, prefix, identifier, taxonomy):
    output_row[f"{prefix}_accession"] = identifier
    output_row[f"{prefix}_name"] = name_for(identifier, taxonomy)
    for rank in RANKS:
        output_row[f"{prefix}_{rank}"] = taxonomy[rank]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--lineages", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    wanted = collect_hit_accessions(args.input)
    taxonomy_by_accession = read_lineages(args.lineages, wanted)
    query_taxonomy = read_query_taxonomy(args.query_fasta, wanted)
    for identifier, taxonomy in query_taxonomy.items():
        taxonomy_by_accession.setdefault(identifier, taxonomy)

    output_columns = [
        "query_group",
        "query_accession",
        "query_name",
        "query_class",
        "query_order",
        "query_family",
        "query_genus",
        "query_species",
        "target_accession",
        "target_name",
        "target_class",
        "target_order",
        "target_family",
        "target_genus",
        "target_species",
        *ALIGNMENT_COLUMNS,
    ]

    with open(args.input) as input_handle, open(
        args.output, "w", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=output_columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            query_accession = accession(row["query"])
            target_accession = accession(row["target"])
            empty_taxonomy = {rank: "" for rank in RANKS}
            query_taxonomy = taxonomy_by_accession.get(
                query_accession, empty_taxonomy
            )
            target_taxonomy = taxonomy_by_accession.get(
                target_accession, empty_taxonomy
            )

            output_row = {
                "query_group": query_taxonomy["family"]
                or "unclassified_family"
            }
            add_sequence_columns(
                output_row, "query", query_accession, query_taxonomy
            )
            add_sequence_columns(
                output_row, "target", target_accession, target_taxonomy
            )
            for column in ALIGNMENT_COLUMNS:
                output_row[column] = row[column]
            writer.writerow(output_row)


if __name__ == "__main__":
    main()
