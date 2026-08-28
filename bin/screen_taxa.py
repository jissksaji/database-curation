#!/usr/bin/env python3

import argparse
import re
from collections import Counter
from pathlib import Path


def read_fasta(fasta_path):
    header = None
    sequence = []

    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    yield header, sequence

                header = line
                sequence = []
            else:
                sequence.append(line)

        if header is not None:
            yield header, sequence


def get_rank(header, rank):
    marker = ";" + rank + ":"

    if marker not in header:
        return ""

    value = header.split(marker, 1)[1]
    value = value.split(";", 1)[0]
    value = value.split("\t", 1)[0]
    return value.strip()


def get_taxon(header):
    species = get_rank(header, "s")
    genus = get_rank(header, "g")

    if species:
        return "species", species
    if genus:
        return "genus", genus
    return "", ""


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def write_record(output_file, header, sequence):
    output_file.write(header + "\n")
    output_file.write("\n".join(sequence) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--minimum", type=int, default=15)
    parser.add_argument("--low-output", required=True)
    parser.add_argument("--family-dir", required=True)
    args = parser.parse_args()

    counts = Counter()

    for header, sequence in read_fasta(args.input):
        rank, taxon = get_taxon(header)
        if taxon:
            counts[(rank, taxon)] += 1

    family_dir = Path(args.family_dir)
    family_dir.mkdir(parents=True, exist_ok=True)
    family_dirs = set()
    low_records = 0

    with open(args.low_output, "w") as low_output:
        for header, sequence in read_fasta(args.input):
            rank, taxon = get_taxon(header)
            group = (rank, taxon)

            if not taxon or counts[group] < args.minimum:
                write_record(low_output, header, sequence)
                low_records += 1
                continue

            family = get_rank(header, "f") or "unclassified_family"
            output_dir = family_dir / safe_name(family)
            output_dir.mkdir(parents=True, exist_ok=True)
            taxon_path = output_dir / (safe_name(taxon) + ".fasta")

            with open(taxon_path, "a") as taxon_output:
                write_record(taxon_output, header, sequence)

            family_dirs.add(output_dir)

    print("Taxon groups:", len(counts))
    print("Records below minimum:", low_records)
    print("Family groups for VSEARCH:", len(family_dirs))


if __name__ == "__main__":
    main()
