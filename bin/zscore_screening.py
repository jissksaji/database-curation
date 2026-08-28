#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path


SIZE_PATTERN = re.compile(r";size=(\d+);?")


def sequence_id(value):
    """Return the accession before taxonomy and VSEARCH annotations."""
    value = SIZE_PATTERN.sub("", value.strip()).rstrip(";")
    return value.split(";", 1)[0]


def read_pairs(pairwise_file, default_taxon):
    """Read three-column or taxon-prefixed VSEARCH pairwise results."""
    pairs_by_taxon = defaultdict(dict)

    with open(pairwise_file) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")

            if "query" in [field.lower() for field in fields[:3]]:
                continue

            if len(fields) >= 4:
                taxon, query, target, identity_text = fields[:4]
            elif len(fields) == 3:
                taxon = default_taxon or Path(pairwise_file).stem
                query, target, identity_text = fields
            else:
                raise ValueError(f"Malformed row {line_number}: {line.rstrip()}")

            try:
                identity = float(identity_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid identity at row {line_number}: {identity_text}"
                ) from error

            if identity < 0 or identity > 100:
                raise ValueError(
                    f"Identity outside 0-100 at row {line_number}: {identity}"
                )

            query = sequence_id(query)
            target = sequence_id(target)

            if query == target:
                continue

            pair = tuple(sorted((query, target)))
            old_identity = pairs_by_taxon[taxon].get(pair)

            if old_identity is not None and old_identity != identity:
                raise ValueError(
                    f"Conflicting duplicate pair at row {line_number}: "
                    f"{query} / {target}"
                )

            pairs_by_taxon[taxon][pair] = identity

    return pairs_by_taxon


def read_fasta(path):
    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")

            if line.startswith(">"):
                if header is not None:
                    yield header, sequence
                header = line
                sequence = []
            elif header is not None:
                sequence.append(line)

    if header is not None:
        yield header, sequence


def read_fasta_directory(fasta_directory):
    """Read the original records, keeping every header unchanged."""
    records_by_taxon = defaultdict(list)

    for fasta_file in sorted(Path(fasta_directory).glob("*.fasta")):
        taxon = fasta_file.stem

        for header, sequence in read_fasta(fasta_file):
            full_id = header[1:].split()[0]
            accession = sequence_id(full_id)
            size_match = SIZE_PATTERN.search(header)
            abundance = int(size_match.group(1)) if size_match else 1

            records_by_taxon[taxon].append(
                {
                    "raw_id": accession,
                    "accession": accession,
                    "abundance": abundance,
                    "header": header,
                    "sequence": sequence,
                }
            )

    return records_by_taxon


def calculate_scores(pairs_by_taxon, records_by_taxon, threshold, minimum):
    rows = []
    status_by_record = {}

    for taxon, records in sorted(records_by_taxon.items()):
        pairs = pairs_by_taxon.get(taxon, {})
        identifiers = [record["accession"] for record in records]
        distances = defaultdict(list)

        for (query, target), identity in pairs.items():
            distance = 1.0 - identity / 100.0
            distances[query].append(distance)
            distances[target].append(distance)

        expected_pairs = len(identifiers) * (len(identifiers) - 1) // 2
        complete = len(pairs) == expected_pairs
        can_score = len(identifiers) >= minimum and complete

        medians = {}
        for identifier in identifiers:
            if distances[identifier]:
                medians[identifier] = statistics.median(distances[identifier])

        mean_distance = statistics.mean(medians.values()) if can_score else None
        sample_sd = (
            statistics.stdev(medians.values())
            if can_score and len(medians) > 1
            else None
        )

        for record in records:
            identifier = record["accession"]
            median_distance = medians.get(identifier)
            z_score = None

            if len(identifiers) < minimum:
                status = "insufficient_sample_size"
            elif not complete:
                status = "incomplete_pair_table"
            elif sample_sd == 0:
                status = "zero_standard_deviation"
            else:
                z_score = (median_distance - mean_distance) / sample_sd
                status = "OUTLIER" if z_score > threshold else "not_outlier"

            status_by_record[(taxon, identifier)] = status
            rows.append(
                {
                    "taxon": taxon,
                    "raw_sequence_id": record["raw_id"],
                    "accession": identifier,
                    "dereplication_abundance": record["abundance"],
                    "pair_count": len(distances[identifier]),
                    "median_distance": (
                        "" if median_distance is None else f"{median_distance:.8f}"
                    ),
                    "taxon_mean_median_distance": (
                        "" if mean_distance is None else f"{mean_distance:.8f}"
                    ),
                    "taxon_sample_sd_median_distance": (
                        "" if sample_sd is None else f"{sample_sd:.8f}"
                    ),
                    "z_score": "" if z_score is None else f"{z_score:.6f}",
                    "z_threshold": threshold,
                    "status": status,
                }
            )

    rows.sort(
        key=lambda row: (
            row["taxon"],
            -(float(row["z_score"]) if row["z_score"] else float("-inf")),
        )
    )
    return rows, status_by_record


def write_table(path, rows):
    columns = [
        "taxon",
        "raw_sequence_id",
        "accession",
        "dereplication_abundance",
        "pair_count",
        "median_distance",
        "taxon_mean_median_distance",
        "taxon_sample_sd_median_distance",
        "z_score",
        "z_threshold",
        "status",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_fastas(prefix, records_by_taxon, status_by_record):
    outlier_path = prefix + ".outliers.fasta"
    non_outlier_path = prefix + ".non_outliers.fasta"

    with open(outlier_path, "w") as outliers, open(
        non_outlier_path, "w"
    ) as non_outliers:
        for taxon, records in sorted(records_by_taxon.items()):
            for record in records:
                status = status_by_record[(taxon, record["accession"])]
                output = outliers if status == "OUTLIER" else non_outliers
                output.write(record["header"] + "\n")
                output.write("\n".join(record["sequence"]) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairwise", required=True)
    parser.add_argument("--fasta-directory", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--species")
    parser.add_argument("--z-threshold", type=float, default=2.5)
    parser.add_argument("--min-sequences", type=int, default=10)
    args = parser.parse_args()

    if args.min_sequences < 3:
        parser.error("--min-sequences must be at least 3")

    pairs_by_taxon = read_pairs(args.pairwise, args.species)
    records_by_taxon = read_fasta_directory(args.fasta_directory)
    rows, statuses = calculate_scores(
        pairs_by_taxon,
        records_by_taxon,
        args.z_threshold,
        args.min_sequences,
    )

    write_table(args.output_prefix + ".zscore.tsv", rows)
    write_fastas(args.output_prefix, records_by_taxon, statuses)

    outlier_count = sum(row["status"] == "OUTLIER" for row in rows)
    print("Sequences:", len(rows))
    print("Outliers:", outlier_count)


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print("ERROR:", error)
        raise SystemExit(2)
