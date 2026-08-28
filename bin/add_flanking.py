#!/usr/bin/env python3
"""Add bases around ITS2 using ITSx coordinates and the source FASTA."""

import argparse
import re


ITS2_COORDINATES = re.compile(r"(?:^|\t)ITS2: (\d+)-(\d+)(?:\t|$)")
ACCESSION_VERSION = re.compile(r"(?:\.\d+)+$")
COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def read_fasta(path):
    """Read one FASTA record at a time."""
    header = None
    sequence_lines = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)

    if header is not None:
        yield header, "".join(sequence_lines)


def record_id(header):
    """Return an accession without its version suffix."""
    accession = header.split(";", 1)[0].split()[0]
    return ACCESSION_VERSION.sub("", accession)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def write_fasta_record(handle, header, sequence):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def read_positions(path, wanted_ids):
    """Read the ITS2 coordinates for records supplied to this module."""
    positions = {}

    with open(path) as handle:
        for line in handle:
            identifier = record_id(line.split("\t", 1)[0])
            if identifier not in wanted_ids:
                continue

            match = ITS2_COORDINATES.search(line.rstrip("\r\n"))
            if match:
                positions[identifier] = (int(match.group(1)), int(match.group(2)))
            else:
                positions[identifier] = None

    return positions


def read_source_sequences(path, wanted_ids):
    """Keep only source sequences that were supplied to this module."""
    sequences = {}
    for header, sequence in read_fasta(path):
        identifier = record_id(header)
        if identifier in wanted_ids:
            sequences[identifier] = sequence
    return sequences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--its2-fasta", required=True)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--source-fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--missing", required=True)
    parser.add_argument("--forward", type=int, default=150)
    parser.add_argument("--reverse", type=int, default=150)
    args = parser.parse_args()

    if args.forward < 0 or args.reverse < 0:
        raise ValueError("Forward and reverse flank lengths must be zero or greater")

    input_records = list(read_fasta(args.its2_fasta))
    wanted_ids = {record_id(header) for header, sequence in input_records}
    positions = read_positions(args.positions, wanted_ids)
    source_sequences = read_source_sequences(args.source_fasta, wanted_ids)

    written = 0
    missing = 0

    with open(args.output, "w") as output_handle, open(
        args.missing, "w"
    ) as missing_handle:
        missing_handle.write("accession\treason\n")

        for header, its2_sequence in input_records:
            identifier = record_id(header)

            if identifier not in positions:
                missing_handle.write(
                    f"{identifier}\taccession not found in ITSx positions file\n"
                )
                missing += 1
                continue

            if positions[identifier] is None:
                missing_handle.write(
                    f"{identifier}\tnumeric ITS2 coordinates not found\n"
                )
                missing += 1
                continue

            if identifier not in source_sequences:
                missing_handle.write(
                    f"{identifier}\taccession not found in source FASTA\n"
                )
                missing += 1
                continue

            start, end = positions[identifier]
            source_sequence = source_sequences[identifier]

            if start < 1 or end < start or end > len(source_sequence):
                missing_handle.write(
                    f"{identifier}\tinvalid ITS2 coordinates {start}-{end} "
                    f"for a {len(source_sequence)} bp source sequence\n"
                )
                missing += 1
                continue

            coordinate_length = end - start + 1
            if coordinate_length != len(its2_sequence):
                missing_handle.write(
                    f"{identifier}\tITSx coordinates describe {coordinate_length} bp "
                    f"but the incoming ITS2 sequence is {len(its2_sequence)} bp\n"
                )
                missing += 1
                continue

            source_its2 = source_sequence[start - 1 : end]
            if source_its2.upper() != its2_sequence.upper():
                reversed_source = reverse_complement(source_sequence)
                reversed_its2 = reversed_source[start - 1 : end]
                if reversed_its2.upper() == its2_sequence.upper():
                    source_sequence = reversed_source
                else:
                    missing_handle.write(
                        f"{identifier}\tincoming ITS2 sequence does not match "
                        "the source sequence at the ITSx coordinates\n"
                    )
                    missing += 1
                    continue

            flanked_start = max(0, start - 1 - args.forward)
            flanked_end = min(len(source_sequence), end + args.reverse)
            flanked_sequence = source_sequence[flanked_start:flanked_end]

            write_fasta_record(output_handle, header, flanked_sequence)
            written += 1

    print(f"Flanked records: {written}")
    print(f"Missing records: {missing}")


if __name__ == "__main__":
    main()
