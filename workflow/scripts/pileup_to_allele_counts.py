#!/usr/bin/env python3
"""Convert filtered samtools mpileup text to explicit nucleotide counts."""

from __future__ import annotations

import csv
import sys


FIELDS = (
    "contig",
    "position_1_based",
    "reference_base",
    "mpileup_depth",
    "A",
    "C",
    "G",
    "T",
    "N",
    "deletion_placeholders",
    "reference_skips",
    "insertion_events",
    "deletion_events",
    "usable_base_depth",
    "minor_allele_fraction",
)


def count_bases(reference, bases):
    counts = {base: 0 for base in "ACGTN"}
    deletion_placeholders = 0
    reference_skips = 0
    insertion_events = 0
    deletion_events = 0
    index = 0
    while index < len(bases):
        symbol = bases[index]
        if symbol == "^":
            index += 2
            continue
        if symbol == "$":
            index += 1
            continue
        if symbol in "+-":
            is_insertion = symbol == "+"
            index += 1
            start = index
            while index < len(bases) and bases[index].isdigit():
                index += 1
            if start == index:
                raise ValueError("indel marker lacks a length")
            length = int(bases[start:index])
            index += length
            if is_insertion:
                insertion_events += 1
            else:
                deletion_events += 1
            continue
        if symbol in ".,":
            counts[reference] = counts.get(reference, 0) + 1
        elif symbol.upper() in counts:
            counts[symbol.upper()] += 1
        elif symbol in "*#":
            deletion_placeholders += 1
        elif symbol in "<>":
            reference_skips += 1
        else:
            raise ValueError(f"unsupported mpileup base symbol {symbol!r}")
        index += 1
    return counts, deletion_placeholders, reference_skips, insertion_events, deletion_events


def main():
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(FIELDS)
    for line_number, line in enumerate(sys.stdin, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            raise ValueError(f"mpileup line {line_number} has fewer than five columns")
        contig, position, reference, depth_text, bases = fields[:5]
        reference = reference.upper()
        if reference not in "ACGTN":
            reference = "N"
        counts, deletions, skips, insertions, deletion_events = count_bases(reference, bases)
        usable = sum(counts.values())
        ordered = sorted((counts[base] for base in "ACGT"), reverse=True)
        minor_fraction = ordered[1] / usable if usable and len(ordered) > 1 else 0.0
        writer.writerow(
            (
                contig,
                int(position),
                reference,
                int(depth_text),
                *(counts[base] for base in "ACGTN"),
                deletions,
                skips,
                insertions,
                deletion_events,
                usable,
                f"{minor_fraction:.12g}",
            )
        )


if __name__ == "__main__":
    main()
