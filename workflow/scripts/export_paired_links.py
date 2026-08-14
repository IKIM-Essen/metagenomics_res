#!/usr/bin/env python3
"""Export one primary first-mate record per pair and aggregate contig links."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter


PAIR_FIELDS = (
    "query_name",
    "read_contig",
    "read_position_1_based",
    "mate_contig",
    "mate_position_1_based",
    "read_orientation",
    "mate_orientation",
    "mapping_quality",
    "template_length",
    "inter_contig",
)
SUMMARY_FIELDS = (
    "read_contig",
    "mate_contig",
    "read_orientation",
    "mate_orientation",
    "pair_count",
)


def open_output(path):
    return gzip.open(path, "wt", encoding="utf-8", newline="") if path.endswith(".gz") else open(
        path, "w", encoding="utf-8", newline=""
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    counts = Counter()
    with open_output(args.pairs) as pair_handle:
        writer = csv.writer(pair_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(PAIR_FIELDS)
        for line_number, line in enumerate(sys.stdin, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"SAM line {line_number} has fewer than nine columns")
            query, flag_text, contig, position, mapq, _cigar, mate_contig, mate_position, tlen = fields[:9]
            flag = int(flag_text)
            if not flag & 0x1 or not flag & 0x40:
                raise ValueError("input must contain primary first-mate paired records only")
            if mate_contig == "=":
                mate_contig = contig
            read_orientation = "R" if flag & 0x10 else "F"
            mate_orientation = "R" if flag & 0x20 else "F"
            key = (contig, mate_contig, read_orientation, mate_orientation)
            counts[key] += 1
            writer.writerow(
                (
                    query,
                    contig,
                    int(position),
                    mate_contig,
                    int(mate_position),
                    read_orientation,
                    mate_orientation,
                    int(mapq),
                    int(tlen),
                    str(contig != mate_contig).lower(),
                )
            )
    with open_output(args.summary) as summary_handle:
        writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SUMMARY_FIELDS)
        for key, count in sorted(counts.items()):
            writer.writerow((*key, count))


if __name__ == "__main__":
    main()
