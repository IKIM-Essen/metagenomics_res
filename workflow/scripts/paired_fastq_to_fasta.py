#!/usr/bin/env python3
"""Stage synchronized paired FASTQ as one mate-labelled FASTA for DeepARG-SS."""

from __future__ import annotations

import argparse
import gzip


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(
        path, encoding="utf-8"
    )


def record(handle, path):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().strip()
    plus = handle.readline()
    quality = handle.readline().strip()
    if not sequence or not plus or not quality or not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"invalid or truncated FASTQ record in {path}")
    if len(sequence) != len(quality):
        raise ValueError(f"sequence/quality length mismatch in {path}")
    name = header.split()[0].removeprefix("@")
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    return name, sequence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with open_text(args.r1) as left, open_text(args.r2) as right, open(
        args.output, "w", encoding="utf-8"
    ) as output:
        while True:
            r1 = record(left, args.r1)
            r2 = record(right, args.r2)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None:
                raise ValueError("paired FASTQs contain different numbers of records")
            if r1[0] != r2[0]:
                raise ValueError(f"paired FASTQ names are not synchronized: {r1[0]!r} != {r2[0]!r}")
            output.write(f">{r1[0]}/1\n{r1[1]}\n>{r2[0]}/2\n{r2[1]}\n")


if __name__ == "__main__":
    main()
