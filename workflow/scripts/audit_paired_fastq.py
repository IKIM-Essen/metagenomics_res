#!/usr/bin/env python3
"""Stream a paired FASTQ once and emit synchronization and identity evidence."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import json
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(
        path, encoding="utf-8"
    )


def normalized_name(header):
    name = header.strip().split()[0].removeprefix("@")
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    return name


def read_record(handle, path):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\n\r")
    plus = handle.readline()
    quality = handle.readline().rstrip("\n\r")
    if not sequence or not plus or not quality:
        raise ValueError(f"truncated FASTQ record in {path}")
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"invalid FASTQ structure in {path}")
    if len(sequence) != len(quality):
        raise ValueError(f"sequence/quality length mismatch in {path}")
    return header, sequence, quality


def audit(r1_path, r2_path, sketch_size):
    name_digest = hashlib.sha256()
    sequence_digest = hashlib.sha256()
    sketch = []
    pair_count = 0
    r1_bases = 0
    r2_bases = 0

    with open_text(r1_path) as r1, open_text(r2_path) as r2:
        while True:
            left = read_record(r1, r1_path)
            right = read_record(r2, r2_path)
            if left is None and right is None:
                break
            if left is None or right is None:
                raise ValueError("paired FASTQs contain different numbers of records")
            left_name = normalized_name(left[0])
            right_name = normalized_name(right[0])
            if left_name != right_name:
                raise ValueError(
                    f"paired FASTQ names are not synchronized: {left_name!r} != {right_name!r}"
                )
            name_digest.update(left_name.encode("utf-8") + b"\n")
            pair_bytes = (left[1] + "\0" + right[1]).encode("ascii")
            sequence_digest.update(pair_bytes + b"\n")
            fingerprint = int.from_bytes(hashlib.sha256(pair_bytes).digest()[:8], "big")
            if len(sketch) < sketch_size:
                heapq.heappush(sketch, -fingerprint)
            elif fingerprint < -sketch[0]:
                heapq.heapreplace(sketch, -fingerprint)
            pair_count += 1
            r1_bases += len(left[1])
            r2_bases += len(right[1])

    return {
        "schema": "resmag-paired-fastq-audit-v1",
        "r1": str(r1_path),
        "r2": str(r2_path),
        "synchronized": True,
        "pair_count": pair_count,
        "read_count": pair_count * 2,
        "r1_bases": r1_bases,
        "r2_bases": r2_bases,
        "normalized_name_sha256": "sha256:" + name_digest.hexdigest(),
        "paired_sequence_sha256": "sha256:" + sequence_digest.hexdigest(),
        "bottom_hash_sketch": [f"{value:016x}" for value in sorted(-item for item in sketch)],
        "sketch_size": sketch_size,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sketch-size", type=int, default=128)
    args = parser.parse_args()
    if args.sketch_size <= 0:
        raise ValueError("sketch-size must be positive")
    payload = audit(args.r1, args.r2, args.sketch_size)
    Path(args.output).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
