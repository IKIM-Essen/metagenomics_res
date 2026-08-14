"""Validate and export the culture-free sample/control design from PEP."""

import csv
import json
from pathlib import Path


FIELDS = (
    "sample_name",
    "sample_role",
    "case_id",
    "material_stage",
    "matched_control",
    "expected_taxon",
)
ROLES = {
    "raw_specimen",
    "plate_sweep",
    "picked_isolate",
    "negative_control",
    "positive_control",
    "unspecified",
}
STAGES = {"raw_specimen", "plate_sweep", "picked_isolates", "control", "unspecified"}

records = [dict(item) for item in snakemake.params.records]
if not records:
    raise ValueError("culture-free sample design contains no records")
by_name = {item["sample_name"]: item for item in records}
if len(by_name) != len(records):
    raise ValueError("culture-free sample design contains duplicate sample_name values")

for item in records:
    missing = [field for field in FIELDS[:4] if not str(item.get(field, "")).strip()]
    if missing:
        raise ValueError(f"sample {item.get('sample_name')!r} is missing {', '.join(missing)}")
    if item["sample_role"] not in ROLES or (
        snakemake.params.enforce_complete and item["sample_role"] == "unspecified"
    ):
        raise ValueError(f"sample {item['sample_name']!r} has invalid sample_role")
    if item["material_stage"] not in STAGES or (
        snakemake.params.enforce_complete and item["material_stage"] == "unspecified"
    ):
        raise ValueError(f"sample {item['sample_name']!r} has invalid material_stage")
    control = str(item.get("matched_control", "")).strip()
    if snakemake.params.require_matched_negative and item["sample_role"] in {
        "raw_specimen",
        "plate_sweep",
    }:
        if not control:
            raise ValueError(f"sample {item['sample_name']!r} lacks a matched negative control")
        if control not in by_name or by_name[control]["sample_role"] != "negative_control":
            raise ValueError(
                f"sample {item['sample_name']!r} matched_control must name a negative_control"
            )

with open(snakemake.output.tsv, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for item in sorted(records, key=lambda row: row["sample_name"]):
        writer.writerow({field: item.get(field, "") for field in FIELDS})

payload = {
    "schema": "resmag-culture-free-sample-design-v1",
    "require_matched_negative": bool(snakemake.params.require_matched_negative),
    "complete_metadata_enforced": bool(snakemake.params.enforce_complete),
    "sample_count": len(records),
    "role_counts": {
        role: sum(item["sample_role"] == role for item in records) for role in sorted(ROLES)
    },
    "valid": True,
}
Path(snakemake.output.json).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
