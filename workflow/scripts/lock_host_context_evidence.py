"""Bind optional ARG-host context to one immutable ResMAG sample manifest."""

import hashlib
import json
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


with open(snakemake.input.sample_manifest, encoding="utf-8") as handle:
    sample_manifest = json.load(handle)
if sample_manifest.get("schema") != "resmag-culture-free-sample-evidence-lock-v1":
    raise ValueError("host context requires a ResMAG sample evidence lock v1")
if sample_manifest.get("sample_name") != snakemake.wildcards.sample:
    raise ValueError("host context sample does not match the parent sample manifest")

paths = {name: str(path) for name, path in snakemake.input.items()}
artifacts = {
    name: {
        "path": path,
        "size_bytes": Path(path).stat().st_size,
        "sha256": sha256_file(path),
    }
    for name, path in sorted(paths.items())
}

parent_bound = ("assembly_unicard", "coverage", "gfa", "paired_link_summary", "paired_links")
parent_artifacts = sample_manifest.get("artifacts", {})
for name in parent_bound:
    parent = parent_artifacts.get(name)
    if not isinstance(parent, dict) or parent.get("sha256") != artifacts[name]["sha256"]:
        raise ValueError(f"host context artifact {name} does not match the parent sample manifest")

payload = {
    "schema": "resmag-culture-free-host-context-lock-v1",
    "sample_name": snakemake.wildcards.sample,
    "parent_sample_evidence": {
        "path": str(snakemake.input.sample_manifest),
        "sha256": artifacts["sample_manifest"]["sha256"],
        "schema": sample_manifest["schema"],
    },
    "tool_versions": dict(snakemake.params.tool_versions),
    "evidence_roles": {
        "arg_locus_sources": ["annotations", "assembly_unicard", "proteins"],
        "candidate_host_sources": ["contig_to_bin", "bin_summary", "bin_taxonomy"],
        "mobile_element_sources": ["plasmid_summary"],
        "physical_link_sources": [
            "coverage",
            "gfa",
            "paired_link_summary",
            "paired_links",
        ],
    },
    "artifacts": artifacts,
    "claim_boundary": [
        "bin and plasmid annotations are candidate host context, not strain ownership",
        "physical links require calibrated scoring and conflict checks downstream",
        "an unlinked ARG must not be assigned to every organism",
        "host context cannot establish susceptibility or MIC",
    ],
}
Path(snakemake.output.json).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
