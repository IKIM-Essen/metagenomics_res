"""Create a content-addressed per-sample evidence manifest."""

import hashlib
import json
import subprocess
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def flatten_inputs():
    flattened = {}
    for name, value in snakemake.input.items():
        paths = list(value) if isinstance(value, (list, tuple)) else [value]
        for index, path in enumerate(paths):
            key = name if len(paths) == 1 else f"{name}_{index + 1}"
            flattened[key] = str(path)
    return flattened


def git_value(*arguments):
    try:
        return subprocess.run(
            ["git", *arguments],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unavailable"


paths = flatten_inputs()
artifacts = {
    name: {
        "path": path,
        "size_bytes": Path(path).stat().st_size,
        "sha256": sha256_file(path),
    }
    for name, path in sorted(paths.items())
}
with open(snakemake.input.read_audit, encoding="utf-8") as handle:
    read_audit = json.load(handle)
if not read_audit.get("synchronized") or read_audit.get("pair_count", 0) <= 0:
    raise ValueError("paired-read audit is not a passing non-empty synchronized pair set")

revision = git_value("rev-parse", "HEAD")
dirty = bool(git_value("status", "--porcelain") not in ("", "unavailable"))
payload = {
    "schema": "resmag-culture-free-sample-evidence-lock-v1",
    "sample_name": snakemake.wildcards.sample,
    "sample_metadata": dict(snakemake.params.sample_metadata),
    "resmag_revision": revision,
    "resmag_worktree_dirty_at_execution": dirty,
    "evidence_filters": dict(snakemake.params.evidence_filters),
    "tool_versions": dict(snakemake.params.tool_versions),
    "paired_read_identity": {
        key: read_audit[key]
        for key in (
            "pair_count",
            "read_count",
            "r1_bases",
            "r2_bases",
            "normalized_name_sha256",
            "paired_sequence_sha256",
            "bottom_hash_sketch",
        )
    },
    "artifacts": artifacts,
    "claim_boundary": [
        "outputs are sequence and physical-link evidence, not susceptibility calls",
        "unlinked AMR evidence must not be assigned to every organism",
        "absence of evidence cannot establish susceptibility",
    ],
}
Path(snakemake.output.json).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
