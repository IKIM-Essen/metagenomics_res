import hashlib
import json
import re
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


if not re.fullmatch(r"[0-9a-fA-F]{40}", snakemake.params.builder_commit):
    raise ValueError("UniCARD builder-commit must be an exact 40-character Git commit")
if not snakemake.params.uniref_release or "REPLACE" in snakemake.params.uniref_release:
    raise ValueError("uniref-release must identify the exact UniRef database release")
with open(snakemake.input.card_json, encoding="utf-8") as handle:
    card_version = json.load(handle).get("_version")
if not isinstance(card_version, str) or not card_version:
    raise ValueError("CARD card.json is missing _version")
if snakemake.params.card_version.removeprefix("v") != card_version.removeprefix("v"):
    raise ValueError("configured CARD version does not match card.json")

artifacts = {
    name: {"path": str(path), "sha256": sha256_file(path)}
    for name, path in {
        "card_json": snakemake.input.card_json,
        "unicard_fasta": snakemake.input.unicard_fasta,
        "unicard_dmnd": snakemake.input.unicard_dmnd,
        "legacy_unicard_hierarchy": snakemake.input.hierarchy,
        "structured_aro_categories": snakemake.input.categories,
    }.items()
}
payload = {
    "schema": "resmag-culture-free-amr-evidence-lock-v2",
    "card_version": snakemake.params.card_version,
    "unicard_builder_commit": snakemake.params.builder_commit,
    "uniref_release": snakemake.params.uniref_release,
    "direct_search": {
        "diamond_version": snakemake.params.diamond_version,
        "mode": "blastx",
        "sensitivity": snakemake.params.sensitivity,
        "evalue": snakemake.params.evalue,
        "max_target_seqs": snakemake.params.max_target_seqs,
        "outfmt_fields": [
            "qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue",
            "bitscore", "qlen", "slen", "qcovhsp", "scovhsp", "qframe"
        ],
    },
    "assembly_protein_search": {
        "diamond_version": snakemake.params.diamond_version,
        "mode": "blastp",
        "sensitivity": snakemake.params.sensitivity,
        "evalue": snakemake.params.evalue,
        "max_target_seqs": snakemake.params.max_target_seqs,
        "outfmt_fields": [
            "qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue",
            "bitscore", "qlen", "slen", "qcovhsp", "scovhsp"
        ],
    },
    "artifacts": artifacts,
    "warnings": [
        "DIAMOND-UniCARD hits are candidate sequence evidence, not phenotype calls",
        "direct-read and assembly-protein searches are separate evidence views",
        "legacy UniCARD hierarchy category labels are display-only; use structured_aro_categories",
        "absence of a hit cannot establish susceptibility",
    ],
}
Path(snakemake.output.json).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
