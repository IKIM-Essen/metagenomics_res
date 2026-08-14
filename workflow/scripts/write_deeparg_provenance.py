"""Lock configured DeepARG runtime/model/database identities."""

import json
import re
from pathlib import Path


digest = snakemake.params.database_digest
if not re.fullmatch(r"sha256:[0-9a-f]{64}", digest):
    raise ValueError("DeepARG database-digest must be an exact sha256 digest")
container = snakemake.params.container
if "@sha256:" not in container and not re.search(r"\.sif$", container):
    raise ValueError("DeepARG container must be a digest-pinned URI or an explicit .sif path")
if not snakemake.params.model_release:
    raise ValueError("DeepARG model-release must not be blank")

payload = {
    "schema": "resmag-deeparg-evidence-lock-v1",
    "runtime_version": snakemake.params.runtime_version,
    "model": "SS",
    "model_release": snakemake.params.model_release,
    "database_path": snakemake.params.database_path,
    "database_digest": digest,
    "container": container,
    "parameters": dict(snakemake.params.parameters),
    "claim_boundary": [
        "DeepARG is secondary broad-class evidence",
        "DeepARG output does not establish an exact allele, host, regulatory mechanism, or MIC",
        "potential.ARG rows are retained and must not be silently promoted",
    ],
}
Path(snakemake.output.json).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
