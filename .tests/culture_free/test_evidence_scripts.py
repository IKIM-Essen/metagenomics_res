import importlib.util
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]


def load_script(name):
    path = ROOT / "workflow" / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.removesuffix(".py"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit_module = load_script("audit_paired_fastq.py")
pileup_module = load_script("pileup_to_allele_counts.py")


class EvidenceScriptTests(unittest.TestCase):
    def test_bundled_fastq_pair_is_synchronized(self):
        result = audit_module.audit(
            ROOT / ".tests/culture_free/reads/sample1_R1.fastq",
            ROOT / ".tests/culture_free/reads/sample1_R2.fastq",
            8,
        )
        self.assertTrue(result["synchronized"])
        self.assertEqual(result["pair_count"], 1)
        self.assertEqual(result["read_count"], 2)
        self.assertTrue(result["paired_sequence_sha256"].startswith("sha256:"))

    def test_mpileup_parser_retains_indel_and_skip_evidence(self):
        counts, deletions, skips, insertions, deletion_events = pileup_module.count_bases(
            "A", ".,Cc+2tt-1a*<>"
        )
        self.assertEqual(counts, {"A": 2, "C": 2, "G": 0, "T": 0, "N": 0})
        self.assertEqual((deletions, skips, insertions, deletion_events), (1, 2, 1, 1))

    def test_deeparg_staging_preserves_mate_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "paired.fasta"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "workflow/scripts/paired_fastq_to_fasta.py"),
                    "--r1",
                    str(ROOT / ".tests/culture_free/reads/sample1_R1.fastq"),
                    "--r2",
                    str(ROOT / ".tests/culture_free/reads/sample1_R2.fastq"),
                    "--output",
                    str(output),
                ],
                check=True,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn(">read1/1", text)
            self.assertIn(">read1/2", text)

    def test_sample_manifest_hashes_and_binds_read_identity(self):
        class NamedInputs(dict):
            __getattr__ = dict.__getitem__

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            evidence = directory / "evidence.tsv"
            evidence.write_text("evidence\n", encoding="utf-8")
            audit = directory / "audit.json"
            audit.write_text(
                json.dumps(
                    {
                        "synchronized": True,
                        "pair_count": 1,
                        "read_count": 2,
                        "r1_bases": 16,
                        "r2_bases": 16,
                        "normalized_name_sha256": "sha256:" + "a" * 64,
                        "paired_sequence_sha256": "sha256:" + "b" * 64,
                        "bottom_hash_sketch": ["0123456789abcdef"],
                    }
                ),
                encoding="utf-8",
            )
            output = directory / "manifest.json"
            fake = SimpleNamespace(
                input=NamedInputs(read_audit=str(audit), evidence=str(evidence)),
                output=SimpleNamespace(json=str(output)),
                params=SimpleNamespace(
                    sample_metadata={"sample_role": "raw_specimen"},
                    evidence_filters={"minimum_mapping_quality": 20},
                    tool_versions={"samtools": "1.23"},
                ),
                wildcards=SimpleNamespace(sample="sample1"),
            )
            runpy.run_path(
                ROOT / "workflow/scripts/lock_sample_evidence.py",
                init_globals={"snakemake": fake},
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_name"], "sample1")
            self.assertEqual(payload["paired_read_identity"]["pair_count"], 1)
            self.assertEqual(
                payload["artifacts"]["evidence"]["sha256"],
                "sha256:bdcf4c994585af6dd6cb1cfbff78bcc73ab27dc30a299db5bb83766ca05b5de4",
            )


if __name__ == "__main__":
    unittest.main()
