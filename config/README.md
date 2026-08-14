# Configuration

`config.yaml` selects the PEP project, databases, QC thresholds, and optional culture-free evidence
profile. `pep/config.yaml` points to the sample table.

The basic PEP columns are `sample_name`, `fq1`, and `fq2`. If
`culture-free-evidence.control-aware.enabled` is true, also provide:

| Column | Values and meaning |
|---|---|
| `sample_role` | `raw_specimen`, `plate_sweep`, `picked_isolate`, `negative_control`, or `positive_control` |
| `case_id` | Deidentified patient episode or control-batch identity |
| `material_stage` | `raw_specimen`, `plate_sweep`, `picked_isolates`, or `control` |
| `matched_control` | Sample name of the matched negative control for every raw specimen and plate sweep |
| `expected_taxon` | Optional controlled taxon such as `NCBITaxon:573` for a positive control |

The culture-free profile is fail-closed where provenance affects interpretation:

- set an exact UniCARD builder commit and UniRef release;
- choose explicit mapping/base-quality and pileup-depth limits;
- pin DeepARG's container, runtime/model release, and database digest before enabling it;
- retain the generated per-sample evidence manifest with downstream results.

The supplied paths and `REPLACE_...` values are templates. They are not usable database identities.
