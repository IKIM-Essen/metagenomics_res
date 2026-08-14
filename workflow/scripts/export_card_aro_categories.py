import csv
import json
import re


ARO = re.compile(r"^\d{7}$")
FIELDS = (
    "determinant_aro",
    "determinant_name",
    "category_aro",
    "category_name",
    "category_type",
    "card_version",
)

with open(snakemake.input.card_json, encoding="utf-8") as handle:
    card = json.load(handle)

version = card.get("_version")
if not isinstance(version, str) or not version:
    raise ValueError("CARD card.json is missing _version")

rows = set()
for model_id, model in card.items():
    if str(model_id).startswith("_"):
        continue
    determinant = str(model.get("ARO_accession", ""))
    determinant_name = model.get("ARO_name")
    if not ARO.fullmatch(determinant) or not isinstance(determinant_name, str):
        raise ValueError(f"CARD model {model_id!r} has an invalid ARO identity")
    categories = model.get("ARO_category", {})
    if not isinstance(categories, dict):
        raise ValueError(f"CARD model {model_id!r} has invalid ARO_category")
    for category_id, category in categories.items():
        category_aro = str(category.get("category_aro_accession", ""))
        category_name = category.get("category_aro_name")
        category_type = category.get("category_aro_class_name")
        if not ARO.fullmatch(category_aro):
            raise ValueError(
                f"CARD model {model_id!r} category {category_id!r} has invalid ARO accession"
            )
        if not isinstance(category_name, str) or not isinstance(category_type, str):
            raise ValueError(
                f"CARD model {model_id!r} category {category_id!r} has missing labels"
            )
        rows.add(
            (
                f"ARO:{determinant}",
                determinant_name,
                f"ARO:{category_aro}",
                category_name,
                category_type,
                version,
            )
        )

with open(snakemake.output.tsv, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(FIELDS)
    writer.writerows(sorted(rows))
