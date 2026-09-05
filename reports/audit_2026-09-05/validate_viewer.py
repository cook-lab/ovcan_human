"""Independently reconcile every displayed value against the processed CSVs."""
import base64
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
encoded = re.search(r'window\.OVCAN_B64="([^"\n]+)"', (ROOT / "app/data.js").read_text()).group(1)
payload = json.loads(gzip.decompress(base64.b64decode(encoded)))
assert payload["ensembl_release"] == 93

gene_symbols = {}
with (ROOT / "output/tx2gene_matched.csv").open() as handle:
    for row in csv.DictReader(handle):
        gene_symbols.setdefault(row["ensembl_gene_id"], row["external_gene_name"].strip())

expected, contributing = {}, {}
with (ROOT / "output/rna_tpm.csv").open() as handle:
    reader = csv.DictReader(handle)
    assert reader.fieldnames[1:] == payload["rna_lines"]
    for row in reader:
        gene = row[reader.fieldnames[0]]
        symbol = gene_symbols.get(gene, "")
        if not symbol or re.fullmatch(r"(nan|null|na(\.\d+)?)", symbol, re.I):
            continue
        expected.setdefault(symbol, [0.0] * len(payload["rna_lines"]))
        contributing.setdefault(symbol, []).append(gene)
        for i, model in enumerate(payload["rna_lines"]):
            expected[symbol][i] += float(row[model])
assert set(expected) == set(payload["rna"])
max_relative_error = 0.0
for symbol, values in expected.items():
    assert payload["rna_gene_id"][symbol].split(", ") == contributing[symbol]
    for original, displayed in zip(values, payload["rna"][symbol]):
        assert math.isfinite(displayed) and displayed >= 0
        error = abs(original - displayed) / max(abs(original), 1e-300)
        max_relative_error = max(error, max_relative_error)
        # Four-significant-figure display, with floating-point summation tolerance.
        assert error <= 0.000501, (symbol, original, displayed)

protein_values = 0
protein_missing = 0
with (ROOT / "output/prot_abundance_matrix.csv").open() as handle:
    reader = csv.DictReader(handle)
    assert reader.fieldnames[1:] == payload["prot_lines"]
    identifiers = set()
    for row in reader:
        feature = row[reader.fieldnames[0]]
        identifiers.add(feature)
        for i, model in enumerate(payload["prot_lines"]):
            displayed = payload["prot"][feature][i]
            if row[model] in {"", "NA", "NaN"}:
                assert displayed is None
                protein_missing += 1
            else:
                assert abs(float(row[model]) - displayed) <= 0.000501
                protein_values += 1
    assert identifiers == set(payload["prot"])

for name, digest in payload["source_sha256"].items():
    assert hashlib.sha256((ROOT / "output" / name).read_bytes()).hexdigest() == digest
assert hashlib.sha256((ROOT / "metadata/samples.csv").read_bytes()).hexdigest() == payload["metadata_sha256"]
assert len(payload["rna_lines"]) == len(payload["prot_lines"]) == 31
assert len(set(payload["rna_lines"]) & set(payload["prot_lines"])) == 30
assert set(payload["all_lines"]) == set(payload["rna_lines"]) | set(payload["prot_lines"])
assert len(payload["all_lines"]) == 32

summary = {"validation": "passed", "rna_symbols": len(expected),
           "rna_display_values": sum(len(x) for x in expected.values()),
           "rna_max_relative_rounding_error": max_relative_error,
           "protein_features": len(identifiers), "protein_display_values": protein_values,
           "protein_missing_values_preserved": protein_missing,
           "source_hashes": "all match", "display_precision": "RNA four significant figures; protein three decimal places"}
(HERE / "viewer_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
