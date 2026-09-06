#!/usr/bin/env python3
"""Build the processed-data release from current, validated analysis outputs.

Uses only Python's standard library. Run after run_all.sh (or after all affected
analysis scripts). The archival source inventory is never edited by this builder.
"""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "release"
OUT = ROOT / "output"


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def truth(value):
    return str(value).upper() in {"TRUE", "Y", "1"}


HISTOTYPES = {
    "HGS": "High-grade serous carcinoma", "LGS": "Low-grade serous carcinoma",
    "CC": "Clear cell carcinoma", "EC": "Endometrioid carcinoma (historical grouping)",
    "MC": "Mucinous carcinoma", "MMMT": "Carcinosarcoma",
    "SCCOHT": "Small cell carcinoma of the ovary, hypercalcaemic type",
}
SITES = {"Mes-Masson": "CHUM", "Huntsman": "BC Cancer-OVCARE",
         "Huntsman/Vanderhyden": "OHRI"}


def main():
    master = read(ROOT / "metadata/samples.csv")
    master = [r for r in master if r["provenance"] == "generated" and r["analysis_include"] == "Y"]
    families = {r["cell_line"]: r for r in read(ROOT / "metadata/line_family_map.csv")}
    quality = {r["cell_line"]: r for r in read(OUT / "supplement_per_line.csv")}
    assert len(master) == 42 and len({r["cell_line"] for r in master}) == 42
    assert set(families) == {r["cell_line"] for r in master}
    models = []
    for row in master:
        f = families[row["cell_line"]]
        m = {"cell_line": row["cell_line"], "histotype_code": row["subtype"],
             "histotype_label": HISTOTYPES[row["subtype"]],
             "histotype_note": ("Published reassignment to dedifferentiated ovarian carcinoma; retained in the historical EC group for analysis."
                                if row["cell_line"] == "TOV112D" else row["subtype_status"]),
             "contributing_centre": SITES[row["source_site"]], "patient_id": f["patient_id"],
             "models_from_same_patient": f["n_lines_in_family"],
             "selected_patient_model": f["patient_representative"],
             "rna_available": f["has_rna"], "protein_available": f["has_prot"],
             "copy_number_available": f["has_wes_cnv"], "variants_available": f["has_wes_maf"],
             "rna_sample_id": row["rna_sample_id"], "rna_passage": row["rna_passage"],
             "wes_passage": row["wes_passage"], "tmt_plex": row["tmt_plex"],
             "tmt_channel": row["tmt_channel"],
             "stock_str_documentation": "Awaiting contributing-laboratory confirmation",
             "stock_mycoplasma_documentation": "Awaiting contributing-laboratory confirmation"}
        for field in ("rna_pseudoalign_pct", "rna_sequenced_fragments_M", "rna_assigned_gene_counts_M",
                      "rna_genes_detected", "prot_n_detected", "prot_pct_missing", "wes_n_coding",
                      "fga_autosome", "cellosaurus_accession", "str_profile_documented", "problematic_flag"):
            m[field] = quality[row["cell_line"]][field]
        for key, value in m.items():
            if value in {"-", "NA", "N/A"}:
                m[key] = ""
        models.append(m)
    assert len({r["patient_id"] for r in models}) == 34
    assert sum(truth(r["selected_patient_model"]) for r in models) == 34
    write(DEST / "metadata/models.csv", models)
    # Same cohort-only table is discoverable beside the internal source inventory.
    write(ROOT / "metadata/resource_models.csv", models)

    records = []
    def add(source, relpath, unit, description):
        assert source.is_file(), f"Missing release input: {source}"
        target = DEST / relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        records.append({"file": relpath, "source": str(source.relative_to(ROOT)),
                        "unit": unit, "description": description})

    records.append({"file": "metadata/models.csv", "source": "metadata/samples.csv; metadata/line_family_map.csv",
                    "unit": "model", "description": "The complete 42-model resource; no cohort filter is required."})
    for name, unit, desc in [
        ("rna_tpm.csv", "gene by model", "Gene TPM from the matched transcript-to-gene reference; unfiltered genes."),
        ("rna_counts.csv", "gene by model", "Filtered tximport estimated counts; not raw sequenced-fragment counts."),
        ("rna_vst.rds", "gene by model", "DESeq2 variance-stabilised expression matrix, in R format."),
        ("rna_qc_metrics.csv", "RNA model", "Pseudoalignment, processed fragments, assigned gene counts, and gene detection."),
        ("prot_abundance_matrix.csv", "protein feature by model", "Supplied log2 protein abundance normalised using the pooled TMT standard; missing measurements retained."),
        ("prot_qc.csv", "protein feature", "Row-matched protein identity, search support, and measurement-presence annotations."),
        ("prot_block_missingness.csv", "protein feature", "Feature detection across models and TMT plexes."),
        ("prot_sample_qc.csv", "protein model", "Per-model protein detection and missingness."),
        ("prot_feature_accounting.csv", "feature-set definition", "Counts and definitions of each protein feature set."),
        ("wes_mutations_filtered.csv", "candidate variant by model", "PASS, rare, coding tumour-only candidates; somatic origin unconfirmed."),
        ("wes_driver_tiers.csv", "prioritised candidate by model", "Heuristic evidence tiers for candidate alterations in selected genes."),
        ("wes_cnv_segments.csv", "copy-number segment by model", "Relative total-copy segments; genomic intervals are zero-based, half-open."),
        ("wes_cnv_fga.csv", "copy-number model", "Fraction of assessed segment length altered, with autosomal estimates preferred."),
        ("wes_cnv_arm_freq_patient.csv", "chromosome arm", "HGSC arm-call frequencies; denominator and patient aggregation recorded per row."),
    ]:
        add(OUT / name, "processed/" + name, unit, desc)
    # Recoverable reference records have stable names and are included when present.
    for name in ("tx2gene_matched.csv", "rna_gene_annotation.csv", "rna_reference_reconciliation.csv",
                 "wes_input_manifest.csv", "wes_filter_cascade.csv",
                 "wes_cnv_arm_boundaries.csv", "wes_vcf_header_provenance.csv",
                 "wes_cnvkit_target_intervals.bed", "wes_target_trinucleotide_opportunities.csv",
                 "wes_cosmic_target_normalized.csv", "external_reference_provenance.csv",
                 "external_reference_provenance.json"):
        if (OUT / name).is_file():
            add(OUT / name, "reference/" + name, "reference/provenance record", "Exact reference or input-accounting record produced by the analysis.")

    add(ROOT / "data/reference/depmap_24Q4_figshare_v1_metadata.json",
        "reference/depmap_24Q4_figshare_v1_metadata.json", "external reference metadata",
        "Official Figshare version-1 metadata for DepMap 24Q4; source file MD5 and sizes match the pinned local inputs.")
    for record in read(OUT / "external_reference_provenance.csv"):
        source = ROOT / record["file"]
        assert sha256(source) == record["sha256"], "External reference manifest is stale: " + record["file"]

    for name, desc in [
        ("rna_reference_sensitivity.csv", "Before/after reference and locus-selection sensitivity."),
        ("integ_rnaprot_cor.csv", "Model-level descriptive RNA–protein concordance."),
        ("integ_rnaprot_patientrep_cor.csv", "Primary gene-wise associations using one paired model per patient."),
        ("integ_rnaprot_patientrep_sensitivity.csv", "Patient-selection concordance sensitivity."),
        ("rna_marker_effectsizes.csv", "Patient-based marker effects; exploratory and partly centre-confounded."),
        ("rna_marker_recovery_permutation.csv", "Joint marker-ordering permutation checks; centre-restricted result is primary."),
        ("sensitivity_patient_reps_pca.csv", "Comparable model-level and patient-based PCA/commonality decompositions."),
        ("sensitivity_patient_reps_pc_permutation.csv", "Patient-based reduced-model PC permutation checks."),
        ("prot_bridge_cor.csv", "Observed agreement for the four repeated cross-plex bridge lines."),
        ("prot_bridge_agreement.csv", "Bland–Altman agreement and conditional repeatability estimates."),
        ("wes_cnv_arm_calls.csv", "Per-model arm calls, overlaps, and full-arm denominators."),
        ("external_ccle_concordance.csv", "External expression-reference identity consistency checks."),
    ]:
        add(OUT / name, "validation/" + name, "technical validation or sensitivity", desc)
    for name in ("wes_sbs_context_target_restricted.csv", "wes_signature_target_sensitivity.csv"):
        add(OUT / name, "exploratory/" + name, "exploratory mutational-signature sensitivity",
            "Tumour-only candidate signatures and target-opportunity sensitivity; not an MSI diagnosis or clinical classifier.")

    add(ROOT / "data/reference/rna_reference_provenance.json", "reference/rna_reference_provenance.json",
        "reference provenance", "Official Ensembl reference URLs and hashes; target-set verification and mapping policy.")
    # The public-facing provenance uses the packaged reference, not a local recovery location.
    provenance_path = DEST / "reference/rna_reference_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    provenance.pop("recovered_from", None)
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    import base64, gzip, re
    payload = json.loads(gzip.decompress(base64.b64decode(
        re.search(r'window\.OVCAN_B64="([^"\n]+)"', (ROOT / "app/data.js").read_text()).group(1))))
    assert payload["ensembl_release"] == 93
    for name, digest in payload["source_sha256"].items():
        assert sha256(OUT / name) == digest, "Viewer is stale; rerun app/build_payload.py and app/build_single.py"
    assert payload["metadata_sha256"] == sha256(ROOT / "metadata/samples.csv"), "Viewer metadata are stale"
    expected_viewer = (ROOT / "app/index.html").read_text().replace(
        '<script src="data.js"></script>', '<script>\n' + (ROOT / "app/data.js").read_text() + '\n</script>')
    assert (ROOT / "app/ovcan_viewer_standalone.html").read_text() == expected_viewer, "Standalone viewer is stale; rerun app/build_single.py"
    add(ROOT / "app/ovcan_viewer_standalone.html", "viewer/index.html", "gene-symbol lookup",
        "Local standalone viewer rebuilt against the same processed matrices; public deployment is separate.")

    # Ensure availability refers to actual current matrices, not inherited flags.
    roster = {r["cell_line"] for r in models}
    for name, flag in [("rna_tpm.csv", "rna_available"), ("rna_counts.csv", "rna_available"),
                       ("prot_abundance_matrix.csv", "protein_available")]:
        with (OUT / name).open(newline="") as handle:
            columns = next(csv.reader(handle))[1:]
        assert set(columns) == {r["cell_line"] for r in models if truth(r[flag])}, name
    for name, flag in [("wes_mutations_filtered.csv", "variants_available"), ("wes_cnv_segments.csv", "copy_number_available")]:
        assert {r["cell_line"] for r in read(OUT / name)} == {r["cell_line"] for r in models if truth(r[flag])}, name
    prot = read(OUT / "prot_abundance_matrix.csv")
    pq = read(OUT / "prot_qc.csv")
    assert len(prot) == len(pq) == 8427
    for p, q in zip(prot, pq):
        assert p["protein"] == q["row"], "Protein annotations are not row-aligned"
        zero = all(v in {"", "NA", "NaN"} for k, v in p.items() if k != "protein")
        assert zero == truth(q["zero_plex"]), "Protein missingness annotation mismatch"

    write(DEST / "file_catalog.csv", records)
    # A field inventory records schema and missingness. Semantic definitions are
    # maintained separately in the human-readable data dictionary.
    fields = []
    for record in records:
        path = DEST / record["file"]
        if path.suffix != ".csv":
            continue
        rows = read(path)
        if not rows:
            continue
        matrix = " by model" in record["unit"] and record["unit"] in {"gene by model", "protein feature by model"}
        for column in rows[0]:
            if matrix and column in roster:
                continue
            values = [r[column] for r in rows]
            missing = sum(v in {"", "NA", "NaN"} for v in values)
            fields.append({"file": record["file"], "column": column, "rows": len(rows),
                           "missing_values": missing, "description": FIELD_DEFINITIONS.get(column, "See file-level definition and the generating script in file_catalog.csv; source field retained without reinterpretation.")})
        if matrix:
            fields.append({"file": record["file"], "column": "<cell_line>", "rows": len(rows),
                           "missing_values": "per-column; empty/NA denotes unmeasured", "description": record["description"]})
    write(DEST / "field_inventory.csv", fields)
    shutil.copy2(ROOT / "docs/release_data_dictionary.md", DEST / "DATA_DICTIONARY.md")
    shutil.copy2(ROOT / "docs/release_readme.md", DEST / "README.md")
    summary = {"models": len(models), "patients": 34,
               **{k: sum(truth(r[k]) for r in models) for k in
                  ("rna_available", "protein_available", "copy_number_available", "variants_available")},
               "protein_features": len(prot), "validation": "passed"}
    (DEST / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    files = sorted(p for p in DEST.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (DEST / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(DEST)}\n" for p in files))
    print(json.dumps(summary, indent=2))


FIELD_DEFINITIONS = {
    "cell_line": "Stable model identifier joining matrices, metadata, and annotations.",
    "gene_id": "Ensembl gene identifier from the matched transcript reference.",
    "protein": "Protein feature identifier; representative symbol or SYMBOL|UNIPROT for secondary rows.",
    "row": "Protein feature identifier matching protein in the abundance matrix.",
    "patient_id": "Patient-of-origin grouping, not an independently verified participant identifier.",
    "selected_patient_model": "TRUE selects one model per patient by assay coverage and lexical tie-break; analysis-specific selections may differ.",
    "models_from_same_patient": "Number of resource models originating from this patient.",
    "histotype_code": "Historical analysis group: HGS, LGS, CC, EC, MC, MMMT, or SCCOHT; see histotype_note.",
    "histotype_note": "Annotation qualification; TOV112D has a published dedifferentiated carcinoma reassignment.",
    "rna_pseudoalign_pct": "Kallisto pseudoalignment percentage, 0–100.",
    "rna_sequenced_fragments_M": "Processed paired-end fragments in millions.",
    "rna_assigned_gene_counts_M": "Sum of estimated gene counts after the global expression filter, in millions.",
    "rna_genes_detected": "Globally retained genes with estimated counts >=1 in this model.",
    "prot_n_detected": "Protein features with measured abundance in this model.",
    "prot_pct_missing": "Percentage of the 8,427 protein features unmeasured in this model.",
    "wes_n_coding": "Number of retained tumour-only coding/splice candidates; not a mutation burden per megabase.",
    "fga_autosome": "Fraction of assessed autosomal segment length with absolute centred log2 ratio >0.20.",
    "cellosaurus_accession": "External Cellosaurus accession matched to the model name.",
    "str_profile_documented": "External STR reference profile exists in Cellosaurus; does not establish matching of study stocks.",
    "problematic_flag": "External Cellosaurus problematic/misidentification annotation; consult the source record.",
    "pseudoalign": "Kallisto pseudoalignment percentage, 0–100.",
    "assigned_gene_counts": "Sum of tximport estimated counts assigned to genes, after the global expression filter.",
    "n_processed_fragments": "Kallisto fragments processed; paired reads count as one fragment.",
    "detected": "Number of globally retained genes with estimated counts >=1 in a model.",
    "n_detected": "Number of protein features with a measured abundance in a model.",
    "pct_missing": "Percentage of the 8,427 protein features unmeasured in this model.",
    "symbol": "Source-table gene symbol; not all rows with the same symbol are interchangeable.",
    "symbol_representative": "Representative row chosen by completeness, peptide support, and q-value.",
    "uniprot": "Source-table UniProt accession(s).",
    "n_peptides": "Number of identified peptides reported in the source protein table.",
    "n_unique_peptides": "Number of unique peptides reported in the source protein table.",
    "npeptides_quant": "Number of quantified peptides reported in the source protein table.",
    "qvalue": "Source protein-search q-value; upstream search configuration requires laboratory confirmation.",
    "cv_replicates": "Source-reported replicate CV in percent; exact upstream definition requires confirmation.",
    "isoDoping": "Source-table isotope-doping annotation, preserved without reinterpretation.",
    "present_n_lines": "Number of resource models with measured protein abundance.",
    "present_n_plex": "Number of TMT plexes with measured protein abundance, 0–5.",
    "plexes_present": "TMT plexes in which this feature was measured.",
    "plexes_absent": "TMT plexes in which this feature was unmeasured.",
    "absent_from_any_plex": "TRUE when feature is missing from at least one TMT plex.",
    "pct_lines_missing": "Percentage of 31 protein models lacking this feature.",
    "pass_presence50": "TRUE for features measured in at least 16 of 31 models.",
    "complete_case": "TRUE for features measured in all 31 resource protein models.",
    "zero_plex": "TRUE for identified features with no quantified value in any resource plex.",
    "Chromosome": "GRCh38 chromosome for a sequence variant.",
    "Start_Position": "One-based MAF variant start coordinate.",
    "End_Position": "One-based MAF variant end coordinate; see indel representation in generating script.",
    "chromosome": "GRCh38 chromosome for a copy-number segment.",
    "start": "Zero-based inclusive segment start coordinate.",
    "end": "Zero-based exclusive segment end coordinate.",
    "log2_raw": "Target-only CNVkit CBS segment log2 ratio before project recentring (script 29).",
    "log2c": "Legacy all-chromosome probe-weighted median-centred log2 ratio.",
    "log2c_auto": "Preferred autosomal probe-weighted median-centred log2 ratio.",
    "probes": "Number of bins used for CNVkit segment estimation.",
    "weight": "Sum of CNVkit bin weights overlapping the target-only segment; can include breakpoint-filtered bins.",
    "pop_af_max": "Maximum available population allele frequency; absent annotations are not evidence of confirmed rarity.",
    "vaf": "Alternate allele fraction among supporting tumour reads; not a direct measure of LOH or somatic origin.",
    "tier": "Heuristic prioritisation tier, not a calibrated somatic probability or clinical classification.",
    "rationale": "Rule used to assign the heuristic evidence tier.",
    "context": "Source or interpretation context for a prioritised candidate.",
    "FILTER": "Archived Mutect2 filter value; released candidates passed upstream filtering.",
    "patient_rule": "Exact rule for aggregation of related models at patient level.",
    "fraction_denominator": "Full annotated non-centromeric arm length used as the arm-majority denominator.",
    "arm_reference": "Reference used to define chromosome-arm boundaries.",
}

if __name__ == "__main__":
    main()
