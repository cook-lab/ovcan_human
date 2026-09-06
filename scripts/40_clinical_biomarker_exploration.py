#!/usr/bin/env python3
"""Exploratory model-selection annotations; never clinical positive/negative calls.

Uses current corrected CNV segments, hash-verified original target CNRs, and the
unchanged RNA/protein matrices. Outputs are separate from the manuscript release.
Requires numpy/pandas. Run from the repository root; OVCAN_PROJ/DATA supported.
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
DEST = ROOT / "reports/clinical_classification_2026-09-06"
DEST.mkdir(parents=True, exist_ok=True)
INPUTS = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, **kwargs):
    p = ROOT / path
    INPUTS[path] = sha(p)
    return pd.read_csv(p, **kwargs)


def save(frame, name):
    frame.to_csv(DEST / name, index=False, na_rep="NA")


def weighted_median(values, weights):
    """Lower inverse-CDF weighted median; defined here, not a new CNV center."""
    values, weights = np.asarray(values), np.asarray(weights)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[ok], weights[ok]
    order = np.argsort(values, kind="stable")
    values, weights = values[order], weights[order]
    return float(values[np.searchsorted(np.cumsum(weights), weights.sum()/2)])


models = read("metadata/resource_models.csv")
assert len(models) == 42 and models.patient_id.nunique() == 34
rna = read("output/rna_tpm.csv").set_index("gene_id")
anno = read("output/rna_gene_annotation.csv")
protein = read("output/prot_abundance_matrix.csv").set_index("protein")
protein_qc = read("output/prot_qc.csv")
segments = read("output/wes_cnv_segments.csv")
cnr_map = read("output/wes_recovered_provenance_cnv_support.csv")
target_manifest = read("output/wes_cnv_target_only/manifest.csv")
arms = read("output/wes_cnv_arm_boundaries.csv")
assert segments.cell_line.nunique() == len(cnr_map) == len(target_manifest) == 23
for item in target_manifest.itertuples():
    assert sha(ROOT / item.cns_path) == item.cns_sha256

# NCBI Gene GRCh38.p14 RS_2025_08, verified 2026-09-06. Convert 1-based closed
# gene coordinates to 0-based half-open intervals before any overlap operation.
loci = pd.DataFrame([
    ["CCNE1", "chr19", 29811990, 29824312, "19q", "898"],
    ["ERBB2", "chr17", 39688093, 39728658, "17q", "2064"],
], columns=["gene", "chromosome", "start", "end", "arm", "ncbi_gene_id"])
loci["assembly"] = "GRCh38.p14"
loci["coordinate_convention"] = "0-based half-open"
loci["source"] = "https://www.ncbi.nlm.nih.gov/gene/" + loci.ncbi_gene_id
loci["reference_annotation"] = "RS_2025_08; accessed 2026-09-06"
save(loci, "locus_definitions.csv")

panel = "CCNE1 ERBB2 FOLR1 ESR1 PGR MSLN TACSTD2 CD274 BRCA1 BRCA2 PALB2 RAD51C RAD51D BRIP1 ATM ATR CHEK2 CDK12 MLH1 MSH2 MSH6 PMS2 POLE POLD1 ARID1A SMARCA4 SMARCA2 RB1 PTEN KRAS NRAS BRAF PIK3CA TP53 FBXW7 MYC".split()
expression = []
for gene in panel:
    d = models[["cell_line", "histotype_code", "patient_id", "rna_passage", "wes_passage", "tmt_plex"]].copy()
    d["gene"] = gene
    ids = anno.loc[(anno.external_gene_name == gene) & anno.primary_assembly, "ensembl_gene_id"]
    ids = ids[ids.isin(rna.index)]
    d["rna_gene_ids"] = ";".join(ids)
    # Sum TPM only for multiple primary-assembly Ensembl loci sharing a symbol.
    rr = rna.loc[ids].sum(axis=0, min_count=1) if len(ids) else pd.Series(dtype=float)
    d["rna_tpm"] = d.cell_line.map(rr)
    d["rna_log2_tpm_plus1"] = np.log2(d.rna_tpm + 1)
    d["rna_rank_all_models_desc"] = d.rna_tpm.rank(method="min", ascending=False)
    d["rna_n_measured_models"] = d.rna_tpm.notna().sum()
    d["rna_rank_within_histotype_desc"] = d.groupby("histotype_code").rna_tpm.rank(method="min", ascending=False)
    d["rna_n_measured_within_histotype"] = d.groupby("histotype_code").rna_tpm.transform("count")
    # Equal patient weighting for reference distributions; no subline borrowing.
    patients = d.groupby(["histotype_code", "patient_id"]).rna_log2_tpm_plus1.mean().reset_index()
    for hist, ix in d.groupby("histotype_code").groups.items():
        values = patients.loc[patients.histotype_code.eq(hist), "rna_log2_tpm_plus1"].dropna()
        median = values.median()
        d.loc[ix, "rna_histotype_patient_median_log2_tpm_plus1"] = median
        d.loc[ix, "rna_histotype_measured_patients"] = len(values)
        d.loc[ix, "rna_log2_tpm_plus1_delta_histotype_patient_median"] = d.loc[ix, "rna_log2_tpm_plus1"] - median
    pr = protein_qc.loc[protein_qc.symbol.eq(gene) & protein_qc.symbol_representative]
    assert len(pr) <= 1, (gene, "multiple representative protein rows")
    prow = pr.iloc[0] if len(pr) else None
    pp = protein.loc[prow["row"]] if prow is not None else pd.Series(dtype=float)
    d["protein_row"] = prow["row"] if prow is not None else None
    d["protein_uniprot"] = prow["uniprot"] if prow is not None else None
    d["protein_isodoping"] = bool(prow.isoDoping) if prow is not None else None
    d["protein_unique_peptides"] = prow.n_unique_peptides if prow is not None else None
    d["protein_log2_normalised"] = d.cell_line.map(pp)
    d["protein_rank_all_models_desc"] = d.protein_log2_normalised.rank(method="min", ascending=False)
    d["protein_n_measured_models"] = d.protein_log2_normalised.notna().sum()
    d["interpretation"] = "descriptive abundance; no clinical expression threshold; isodoping flagged where supplied"
    expression.append(d)
expression = pd.concat(expression, ignore_index=True)
save(expression, "marker_expression.csv")

cnv, gene_bins, region_segments = [], [], []
for item in cnr_map.itertuples():
    path = ROOT / item.cnr_source
    external_data = os.environ.get("OVCAN_DATA")
    if external_data:
        path = Path(external_data) / Path(item.cnr_source).relative_to("judy_archive/data")
    assert sha(path) == item.cnr_sha256, f"Source CNR changed: {item.cell_line}"
    INPUTS[item.cnr_source] = item.cnr_sha256
    bins = pd.read_csv(path, sep="\t")
    bins = bins.loc[bins.gene.ne("Antitarget")].copy()
    s = segments.loc[segments.cell_line.eq(item.cell_line)].copy()
    offset = s.log2_raw - s.log2c_auto
    assert offset.max() - offset.min() < 1e-10
    center = float(offset.iloc[0])
    for locus in loci.itertuples():
        hits = s.loc[s.chromosome.eq(locus.chromosome) & s.start.lt(locus.end) & s.end.gt(locus.start)].copy()
        overlap = np.minimum(hits.end, locus.end) - np.maximum(hits.start, locus.start)
        assert len(hits) and overlap.sum() == locus.end-locus.start, (item.cell_line, locus.gene)
        assert (overlap > 0).all()
        log2 = float(np.average(hits.log2c_auto, weights=overlap))
        raw = float(np.average(hits.log2_raw, weights=overlap))
        b = bins.loc[bins.chromosome.eq(locus.chromosome) & bins.start.lt(locus.end) & bins.end.gt(locus.start)].copy()
        b["gene_symbol"] = locus.gene
        b["cell_line"] = item.cell_line
        b["log2c_auto"] = b.log2 - center
        gene_bins.append(b)
        positive = b.loc[b.depth.gt(0) & b.weight.gt(0)]
        arm = arms.loc[arms.arm.eq(locus.arm)].iloc[0]
        arm_seg = s.loc[s.chromosome.eq(locus.chromosome) & s.start.lt(arm.arm_end) & s.end.gt(arm.arm_start)]
        arm_weights = np.minimum(arm_seg.end, arm.arm_end)-np.maximum(arm_seg.start, arm.arm_start)
        arm_med = weighted_median(arm_seg.log2c_auto, arm_weights)
        cnv.append(dict(cell_line=item.cell_line, gene=locus.gene, chromosome=locus.chromosome,
            gene_start=locus.start, gene_end=locus.end, segments_overlapping_gene=len(hits),
            segment_start_min=int(hits.start.min()), segment_end_max=int(hits.end.max()),
            locus_segment_span_mb=float((hits.end.max()-hits.start.min())/1e6),
            gene_segment_log2_raw=raw, gene_segment_log2c_auto=log2,
            ratio_to_autosomal_baseline=2**log2, autosomal_center=center,
            arm_length_weighted_median_log2=arm_med, locus_minus_arm_log2=log2-arm_med,
            target_bins_overlapping_gene=len(b), positive_target_bins=len(positive),
            gene_target_median_log2c_auto=positive.log2c_auto.median(),
            gene_target_min_log2c_auto=positive.log2c_auto.min(),
            gene_target_max_log2c_auto=positive.log2c_auto.max(),
            strong_relative_gain_screen=bool(log2 >= 1),
            clinical_amplification_status="not established: absolute copy number/ploidy unavailable"))
    r = s.loc[s.chromosome.eq("chr19")].copy()
    region_segments.append(r)
cnv = pd.DataFrame(cnv).merge(models[["cell_line", "histotype_code", "patient_id"]], on="cell_line", validate="many_to_one")
save(cnv, "cnv_locus_summary.csv")
save(pd.concat(gene_bins, ignore_index=True), "locus_target_bin_evidence.csv")
save(pd.concat(region_segments, ignore_index=True), "chr19_segments.csv")
ccne = expression.loc[expression.gene.eq("CCNE1")].merge(cnv.drop(columns=["histotype_code", "patient_id"]), on=["cell_line", "gene"], how="left", validate="one_to_one")
ccne["wes_available"] = ccne.gene_segment_log2c_auto.notna()
ccne["clinical_amplification_status"] = ccne.clinical_amplification_status.fillna("not assessed: WES unavailable")
assert len(ccne) == 42 and ccne.rna_tpm.notna().sum() == ccne.protein_log2_normalised.notna().sum() == 31
assert ccne.wes_available.sum() == 23
save(ccne, "ccne1_model_summary.csv")

variants = read("output/wes_mutations_filtered.csv")
tiers = read("output/wes_driver_tiers.csv")
keys = ["cell_line", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
v = variants.loc[variants.Hugo_Symbol.isin(panel)].copy()
assert not tiers.duplicated(keys).any()
v = v.merge(tiers[keys + ["tier", "rationale"]], on=keys, how="left", validate="many_to_one")
v = v.merge(models[["cell_line", "patient_id"]], on="cell_line", validate="many_to_one")
v["clinical_reclassification"] = "not performed; preserved research-prioritisation tier is not an AMP/ASCO/CAP clinical tier"
v["somatic_origin"] = "unconfirmed: tumour-only calling"
v["biallelic_status"] = "not established"
save(v, "biomarker_variant_candidates.csv")
summary = {
    "purpose": "Exploratory cell-model selection; no clinical positive/negative calls",
    "script_sha256": sha(Path(__file__)), "input_sha256": INPUTS,
    "models": len(models), "patients": int(models.patient_id.nunique()),
    "ccne1_cnv_models": int(ccne.wes_available.sum()), "ccne1_rna_models": int(ccne.rna_tpm.notna().sum()),
    "ccne1_protein_models": int(ccne.protein_log2_normalised.notna().sum()),
    "ccne1_target_bin_counts": sorted(cnv.loc[cnv.gene.eq("CCNE1"), "positive_target_bins"].unique().tolist()),
    "ccne1_strong_relative_gain_screen_models": cnv.loc[cnv.gene.eq("CCNE1") & cnv.strong_relative_gain_screen].sort_values("gene_segment_log2c_auto", ascending=False).cell_line.tolist(),
    "screen_definition": "gene-overlap length-weighted corrected segment log2 >= 1; prioritisation only, not absolute CN or clinical amplification",
    "marker_genes": panel, "marker_expression_rows": len(expression), "retained_panel_candidates": len(v),
    "coordinate_convention": "0-based half-open (gene and segment overlaps)",
    "canonical_outputs_modified": False,
}
(DEST / "analysis_validation.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps({k:v for k,v in summary.items() if k not in ["input_sha256", "marker_genes"]}, indent=2))
