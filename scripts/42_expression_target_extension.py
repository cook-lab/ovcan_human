#!/usr/bin/env python3
"""Focused, descriptive RNA/protein model selection; no clinical categories.

Requires numpy/pandas. Reads current matrices only, preserves all 42 model rows,
and writes separate reports/molecular_extension_2026-09-06/expression outputs.
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
DEST = ROOT / "reports/molecular_extension_2026-09-06/expression"
DEST.mkdir(parents=True, exist_ok=True)
SEED, BOOTSTRAPS = 420906, 2000
INPUTS = {}


def read(name, **kwargs):
    p = ROOT / name
    INPUTS[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p, **kwargs)


def save(df, name):
    df.to_csv(DEST / name, index=False, na_rep="NA")


def rho(x, y):
    xy = pd.DataFrame({"x": np.asarray(x), "y": np.asarray(y)}).dropna()
    if len(xy) < 3 or xy.x.nunique() < 2 or xy.y.nunique() < 2:
        return np.nan
    return float(xy.rank(method="average").corr().iloc[0, 1])


def boot_interval(x, y, seed):
    if len(x) < 6:
        return np.nan, np.nan, 0
    rng = np.random.default_rng(seed)
    ix = rng.integers(0, len(x), (BOOTSTRAPS, len(x)))
    a = pd.DataFrame(np.asarray(x)[ix]).rank(axis=1).to_numpy()
    b = pd.DataFrame(np.asarray(y)[ix]).rank(axis=1).to_numpy()
    a -= a.mean(axis=1, keepdims=True)
    b -= b.mean(axis=1, keepdims=True)
    den = np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
    good = den > 0
    vals = (a*b).sum(axis=1)[good]/den[good]
    if not len(vals):
        return np.nan, np.nan, 0
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(lo), float(hi), len(vals)


def percentile(s, direction):
    """0..1 with the preferred end at 1; missing and constant scales stay NA."""
    if s.notna().sum() < 2 or s.nunique() < 2:
        return pd.Series(np.nan, index=s.index)
    ranks = s.rank(method="average", ascending=direction == "higher")
    return (ranks-1)/(s.notna().sum()-1)


models = read("metadata/resource_models.csv")
assert len(models) == models.cell_line.nunique() == 42
assert models.patient_id.nunique() == 34
assert not models.loc[models.selected_patient_model, "patient_id"].duplicated().any()
rna = read("output/rna_tpm.csv", index_col="gene_id")
anno = read("output/rna_gene_annotation.csv")
prot = read("output/prot_abundance_matrix.csv", index_col="protein")
qc = read("output/prot_qc.csv")
old_adc = read("output/adc_expression.csv")
panel = pd.DataFrame([
    ("ERBB2", "Surface", "higher"), ("FOLR1", "Surface", "higher"),
    ("MSLN", "Surface", "higher"), ("TACSTD2", "Surface", "higher"),
    ("SLC34A2", "Surface", "higher"), ("CDH6", "Surface", "higher"),
    ("CD276", "Surface", "higher"), ("VTCN1", "Surface", "higher"),
    ("DPEP3", "Surface", "higher"),
    ("ESR1", "Endocrine", "higher"), ("PGR", "Endocrine", "higher"),
    ("MLH1", "MMR", "lower"), ("PMS2", "MMR", "lower"),
    ("MSH2", "MMR", "lower"), ("MSH6", "MMR", "lower"),
    ("SMARCA4", "SWI/SNF", "lower"), ("SMARCA2", "SWI/SNF", "lower"),
    ("ARID1A", "SWI/SNF", "lower"), ("CCNE1", "Cell cycle", "higher"),
], columns=["gene", "panel_group", "ranking_direction"])
panel["interpretation"] = np.where(panel.ranking_direction.eq("higher"),
    "higher expression follow-up; no surface localisation, receptor activity or response claim",
    "lower expression follow-up; no protein loss, biallelic loss or functional deficiency claim")
save(panel, "panel.csv")

raw, concordance, subjects, ranking, summaries = [], [], [], [], []
for gi, spec in enumerate(panel.itertuples()):
    gene = spec.gene
    d = models.copy()
    d["gene"], d["panel_group"], d["ranking_direction"] = gene, spec.panel_group, spec.ranking_direction
    ids = anno.loc[anno.external_gene_name.eq(gene) & anno.primary_assembly, "ensembl_gene_id"].drop_duplicates()
    ids = ids[ids.isin(rna.index)]
    rv = rna.loc[ids].sum(axis=0, min_count=1) if len(ids) else pd.Series(dtype=float)
    d["rna_gene_ids"] = ";".join(ids)
    d["rna_tpm"] = d.cell_line.map(rv)
    d["rna_log2_tpm_plus1"] = np.log2(d.rna_tpm+1)
    p = qc.loc[qc.symbol.eq(gene) & qc.symbol_representative]
    assert len(p) <= 1
    p = p.iloc[0] if len(p) else None
    pv = prot.loc[p["row"]] if p is not None else pd.Series(dtype=float)
    d["protein_log2_normalised"] = d.cell_line.map(pv)
    for dest, source in [("protein_row", "row"), ("protein_uniprot", "uniprot"),
                         ("protein_isodoping", "isoDoping"), ("protein_unique_peptides_global", "n_unique_peptides"),
                         ("protein_quantified_peptides_global", "npeptides_quant")]:
        d[dest] = p[source] if p is not None else None
    d["rna_status"] = np.select([~d.rna_available, d.rna_tpm.isna(), d.rna_tpm.eq(0)],
        ["assay unavailable", "gene unavailable in filtered matrix", "measured TPM zero"], default="measured")
    d["protein_status"] = np.select([~d.protein_available, d.protein_log2_normalised.isna()],
        ["assay unavailable", "target not quantified in supplied matrix"], default="quantified")
    d["rna_abundance_percentile"] = percentile(d.rna_tpm, "higher")
    d["protein_abundance_percentile"] = percentile(d.protein_log2_normalised, "higher")
    d["rna_abundance_rank_desc"] = d.rna_tpm.rank(method="min", ascending=False)
    d["protein_abundance_rank_desc"] = d.protein_log2_normalised.rank(method="min", ascending=False)
    d["rna_n_measured_models"] = int(d.rna_tpm.notna().sum())
    d["protein_n_measured_models"] = int(d.protein_log2_normalised.notna().sum())
    raw.append(d)
    # Every paired cohort is assembled from SAME-MODEL complete pairs before
    # patient averaging; no RNA/protein borrowing across different sublines.
    pairs = d.dropna(subset=["rna_log2_tpm_plus1", "protein_log2_normalised"]).copy()
    reps = pairs.loc[pairs.selected_patient_model].copy()
    means = pairs.groupby("patient_id", as_index=False).agg(
        cell_line=("cell_line", lambda s: ";".join(sorted(s))),
        histotype_code=("histotype_code", "first"),
        rna_tpm=("rna_tpm", "mean"),
        rna_log2_tpm_plus1=("rna_log2_tpm_plus1", "mean"),
        protein_log2_normalised=("protein_log2_normalised", "mean"))
    assert pairs.groupby("patient_id").histotype_code.nunique().le(1).all()
    cohorts = [("matched_models", pairs), ("selected_patient_representatives", reps),
               ("patient_means_of_matched_models", means),
               ("HGS_selected_representatives", reps.loc[reps.histotype_code.eq("HGS")]),
               ("CC_selected_representatives", reps.loc[reps.histotype_code.eq("CC")])]
    for ci, (cohort, frame) in enumerate(cohorts):
        x, y = frame.rna_log2_tpm_plus1, frame.protein_log2_normalised
        r = rho(x, y)
        if cohort in ["selected_patient_representatives", "HGS_selected_representatives", "CC_selected_representatives"] and np.isfinite(r):
            lo, hi, valid = boot_interval(x, y, SEED+100*gi+ci)
        else:
            lo, hi, valid = np.nan, np.nan, 0
        concordance.append(dict(gene=gene, panel_group=spec.panel_group, cohort=cohort,
            n_observations=len(frame), n_patients=frame.patient_id.nunique(), spearman=r,
            bootstrap_lo=lo, bootstrap_hi=hi, bootstrap_valid=valid,
            interval_method="patient-resampling percentile 95%; descriptive convenience-panel sensitivity" if valid else "not estimated",
            protein_isodoping=p.isoDoping if p is not None else None,
            interpretation="unadjusted association; subtype/centre/plex/passage not controlled"))
        f = frame[["cell_line", "patient_id", "histotype_code", "rna_tpm", "rna_log2_tpm_plus1", "protein_log2_normalised"]].copy()
        f["cohort"], f["gene"] = cohort, gene
        subjects.append(f)
        if not len(f):
            continue
        f["ranking_direction"] = spec.ranking_direction
        f["rna_preference_percentile"] = percentile(f.rna_log2_tpm_plus1, spec.ranking_direction)
        f["protein_preference_percentile"] = percentile(f.protein_log2_normalised, spec.ranking_direction)
        f["joint_mean_score"] = (f.rna_preference_percentile+f.protein_preference_percentile)/2
        f["joint_weaker_assay_score"] = f[["rna_preference_percentile", "protein_preference_percentile"]].min(axis=1, skipna=False)
        f["joint_rna25_score"] = .25*f.rna_preference_percentile+.75*f.protein_preference_percentile
        f["joint_rna75_score"] = .75*f.rna_preference_percentile+.25*f.protein_preference_percentile
        for method in ["mean", "weaker_assay", "rna25", "rna75"]:
            f[f"joint_{method}_rank"] = f[f"joint_{method}_score"].rank(method="min", ascending=False)
        f["rank_best_across_methods"] = f[[f"joint_{m}_rank" for m in ["mean", "weaker_assay", "rna25", "rna75"]]].min(axis=1)
        f["rank_worst_across_methods"] = f[[f"joint_{m}_rank" for m in ["mean", "weaker_assay", "rna25", "rna75"]]].max(axis=1)
        f["absolute_assay_percentile_gap"] = abs(f.rna_preference_percentile-f.protein_preference_percentile)
        f["both_assays_same_direction"] = (f.rna_preference_percentile.ge(.5) & f.protein_preference_percentile.ge(.5))
        f["n_ranked"] = len(f)
        f["protein_isodoping"] = p.isoDoping if p is not None else None
        ranking.append(f)
    summary = dict(gene=gene, panel_group=spec.panel_group, ranking_direction=spec.ranking_direction,
        n_models=42, n_rna_models=int(d.rna_tpm.notna().sum()), n_protein_models=int(d.protein_log2_normalised.notna().sum()),
        n_rna_zero_models=int(d.rna_tpm.eq(0).sum()), n_rna_positive_models=int(d.rna_tpm.gt(0).sum()),
        n_paired_models=len(pairs), n_paired_patient_representatives=len(reps),
        protein_isodoping=p.isoDoping if p is not None else None,
        protein_unique_peptides_global=p.n_unique_peptides if p is not None else None)
    for col in ["rna_tpm", "protein_log2_normalised"]:
        values = d[col].dropna()
        for suffix, value in [("min", values.min()), ("median", values.median()), ("max", values.max())]:
            summary[f"{col}_{suffix}"] = value
    summaries.append(summary)

raw = pd.concat(raw, ignore_index=True)
concordance = pd.DataFrame(concordance)
ranking = pd.concat(ranking, ignore_index=True)
subjects = pd.concat(subjects, ignore_index=True)
assert len(raw) == 42*len(panel) and not raw.duplicated(["gene", "cell_line"]).any()
assert not raw.loc[~raw.rna_available, "rna_tpm"].notna().any()
assert not raw.loc[~raw.protein_available, "protein_log2_normalised"].notna().any()
assert ranking[["rna_log2_tpm_plus1", "protein_log2_normalised"]].notna().all().all()
save(raw, "model_expression.csv")
save(raw.loc[raw.rna_tpm.notna() ^ raw.protein_log2_normalised.notna()], "unpaired_expression.csv")
save(pd.DataFrame(summaries), "target_summary.csv")
save(concordance, "concordance.csv")
save(subjects, "cohort_subjects.csv")
save(ranking, "joint_ranking.csv")

# Present the top three joint ranks per patient cohort, without making a
# positive/negative or response classification. Ties may produce >3 rows.
leads = ranking.loc[ranking.cohort.eq("selected_patient_representatives") & ranking.joint_mean_rank.le(3)].copy()
for cname, prefix in [("patient_means_of_matched_models", "patient_mean"),
                      ("HGS_selected_representatives", "HGS"), ("CC_selected_representatives", "CC")]:
    z = ranking.loc[ranking.cohort.eq(cname), ["gene", "patient_id", "joint_mean_rank", "n_ranked"]]
    z = z.rename(columns={"joint_mean_rank": prefix+"_rank", "n_ranked": prefix+"_n"})
    leads = leads.merge(z, on=["gene", "patient_id"], how="left", validate="one_to_one")
leads = leads.merge(raw[["gene", "cell_line", "rna_passage", "wes_passage", "tmt_plex", "tmt_channel", "protein_unique_peptides_global"]],
                    on=["gene", "cell_line"], validate="one_to_one")
save(leads.sort_values(["gene", "joint_mean_rank", "cell_line"]), "model_leads.csv")
discordance = ranking.loc[ranking.cohort.eq("selected_patient_representatives")].copy()
discordance["discordance_rank"] = discordance.groupby("gene").absolute_assay_percentile_gap.rank(method="min", ascending=False)
save(discordance.loc[discordance.discordance_rank.le(3)].sort_values(["gene", "discordance_rank"]), "discordant_models.csv")

# Independent checkpoint against existing atlas: use its three-decimal precision.
overlap = old_adc.merge(raw, left_on=["symbol", "cell_line"], right_on=["gene", "cell_line"], how="inner")
ov = np.where(overlap.assay.eq("RNA"), overlap.rna_log2_tpm_plus1, overlap.protein_log2_normalised)
delta = abs(ov-overlap.log2_expr)
finite = np.isfinite(delta)
assert finite.any() and float(delta[finite].max()) <= .000501
qa = {"script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "inputs_sha256": INPUTS,
    "n_genes": len(panel), "n_model_expression_rows": len(raw), "all_models": 42, "patients": 34,
    "all_42_models_retained_per_gene": bool(raw.groupby("gene").cell_line.nunique().eq(42).all()),
    "existing_ADC_agreement_comparisons": int(finite.sum()), "existing_ADC_max_abs_difference": float(delta[finite].max()),
    "bootstrap_seed": SEED, "bootstrap_replicates": BOOTSTRAPS,
    "paired_before_patient_aggregation": True, "missing_values_imputed": False,
    "clinical_categories_assigned": False, "canonical_outputs_modified": False,
    "versions": {"numpy": np.__version__, "pandas": pd.__version__}}
(DEST/"analysis_qa.json").write_text(json.dumps(qa, indent=2)+"\n")
print(json.dumps({k: v for k, v in qa.items() if k != "inputs_sha256"}, indent=2))
