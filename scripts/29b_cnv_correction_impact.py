#!/usr/bin/env python3
"""Compare current target-only CNV summaries with the pinned public v8 revision.

Run after scripts 29, 08, 10 and 20. Requires pandas/numpy and local Git history;
does not fetch, mutate Git, or modify scientific outputs. The baseline is fixed
at the reviewed pre-correction commit, rather than whatever HEAD happens to be.
"""
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = "854def3653918034055ac0e8e94907c852da2484"
REPORT = ROOT / "reports/wes_cnv_coverage_2026-09-06"
SOURCES = []


def pair(name):
    old = subprocess.check_output(["git", "show", f"{BASE}:output/{name}"], cwd=ROOT)
    new = (ROOT / "output" / name).read_bytes()
    SOURCES.append(dict(file=name, archived_sha256=hashlib.sha256(old).hexdigest(),
                        current_sha256=hashlib.sha256(new).hexdigest()))
    return pd.read_csv(io.BytesIO(old)), pd.read_csv(io.BytesIO(new))


def changes(old, new, keys, fields):
    merged = old[keys + fields].merge(new[keys + fields], on=keys, suffixes=("_archived", "_target_only"), validate="one_to_one")
    assert len(merged) == len(old) == len(new)
    changed = np.zeros(len(merged), dtype=bool)
    for f in fields:
        changed |= merged[f + "_archived"].fillna("<NA>") != merged[f + "_target_only"].fillna("<NA>")
    return merged.loc[changed]


def main():
    REPORT.mkdir(exist_ok=True)
    old, new = pair("wes_cnv_fga.csv")
    assert len(old) == len(new) == 23
    cols = ["cell_line", "subtype", "patient_id", "fga_auto_0.2", "fga_auto_0.3",
            "total_assessed_autosome_mb", "n_segments_auto", "chrX_median_log2c"]
    fga = old[cols].merge(new[cols], on=["cell_line", "subtype", "patient_id"],
                         suffixes=("_archived", "_target_only"), validate="one_to_one")
    assert len(fga) == 23 and set(fga.cell_line) == set(old.cell_line) == set(new.cell_line)
    fga["fga_auto_0.2_delta"] = fga["fga_auto_0.2_target_only"] - fga["fga_auto_0.2_archived"]
    fga.to_csv(REPORT / "correction_fga_comparison.csv", index=False)
    patient = fga.groupby(["subtype", "patient_id"])[["fga_auto_0.2_archived", "fga_auto_0.2_target_only"]].mean().reset_index()
    patient.to_csv(REPORT / "correction_patient_fga.csv", index=False)
    arm_old, arm_new = pair("wes_cnv_arm_calls.csv")
    arm_changes = changes(arm_old, arm_new, ["cell_line", "arm"], ["call"])
    arm_changes.to_csv(REPORT / "correction_arm_call_changes.csv", index=False)
    af_old, af_new = pair("wes_cnv_arm_freq_patient.csv")
    af_changes = changes(af_old, af_new, ["arm"], ["n_lines_gain", "n_lines_loss", "n_patients_gain", "n_patients_loss"])
    af_changes.to_csv(REPORT / "correction_arm_frequency_changes.csv", index=False)
    sens_old, sens_new = pair("wes_cnv_arm_freq_sensitivity.csv")
    sens_changes = changes(sens_old, sens_new, ["arm", "log2_threshold", "arm_majority_frac"], ["n_patients_gain", "n_patients_loss"])
    sens_changes.to_csv(REPORT / "correction_arm_sensitivity_changes.csv", index=False)
    auth_old, auth_new = pair("auth_perline_table.csv")
    auth_changes = changes(auth_old, auth_new, ["cell_line"], ["cnv_instability", "genomics_consistent"])
    auth_changes.to_csv(REPORT / "correction_authentication_changes.csv", index=False)
    sup_old, sup_new = pair("supplement_per_line.csv")
    sup_changes = changes(sup_old, sup_new, ["cell_line"], ["fga_autosome", "high_fga_flag"])
    sup_changes.to_csv(REPORT / "correction_supplement_changes.csv", index=False)
    seg_old, seg_new = pair("wes_cnv_segments.csv")
    assert set(seg_old.cell_line) == set(seg_new.cell_line) == set(old.cell_line)
    # Point-locus checks use precisely the locations in script 08. This is a
    # sensitivity comparison, not evidence for gene-level focal amplification.
    loci = [("MECOM/PRKCI", "chr3", 169500000, "gain"), ("SOX2", "chr3", 181700000, "gain"),
            ("MYC", "chr8", 127700000, "gain"), ("CCNE1", "chr19", 29800000, "gain"),
            ("20q", "chr20", 50000000, "gain"), ("TP53", "chr17", 7700000, "loss"),
            ("RB1", "chr13", 48300000, "loss"), ("PTEN", "chr10", 87900000, "loss")]
    locus_rows = []
    for label, chrom, pos, direction in loci:
        for line in sorted(old.loc[old.subtype == "HGS", "cell_line"]):
            row = dict(cell_line=line, locus=label, chromosome=chrom, position=pos, direction=direction)
            for name, df in [("archived", seg_old), ("target_only", seg_new)]:
                v = df.loc[(df.cell_line == line) & (df.chromosome == chrom) & (df.start < pos) & (df.end >= pos), "log2c_auto"]
                row[name + "_hit"] = bool((v > 0.2).any() if direction == "gain" else (v < -0.2).any())
            locus_rows.append(row)
    locus = pd.DataFrame(locus_rows)
    locus_changes = locus[locus.archived_hit != locus.target_only_hit]
    locus_changes.to_csv(REPORT / "correction_locus_changes.csv", index=False)
    hgs = patient[patient.subtype == "HGS"]
    summary = dict(baseline_commit=BASE, models=23, patients=int(new.patient_id.nunique()),
                   hgs_patients=len(hgs),
                   hgs_patient_median_fga_archived=float(hgs["fga_auto_0.2_archived"].median()),
                   hgs_patient_median_fga_target_only=float(hgs["fga_auto_0.2_target_only"].median()),
                   maximum_absolute_model_fga_delta=float(fga["fga_auto_0.2_delta"].abs().max()),
                   model_arm_call_changes=len(arm_changes), arm_frequency_rows_changed=len(af_changes),
                   patient_sensitivity_rows_changed=len(sens_changes), authentication_class_changes=len(auth_changes),
                   hgsc_point_locus_hit_changes=len(locus_changes), main_chromosome_segments_archived=len(seg_old),
                   main_chromosome_segments_target_only=len(seg_new), sources=SOURCES)
    (REPORT / "correction_impact.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "sources"}, indent=2))


if __name__ == "__main__":
    main()
