#!/usr/bin/env python3
"""Local depth checks for gene-level gains/depletion; no copy-state calls."""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(os.environ.get("OVCAN_PROJ",Path.cwd())).resolve()
OUT=ROOT/"reports/molecular_extension_2026-09-06/copy_number"
BUNDLE=Path(os.environ.get("OVCAN_COVERAGE_BUNDLE",ROOT/"data/cluster_wes_retrieval/2026-09-06/ovcan_human_wes_cnv_coverage_2026-09-06"))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
genes=["AKT2","CCNE1","CCND1","ERBB2","NF1","CDKN2A","CDKN2B"]
bins=pd.read_csv(OUT/"gene_target_bins.csv")
bins=bins.loc[bins.gene_symbol.isin(genes)]
keys=["chromosome","start","end"]
coords=bins[["gene_symbol"]+keys].drop_duplicates()
assert not coords.duplicated(["gene_symbol"]+keys).any()
sources=pd.read_csv(ROOT/"output/wes_cnv_coverage_sources.csv")
normal_sources=sources.loc[sources.path.str.match(r"normals/SRR\d+\.sorted\.targetcoverage\.cnn$")]
assert len(normal_sources)==5
norm=[]; hashes={}
for row in normal_sources.itertuples():
    p=BUNDLE/row.path
    assert sha(p)==row.sha256
    hashes[row.path]=row.sha256
    d=pd.read_csv(p,sep="\t")
    d=coords.merge(d[keys+["depth"]],on=keys,how="left",validate="many_to_one")
    assert d.depth.notna().all()
    d["sample_id"]=Path(row.path).name.split(".")[0]
    norm.append(d)
norm=pd.concat(norm,ignore_index=True)
norm.to_csv(OUT/"selected_loci_normal_bin_depth.csv",index=False)
normal_summary=norm.groupby(["gene_symbol"]+keys).depth.agg(normal_n="count",normal_min_depth="min",normal_median_depth="median",normal_max_depth="max").reset_index()
combined=bins.merge(normal_summary,on=["gene_symbol"]+keys,validate="many_to_one")
combined.to_csv(OUT/"selected_loci_model_normal_bins.csv",index=False)
rows=[]
for (model,gene),d in combined.groupby(["cell_line","gene_symbol"]):
    rows.append(dict(cell_line=model,gene=gene,n_bins=len(d),model_zero_bins=int(d.depth.eq(0).sum()),
        model_median_depth=d.depth.median(),model_min_depth=d.depth.min(),model_max_depth=d.depth.max(),
        median_reference_normal_bin_depth=d.normal_median_depth.median(),
        min_any_reference_normal_bin_depth=d.normal_min_depth.min(),
        n_bins_positive_all_five_normals=int(d.normal_min_depth.gt(0).sum()),
        positive_model_bin_median_log2=d.loc[d.depth.gt(0)&d.weight.gt(0),"log2c_auto"].median()))
result=pd.DataFrame(rows)
result.to_csv(OUT/"selected_loci_depth_summary.csv",index=False)
check=dict(normal_sha256=hashes,source_gene_bins_sha256=sha(OUT/"gene_target_bins.csv"),genes=genes,
           rows=len(result),models=result.cell_line.nunique(),normal_samples=5,
           note="Same retained target-bin coordinates; normal/model preprocessing may differ. Zero depth is evidence for review, not a validated homozygous deletion.")
(OUT/"depth_review_validation.json").write_text(json.dumps(check,indent=2)+"\n")
print(result.loc[(result.gene.isin(["NF1","CDKN2A"])) & result.model_zero_bins.gt(0)].round(4).to_string(index=False))
