#!/usr/bin/env python3
"""Check selected gene-span depth findings against canonical transcript exons."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import urllib.request

import pandas as pd

ROOT=Path(os.environ.get("OVCAN_PROJ",Path.cwd())).resolve()
OUT=ROOT/"reports/molecular_extension_2026-09-06/copy_number"
GENES={"NF1":"ENSG00000196712","CDKN2A":"ENSG00000147889","AKT2":"ENSG00000105221"}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fetch-exons",action="store_true")
args=parser.parse_args()
lookup=json.loads((OUT/"ensembl_gene_lookup.json").read_text())
cache=OUT/"canonical_transcript_exons.json"
if args.fetch_exons:
    data={}
    for gene,ens in GENES.items():
        tid=lookup[ens]["canonical_transcript"]
        url="https://rest.ensembl.org/lookup/id/"+tid.split(".")[0]+"?expand=1;content-type=application/json"
        with urllib.request.urlopen(url,timeout=30) as r:
            item=json.load(r)
        assert item["assembly_name"]=="GRCh38" and item["Parent"]==ens
        assert item["id"]+"."+str(item["version"])==tid
        data[gene]={"lookup_url":url,"canonical_transcript":tid,"response":item}
    cache.write_text(json.dumps(data,indent=2)+"\n")
    print("Pinned canonical transcript exons for",", ".join(data))
    raise SystemExit
data=json.loads(cache.read_text())
bins=pd.read_csv(OUT/"selected_loci_model_normal_bins.csv")
out=[]; details=[]
for gene,item in data.items():
    exons=item["response"]["Exon"]
    b=bins.loc[bins.gene_symbol.eq(gene)].copy()
    b["canonical_exon_overlap"]=False
    for exon in exons:
        b.loc[b.start.lt(exon["end"]) & b.end.gt(exon["start"]-1),"canonical_exon_overlap"]=True
    b["canonical_transcript"]=item["canonical_transcript"]
    details.append(b)
    for model,d in b.groupby("cell_line"):
        e=d.loc[d.canonical_exon_overlap]
        assert len(e)>0
        out.append(dict(cell_line=model,gene=gene,canonical_transcript=item["canonical_transcript"],
            gene_span_bins=len(d),canonical_exon_bins=len(e),excluded_noncanonical_or_intronic_bins=len(d)-len(e),
            zero_depth_canonical_exon_bins=int(e.depth.eq(0).sum()),canonical_exon_median_depth=e.depth.median(),
            canonical_exon_max_depth=e.depth.max(),min_normal_depth_across_canonical_bins=e.normal_min_depth.min(),
            positive_canonical_exon_bin_median_log2=e.loc[e.depth.gt(0)&e.weight.gt(0),"log2c_auto"].median()))
pd.concat(details,ignore_index=True).to_csv(OUT/"canonical_exon_bin_evidence.csv",index=False)
summary=pd.DataFrame(out)
summary.to_csv(OUT/"canonical_exon_depth_summary.csv",index=False)
(OUT/"exon_review_validation.json").write_text(json.dumps({"reference_sha256":sha(cache),"source_bins_sha256":sha(OUT/"selected_loci_model_normal_bins.csv"),"genes":GENES,"model_gene_rows":len(summary),"scope":"Canonical transcript exons, including UTRs; other isoforms and regulatory regions are not ruled out."},indent=2)+"\n")
print(summary.loc[summary.zero_depth_canonical_exon_bins.gt(0)].round(4).to_string(index=False))
