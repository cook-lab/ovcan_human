#!/usr/bin/env python3
"""Target-supported relative gene CN screen; not absolute amplification/deletion.

Run --fetch-reference once to retrieve and pin the small Ensembl gene lookup.
Normal runs use that saved reference without network. Original data are read-only.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import urllib.request
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
DEST = ROOT / "reports/molecular_extension_2026-09-06/copy_number"
DEST.mkdir(parents=True, exist_ok=True)
GENES = "CCNE1 CCND1 CDK4 CDK6 MYC MYCN MDM2 MDM4 ERBB2 EGFR MET FGFR1 FGFR2 FGFR3 KRAS NRAS BRAF AKT2 PIK3CA PRKCI MECOM BCL2L1 ZNF217 SOX2 NF1 RB1 PTEN CDKN2A CDKN2B BRCA1 BRCA2 PALB2 RAD51C RAD51D BRIP1 ATM CHEK2 CDK12 ARID1A SMARCA4 SMARCA2 MLH1 MSH2 MSH6 PMS2 TP53 FBXW7 POLE POLD1 CTNNB1 ARID1B SMAD4".split()
GAIN_GENES = set(GENES[:24])
INPUTS = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p, **kw):
    INPUTS[p] = sha(ROOT/p)
    return pd.read_csv(ROOT/p, **kw)


def save(df, name):
    df.to_csv(DEST/name, index=False, na_rep="NA")


def source(p):
    return Path(os.environ["OVCAN_DATA"])/Path(p).relative_to("judy_archive/data") if os.environ.get("OVCAN_DATA") else ROOT/p


def wmedian(v, w):
    v, w = np.asarray(v), np.asarray(w)
    good = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[good], w[good]
    if not len(v):
        return np.nan
    order = np.argsort(v, kind="stable")
    return float(v[order][np.searchsorted(np.cumsum(w[order]),w.sum()/2)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch-reference", action="store_true")
    args = ap.parse_args()
    ann = read("output/rna_gene_annotation.csv")
    ids = {}
    for g in GENES:
        z = ann.loc[ann.external_gene_name.eq(g) & ann.primary_assembly, "ensembl_gene_id"]
        assert len(z) == 1, (g, z.tolist())
        ids[g] = z.iloc[0]
    raw_path = DEST/"ensembl_gene_lookup.json"
    if args.fetch_reference:
        payload = json.dumps({"ids":list(ids.values())}).encode()
        url = "https://rest.ensembl.org/lookup/id"
        req = urllib.request.Request(url, data=payload, headers={"Content-Type":"application/json", "Accept":"application/json"})
        with urllib.request.urlopen(req, timeout=45) as response:
            raw = response.read()
        look = json.loads(raw)
        assert set(look) == set(ids.values()) and all(look.values())
        with urllib.request.urlopen("https://rest.ensembl.org/info/data?content-type=application/json", timeout=30) as r:
            release = json.load(r)
        raw_path.write_bytes(raw)
        (DEST/"reference_provenance.json").write_text(json.dumps({
            "lookup_url":url, "query_ids":ids, "query_sha256":hashlib.sha256(payload).hexdigest(),
            "response_sha256":sha(raw_path), "info_data":release,
            "retrieved_utc":datetime.now(timezone.utc).isoformat(),
            "use":"New locus boundaries only; RNA quantification remains matched Ensembl93",
            "preserved_loci":"CCNE1 and ERBB2 use previously pinned NCBI coordinates for continuity"
        }, indent=2)+"\n")
        print("Pinned gene lookup:", len(look), "genes;", release)
        return
    provenance = json.loads((DEST/"reference_provenance.json").read_text())
    assert sha(raw_path) == provenance["response_sha256"]
    lookup = json.loads(raw_path.read_text())
    prior_loci = read("reports/clinical_classification_2026-09-06/locus_definitions.csv").set_index("gene")
    arms = read("output/wes_cnv_arm_boundaries.csv")
    loci = []
    for gene, ens in ids.items():
        obj = lookup[ens]
        assert obj["assembly_name"] == "GRCh38" and obj["display_name"] == gene
        chrom, start, end = "chr"+obj["seq_region_name"], obj["start"]-1, obj["end"]
        url, version = "https://rest.ensembl.org/lookup/id/"+ens, str(provenance["info_data"])
        if gene in prior_loci.index:
            old = prior_loci.loc[gene]
            chrom, start, end, url, version = old.chromosome, int(old.start), int(old.end), old.source, old.reference_annotation
        arm_hits = arms.loc[arms.chromosome.eq(chrom) & arms.arm_start.le((start+end)/2) & arms.arm_end.gt((start+end)/2)]
        assert len(arm_hits) == 1 and end>start
        loci.append(dict(gene=gene, ensembl_gene_id=ens, chromosome=chrom,start=start,end=end,arm=arm_hits.iloc[0].arm,
            coordinate_convention="0-based half-open", assembly="GRCh38",source_url=url,annotation_version=version,
            panel_group="oncogene_gain_screen" if gene in GAIN_GENES else "other_locus_screen"))
    loci=pd.DataFrame(loci)
    save(loci,"gene_loci.csv")
    models=read("metadata/resource_models.csv")
    rna=read("output/rna_tpm.csv").set_index("gene_id")
    prot=read("output/prot_abundance_matrix.csv").set_index("protein")
    pq=read("output/prot_qc.csv")
    segs=read("output/wes_cnv_segments.csv")
    cnr_map=read("output/wes_recovered_provenance_cnv_support.csv")
    manifest=read("output/wes_cnv_target_only/manifest.csv")
    for row in manifest.itertuples():
        assert sha(ROOT/row.cns_path)==row.cns_sha256
    records=[]; evidence=[]
    for sample in cnr_map.itertuples():
        f=source(sample.cnr_source)
        assert sha(f)==sample.cnr_sha256
        INPUTS[sample.cnr_source]=sample.cnr_sha256
        bins=pd.read_csv(f,sep="\t")
        bins=bins.loc[bins.gene.ne("Antitarget")].copy()
        s=segs.loc[segs.cell_line.eq(sample.cell_line)]
        offsets=s.log2_raw-s.log2c_auto
        assert offsets.max()-offsets.min()<1e-10
        center=float(offsets.iloc[0])
        bins["log2c_auto"]=bins.log2-center
        for gene in loci.itertuples():
            hit=s.loc[s.chromosome.eq(gene.chromosome) & s.start.lt(gene.end) & s.end.gt(gene.start)]
            ov=np.minimum(hit.end,gene.end)-np.maximum(hit.start,gene.start)
            bp=float(ov.sum())
            segment=float(np.average(hit.log2c_auto,weights=ov)) if bp else np.nan
            b=bins.loc[bins.chromosome.eq(gene.chromosome) & bins.start.lt(gene.end) & bins.end.gt(gene.start)].copy()
            b["cell_line"]=sample.cell_line; b["gene_symbol"]=gene.gene
            b["target_bin_width"]=b.end-b.start
            b["gene_overlap_bp"]=np.minimum(b.end,gene.end)-np.maximum(b.start,gene.start)
            evidence.append(b)
            usable=b.loc[b.depth.gt(0) & b.weight.gt(0)]
            n=len(usable); med=usable.log2c_auto.median(); minv=usable.log2c_auto.min(); maxv=usable.log2c_auto.max()
            arm=arms.loc[arms.arm.eq(gene.arm)].iloc[0]
            arm_seg=s.loc[s.chromosome.eq(gene.chromosome) & s.start.lt(arm.arm_end) & s.end.gt(arm.arm_start)]
            aw=np.minimum(arm_seg.end,arm.arm_end)-np.maximum(arm_seg.start,arm.arm_start)
            am=wmedian(arm_seg.log2c_auto,aw)
            full=bp==gene.end-gene.start
            # >=3 positive bins and consistent direction between measured bins
            # and segment; these explicit screens prioritise review only.
            eligible=bool(full and n>=3)
            gain=bool(eligible and segment>=1 and med>=1)
            loss=bool(eligible and segment<=-1 and med<=-1)
            records.append(dict(cell_line=sample.cell_line,gene=gene.gene,arm=gene.arm,
                segment_gene_covered_fraction=bp/(gene.end-gene.start),gene_segment_log2=segment,
                gene_segment_min_log2=hit.log2c_auto.min(),gene_segment_max_log2=hit.log2c_auto.max(),
                n_segments_overlapping_gene=len(hit),target_bins=len(b),positive_target_bins=n,
                zero_depth_target_bins=int(b.depth.eq(0).sum()),target_overlap_bp=int(b.gene_overlap_bp.sum()),
                target_bin_median_log2=med,target_bin_min_log2=minv,target_bin_max_log2=maxv,
                target_bin_iqr_log2=usable.log2c_auto.quantile(.75)-usable.log2c_auto.quantile(.25),
                gene_target_length_weighted_mean_depth=float(np.average(b.depth,weights=b.target_bin_width)) if len(b) else np.nan,
                arm_median_log2=am,segment_minus_arm_log2=segment-am,target_median_minus_arm_log2=med-am,
                segment_minus_target_median_log2=segment-med,
                relative_screen_eligible=eligible,
                supported_relative_gain_screen=gain,supported_relative_loss_screen=loss,
                segment_bin_disagreement_flag=bool(np.isfinite(segment) and np.isfinite(med) and abs(segment-med)>=.5),
                target_only_extreme_flag=bool(n>=3 and np.isfinite(med) and abs(med)>=1 and not (gain or loss)),
                zero_depth_followup_flag=bool(len(b)>=3 and b.depth.eq(0).mean()>=.5),
                status="relative gene screen; absolute copy number, allele-specific loss and clinical classification unestablished"))
    cn=pd.DataFrame(records)
    assert len(cn)==23*len(GENES) and not cn.duplicated(["cell_line","gene"]).any()
    allrows=models[["cell_line","patient_id","histotype_code","rna_passage","wes_passage","tmt_plex"]].merge(loci[["gene","ensembl_gene_id","panel_group"]],how="cross")
    allrows=allrows.merge(cn,how="left",on=["cell_line","gene"],validate="one_to_one")
    for gene, ix in allrows.groupby("gene").groups.items():
        vals=rna.loc[ids[gene]]
        allrows.loc[ix,"rna_tpm"]=allrows.loc[ix,"cell_line"].map(vals)
        row=pq.loc[pq.symbol.eq(gene) & pq.symbol_representative]
        assert len(row)<=1
        if len(row):
            allrows.loc[ix,"protein_log2_normalised"]=allrows.loc[ix,"cell_line"].map(prot.loc[row.iloc[0]["row"]])
            allrows.loc[ix,"protein_isodoping"]=bool(row.iloc[0].isoDoping)
        # Ranks are descriptive; patient and histotype denominators remain in the table.
        for field in ["rna_tpm","protein_log2_normalised"]:
            allrows.loc[ix,field+"_rank_desc"]=allrows.loc[ix,field].rank(method="min",ascending=False)
            allrows.loc[ix,field+"_n_models"]=allrows.loc[ix,field].count()
    allrows["wes_measured"]=allrows.positive_target_bins.notna()
    allrows["status"]=allrows.status.fillna("not assessed: WES unavailable")
    save(allrows,"gene_model_evidence.csv")
    save(pd.concat(evidence,ignore_index=True),"gene_target_bins.csv")
    events=allrows.loc[allrows.supported_relative_gain_screen.eq(True) | allrows.supported_relative_loss_screen.eq(True)].copy()
    assert events.relative_screen_eligible.eq(True).all()
    save(events,"supported_relative_events.csv")
    save(allrows.loc[allrows.segment_bin_disagreement_flag.eq(True) | allrows.target_only_extreme_flag.eq(True) | allrows.zero_depth_followup_flag.eq(True)],"locus_review_queue.csv")
    frequencies=[]
    for (gene,hist),d in allrows.groupby(["gene","histotype_code"]):
        w=d.loc[d.wes_measured]
        eligible=w.loc[w.relative_screen_eligible.eq(True)]
        for event,col in [("gain","supported_relative_gain_screen"),("loss","supported_relative_loss_screen")]:
            positive=w.loc[w[col].eq(True)]
            frequencies.append(dict(gene=gene,histotype=hist,event=event,
                models_with_WES=len(w),patients_with_WES=w.patient_id.nunique(),
                models_screen_eligible=len(eligible),patients_any_model_screen_eligible=eligible.patient_id.nunique(),
                models_screen_positive=len(positive),patients_any_model_screen_positive=positive.patient_id.nunique()))
    save(pd.DataFrame(frequencies),"patient_event_counts.csv")
    # Related models: report within-patient contrasts rather than independent cases.
    family=[]
    for (patient,gene),d in allrows.loc[allrows.wes_measured].groupby(["patient_id","gene"]):
        if len(d)<2:continue
        family.append(dict(patient_id=patient,gene=gene,n_wes_models=len(d),models=";".join(d.cell_line),
            segment_log2_range=d.gene_segment_log2.max()-d.gene_segment_log2.min(),
            target_median_log2_range=d.target_bin_median_log2.max()-d.target_bin_median_log2.min(),
            n_supported_gains=int(d.supported_relative_gain_screen.eq(True).sum()),
            n_supported_losses=int(d.supported_relative_loss_screen.eq(True).sum())))
    save(pd.DataFrame(family),"within_patient_contrasts.csv")
    # Prior two-locus gene medians and segmentation must be reproduced exactly.
    prior=read("reports/clinical_classification_2026-09-06/cnv_locus_summary.csv")
    compare=cn.merge(prior,on=["cell_line","gene"],validate="one_to_one")
    assert len(compare)==46
    assert np.allclose(compare.gene_segment_log2,compare.gene_segment_log2c_auto,atol=1e-12)
    assert np.allclose(compare.target_bin_median_log2,compare.gene_target_median_log2c_auto,atol=1e-12)
    summary=dict(script_sha256=sha(__file__),input_sha256=INPUTS,reference_sha256=sha(raw_path),models=42,patients=34,genes=len(GENES),
        model_gene_rows=len(allrows),wes_model_gene_rows=len(cn),screen_eligible_model_gene_rows=int(cn.relative_screen_eligible.sum()),
        positive_bin_supported_gain_rows=int(cn.supported_relative_gain_screen.sum()),
        positive_bin_supported_loss_rows=int(cn.supported_relative_loss_screen.sum()),prior_46_locus_results_reproduced=True,
        definitions=">=3 positive weighted target bins; full gene span in segments; both segment and target median >=1 (gain) or <=-1 (loss); research screens only",
        reference_note="52-locus curated panel, not genome-wide event discovery; new gene coordinates do not alter RNA quantification",
        canonical_outputs_modified=False)
    (DEST/"validation.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k!="input_sha256"},indent=2))


if __name__=="__main__":
    main()
