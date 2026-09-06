#!/usr/bin/env python3
"""NCBI annotation of the explicitly approved 53-allele coordinate-only queue.

Normal runs replay saved responses offline. --fetch performs missing NCBI queries.
No model names, patient identifiers, read evidence or source paths are transmitted.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
BASE = ROOT/"reports/molecular_extension_2026-09-06/variants"
DEST = BASE/"clinvar"
QUEUE = BASE/"external_annotation_request.tsv"
VERSIONS = [11,12,12,12,10,12,14,11,12,11,10,12,11,9,10,10,11,10,10,11,9,11]
ACCESSIONS = {f"chr{i}":f"NC_{i:06d}.{v}" for i,v in enumerate(VERSIONS,1)}
ACCESSIONS.update(chrX="NC_000023.11",chrY="NC_000024.10")
HOSTS = {"api.ncbi.nlm.nih.gov","eutils.ncbi.nlm.nih.gov","ftp.ncbi.nlm.nih.gov"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save_csv(path, rows, fields=None):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction="raise")
        w.writeheader();w.writerows(rows)


def text(node, path):
    return node.findtext(path,default="") if node is not None else ""


def spdi(obj):
    return f'{obj["seq_id"]}:{obj["position"]}:{obj["deleted_sequence"]}:{obj["inserted_sequence"]}'


def preferred_traits(node):
    return sorted(set(e.text for e in node.findall('.//Name/ElementValue[@Type="Preferred"]') if e.text))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--fetch",action="store_true")
    args=ap.parse_args();DEST.mkdir(parents=True,exist_ok=True)
    responses=DEST/"responses";responses.mkdir(exist_ok=True)
    manifest_path=DEST/"requests.json"
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    with QUEUE.open(newline="") as f: queue=list(csv.DictReader(f,delimiter="\t"))
    assert len(queue)==53 and set(queue[0])=={"genome_build","chromosome","position_1based","reference","alternate"}
    assert len({tuple(q.values()) for q in queue})==53
    last_request=0.

    def get(url, kind):
        nonlocal last_request
        assert urllib.parse.urlparse(url).hostname in HOSTS
        key=sha(url.encode())
        path=responses/(key+"."+kind)
        if key in manifest:
            raw=path.read_bytes();assert sha(raw)==manifest[key]["sha256"]
            return raw,manifest[key]
        if not args.fetch:raise RuntimeError("Missing saved response; explicit --fetch required: "+url)
        for attempt in range(4):
            time.sleep(max(0,.42-(time.monotonic()-last_request)))
            last_request=time.monotonic()
            try:
                req=urllib.request.Request(url,headers={"User-Agent":"OvCAN-coordinate-annotation/1.0","Accept":"application/xml" if kind=="xml" else "application/json"})
                with urllib.request.urlopen(req,timeout=55) as response:
                    assert urllib.parse.urlparse(response.url).hostname in HOSTS
                    raw=response.read();status=response.status
                if kind=="json":json.loads(raw)
                if kind=="xml":ET.fromstring(raw)
                break
            except urllib.error.HTTPError as e:
                if e.code not in {429,500,502,503,504} or attempt==3:raise
                time.sleep(2**attempt)
        path.write_bytes(raw)
        item=dict(url=url,http_status=status,retrieved_utc=datetime.now(timezone.utc).isoformat(),
                  bytes=len(raw),sha256=sha(raw),response_file=str(path.relative_to(ROOT)))
        manifest[key]=item;manifest_path.write_text(json.dumps(manifest,indent=2)+"\n")
        return raw,item

    # Verify chromosome accession versions against NCBI's pinned GRCh38.p14 report.
    assembly_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
    raw,_=get(assembly_url,"txt")
    mapping={}
    for line in raw.decode().splitlines():
        if not line or line.startswith("#"):continue
        fields=line.split("\t")
        if fields[1]=="assembled-molecule":mapping[fields[9]]=fields[6]
    assert all(mapping[c]==a for c,a in ACCESSIONS.items())

    results=[];records={};submissions=[];conditions=[];matches=[]
    for index,q in enumerate(queue,1):
        assert q["genome_build"]=="GRCh38"
        input_spdi=f'{ACCESSIONS[q["chromosome"]]}:{int(q["position_1based"])-1}:{q["reference"]}:{q["alternate"]}'
        url="https://api.ncbi.nlm.nih.gov/variation/v0/spdi/"+urllib.parse.quote(input_spdi,safe=":")+"/canonical_representative"
        raw,norm_prov=get(url,"json");obj=json.loads(raw)
        assert "data" in obj and "error" not in obj,obj
        canonical=spdi(obj["data"])
        # Quoting the complete SPDI is essential: unquoted searches split colons.
        params={"db":"clinvar","term":f'"{canonical}"[Canonical SPDI]',"retmode":"json","retmax":100}
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"+urllib.parse.urlencode(params)
        raw,search_prov=get(url,"json");search=json.loads(raw)["esearchresult"]
        ids=search["idlist"];assert len(ids)==int(search["count"]),search
        assert not search.get("errorlist",{}).get("fieldsnotfound"),search
        exact=[];rejected=[]
        for uid in ids:
            if uid not in records:
                url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+urllib.parse.urlencode({"db":"clinvar","rettype":"vcv","id":f"VCV{int(uid):09d}"})
                raw,record_prov=get(url,"xml");tree=ET.fromstring(raw)
                archives=tree.findall("VariationArchive");assert len(archives)==1
                v=archives[0];assert v.attrib["VariationID"]==uid
                records[uid]=(v,record_prov)
            v,record_prov=records[uid]
            c=v.find("ClassifiedRecord");a=c.find("SimpleAllele") if c is not None else None
            # Never apply a haplotype/included-variant classification to one component.
            if a is None or text(a,"CanonicalSPDI")!=canonical or a.get("VariationID")!=uid:
                rejected.append(uid);continue
            exact.append(uid)
            classes=c.find("Classifications")
            summary=dict(**q,input_spdi=input_spdi,canonical_spdi=canonical,variation_id=uid,
                vcv=v.get("Accession")+"."+v.get("Version"),record_status=text(v,"RecordStatus"),
                record_last_updated=v.get("DateLastUpdated",""),name=v.get("VariationName",""),
                gene=";".join(g.get("Symbol","") for g in a.findall("GeneList/Gene")),
                allele_id=a.get("AlleleID",""),source_url=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/",
                record_response_sha256=record_prov["sha256"],retrieved_utc=record_prov["retrieved_utc"],
                identity_match="exact canonical SPDI; directly classified single allele",
                model_variant_origin="unresolved by tumour-only calling",model_biallelic_status="not established",
                current_stock_function="not assessed",ovarian_treatment_classification="not assigned")
            for prefix,tag in [("germline","GermlineClassification"),("somatic_clinical_impact","SomaticClinicalImpact"),("oncogenicity","OncogenicityClassification")]:
                x=classes.find(tag) if classes is not None else None
                summary[prefix+"_description"]=text(x,"Description")
                summary[prefix+"_review_status"]=text(x,"ReviewStatus")
                summary[prefix+"_last_evaluated"]=x.get("DateLastEvaluated","") if x is not None else ""
                summary[prefix+"_n_submissions"]=x.get("NumberOfSubmissions","") if x is not None else ""
                contributing=[];other=[]
                if x is not None:
                    for ts in x.findall("ConditionList/TraitSet"):
                        entry=dict(conditions=preferred_traits(ts),attributes=ts.attrib)
                        (contributing if ts.get("ContributesToAggregateClassification")=="true" else other).append(entry)
                summary[prefix+"_contributing_conditions"]=json.dumps(contributing,sort_keys=True)
                summary[prefix+"_other_conditions"]=json.dumps(other,sort_keys=True)
            # Keep exact framework-specific content for non-scalar impact tiers.
            summary["aggregate_classifications_xml"]=ET.tostring(classes,encoding="unicode") if classes is not None else ""
            mane=[]
            for h in a.findall("HGVSlist/HGVS"):
                n=h.find("NucleotideExpression")
                if n is not None and n.get("MANESelect")=="true":mane.append(text(n,"Expression"))
            summary["mane_select_hgvsc"]=";".join(mane)
            locations=[x.attrib for x in a.findall("Location/SequenceLocation") if x.get("Assembly")=="GRCh38"]
            assert locations,"Direct allele lacks GRCh38 coordinates"
            summary["grch38_locations_json"]=json.dumps(locations,sort_keys=True)
            summary["exact_grch38_vcf_tuple_match"]=any(
                x.get("Accession")==ACCESSIONS[q["chromosome"]] and x.get("positionVCF")==q["position_1based"]
                and x.get("referenceAlleleVCF")==q["reference"] and x.get("alternateAlleleVCF")==q["alternate"]
                for x in locations)
            matches.append(summary)
            for ca in c.findall("ClinicalAssertionList/ClinicalAssertion"):
                acc=ca.find("ClinVarAccession");cl=ca.find("Classification")
                if acc is None or cl is None:continue
                framework_fields=[("germline","GermlineClassification"),("somatic_clinical_impact","SomaticClinicalImpact"),("oncogenicity","OncogenicityClassification")]
                for framework,tag in framework_fields:
                    for x in cl.findall(tag):
                        submissions.append(dict(variation_id=uid,vcv=summary["vcv"],scv=acc.get("Accession","")+"."+acc.get("Version",""),
                            submitter=acc.get("SubmitterName",""),framework=framework,classification=x.text or "",
                            review_status=text(cl,"ReviewStatus"),last_evaluated=cl.get("DateLastEvaluated",""),
                            record_status=text(ca,"RecordStatus"),contributes_to_aggregate=ca.get("ContributesToAggregateClassification",""),
                            assertion_attributes=json.dumps(x.attrib,sort_keys=True),conditions=";".join(preferred_traits(ca.find("TraitSet"))) if ca.find("TraitSet") is not None else "",
                            condition_xml=ET.tostring(ca.find("TraitSet"),encoding="unicode") if ca.find("TraitSet") is not None else ""))
            for rcv in c.findall("RCVList/RCVAccession"):
                rc=rcv.find("RCVClassifications")
                if rc is None:continue
                conditions.append(dict(variation_id=uid,vcv=summary["vcv"],rcv=rcv.get("Accession","")+"."+rcv.get("Version",""),
                    title=rcv.get("Title",""),conditions=";".join(x.text or "" for x in rcv.findall("ClassifiedConditionList/ClassifiedCondition")),
                    classifications_xml=ET.tostring(rc,encoding="unicode")))
        status="exact_match" if exact else "no_exact_indexed_match" if not ids else "returned_record_not_direct_exact_allele"
        results.append(dict(**q,input_spdi=input_spdi,canonical_spdi=canonical,lookup_status=status,
            n_search_hits=len(ids),exact_variation_ids=";".join(exact),rejected_variation_ids=";".join(rejected),
            normalization_sha256=norm_prov["sha256"],search_sha256=search_prov["sha256"],
            search_query=search.get("querytranslation",""),search_warnings=json.dumps(search.get("warninglist",{}),sort_keys=True)))
        print(f"{index}/53: {status}; ClinVar IDs={';'.join(exact) or 'none'}",flush=True)
    assert len(results)==53
    save_csv(DEST/"allele_lookup.csv",results)
    save_csv(DEST/"matched_allele_annotations.csv",matches)
    save_csv(DEST/"submission_assertions.csv",submissions)
    save_csv(DEST/"condition_assertions.csv",conditions)
    local=list(csv.DictReader((BASE/"variant_read_evidence.csv").open()))
    key=lambda r:(r["chromosome"],int(r.get("position_1based",r.get("vcf_pos_1based"))),r.get("reference",r.get("vcf_ref")),r.get("alternate",r.get("vcf_alt")))
    look={key(r):r for r in results};ann={key(r):r for r in matches}
    assert len(ann)==len(matches),"Multiple directly classified VCVs per allele require manual reconciliation"
    joined=[]
    for r in local:
        k=key(r);out=dict(r)
        out["clinvar_local_haplotype_caveat"]=(
            "Any matched ClinVar classification describes the isolated allele; it is not transferred to the combined phased haplotype"
            if r.get("same_phase_group_selected_candidates","") else "")
        for field,val in look[k].items():out["clinvar_"+field]=val
        if k in ann:
            for field,val in ann[k].items():out["clinvar_"+field]=val
        joined.append(out)
    all_fields=list(dict.fromkeys(field for r in joined for field in r))
    save_csv(DEST/"model_annotations.csv",[{field:r.get(field,"") for field in all_fields} for r in joined],all_fields)
    report=dict(queue_rows=len(queue),queue_sha256=sha(QUEUE.read_bytes()),script_sha256=sha(Path(__file__).read_bytes()),
        model_evidence_sha256=sha((BASE/"variant_read_evidence.csv").read_bytes()),
        exact_matched_alleles=len(matches),no_exact_indexed_match=sum(r["lookup_status"]=="no_exact_indexed_match" for r in results),
        source_model_records=len(local),model_records_with_exact_match=sum(key(r) in ann for r in local),
        submission_assertion_rows=len(submissions),condition_assertion_rows=len(conditions),
        requests=len(manifest),all_response_hashes_checked=True,canonical_calls_modified=False,
        exact_grch38_vcf_tuple_matches=sum(r["exact_grch38_vcf_tuple_match"] for r in matches),
        user_authorization="Explicit 6 September 2026 approval: These data are not particularly sensitive—you can submit the query",
        transmitted_fields="reference accession/genomic position/ref/alt plus public ClinVar identifiers; no project model/read/path fields",
        completed_utc=datetime.now(timezone.utc).isoformat())
    (DEST/"validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
