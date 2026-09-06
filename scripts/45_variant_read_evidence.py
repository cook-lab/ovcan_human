#!/usr/bin/env python3
"""Read-only source-VCF evidence extension for retained molecular-panel candidates.

Uses Python's standard library. Writes only the requested report directory; does
not change canonical variants or tiers. Exact source-file hashes and the archived
MAF/VCF record-order linkage are checked, then allele edits are independently
compared after common-prefix/suffix trimming (not reference left-normalization).
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path

csv.field_size_limit(10_000_000)
ROOT = Path(os.environ.get("OVCAN_PROJ", Path(__file__).resolve().parents[1])).resolve()
CORE = "CCNE1 ERBB2 FOLR1 ESR1 PGR MSLN TACSTD2 CD274 BRCA1 BRCA2 PALB2 RAD51C RAD51D BRIP1 ATM ATR CHEK2 CDK12 MLH1 MSH2 MSH6 PMS2 POLE POLD1 ARID1A SMARCA4 SMARCA2 RB1 PTEN KRAS NRAS BRAF PIK3CA TP53 FBXW7 MYC".split()
KEYS = ["cell_line", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
INFO_KEEP = "TLOD NALOD NLOD GERMQ POPAF ROQ STRQ RU RPA ECNT MPOS MBQ MMQ MFRL AS_FilterStatus AS_SB_TABLE".split()


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read(path, delimiter=","):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def save(path, rows):
    assert rows, path
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def key(r, cell=None):
    return tuple((cell if k == "cell_line" and cell is not None else r[k]) for k in KEYS)


def edit(pos0, ref, alt):
    """Minimal substitution at a zero-based interbase position; no genome needed."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos0 = ref[1:], alt[1:], pos0 + 1
    return pos0, ref, alt


def maf_edit(r):
    ref, alt = r["Reference_Allele"].replace("-", ""), r["Tumor_Seq_Allele2"].replace("-", "")
    pos0 = int(r["Start_Position"]) if not ref else int(r["Start_Position"]) - 1
    return edit(pos0, ref, alt)


def numbers(s, cast=float):
    return [None if x in ("", ".", "NA") else cast(x) for x in s.split(",")] if s else []


def val(s, index, cast=float):
    n = numbers(s, cast)
    return n[index] if len(n) > index else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "reports/molecular_extension_2026-09-06/variants")
    parser.add_argument("--additional-genes", nargs="*", default=[])
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    panel = list(dict.fromkeys(CORE + args.additional_genes))
    inputs = {}

    def source(relative):
        p = ROOT / relative
        if relative.startswith("judy_archive/data/") and os.environ.get("OVCAN_DATA"):
            p = Path(os.environ["OVCAN_DATA"]) / Path(relative).relative_to("judy_archive/data")
        inputs[relative] = sha(p)
        return p

    canonical_paths = ["output/wes_mutations_filtered.csv", "output/wes_driver_tiers.csv", "metadata/resource_models.csv", "output/wes_input_manifest.csv", "output/wes_vcf_header_provenance.csv"]
    paths = {p: source(p) for p in canonical_paths}
    all_candidates = read(paths[canonical_paths[0]])
    selected = [r for r in all_candidates if r["Hugo_Symbol"] in panel]
    selected_by_key = {key(r): r for r in selected}
    assert len(selected) == len(selected_by_key)
    tiers = {key(r): r for r in read(paths[canonical_paths[1]])}
    models = {r["cell_line"]: r for r in read(paths[canonical_paths[2]])}
    manifests = read(paths[canonical_paths[3]])
    versions = {r["cell_line"]: r for r in read(paths[canonical_paths[4]])}
    assert len(all_candidates) == 6194 and len(manifests) == 23
    evidence, csq_rows, files, matched = [], [], [], set()
    total_records = 0
    for m in manifests:
        cell = m["cell_line"]
        vp, mp = source(m["source_vcf"]), source(m["source_maf"])
        assert inputs[m["source_vcf"]] == m["vcf_sha256"]
        assert inputs[m["source_maf"]] == m["maf_sha256"]
        opener = gzip.open if vp.suffix == ".gz" else open
        count, kept, csq_fields, sample = 0, 0, None, None
        with opener(vp, "rt") as vf, mp.open() as mf:
            for line in vf:
                if line.startswith("##INFO=<ID=CSQ,"):
                    csq_fields = line.split("Format: ", 1)[1].split('"', 1)[0].strip().split("|")
                if line.startswith("#CHROM"):
                    header = line.rstrip().split("\t")
                    assert len(header) == 10, (cell, "expected one tumour sample")
                    sample = header[9]
                    break
            assert csq_fields and sample
            mr = csv.DictReader((line for line in mf if not line.startswith("#")), delimiter="\t")
            for count, pair in enumerate(itertools.zip_longest(vf, mr), 1):
                line, r = pair
                assert line is not None and r is not None, (cell, "MAF/VCF row-count mismatch")
                v = line.rstrip().split("\t")
                assert v[0] == r["Chromosome"]
                k = key(r, cell)
                if k not in selected_by_key:
                    continue
                assert k not in matched
                candidate = selected_by_key[k]
                alts = v[4].split(",")
                hits = [i for i, alt in enumerate(alts) if edit(int(v[1])-1, v[3], alt) == maf_edit(r)]
                assert len(hits) == 1, (cell, k, "exact minimal allele edit must match uniquely", hits)
                ai = hits[0]
                info = dict((p.split("=", 1) if "=" in p else (p, "true")) for p in v[7].split(";"))
                fmt = dict(zip(v[8].split(":"), v[9].split(":")))
                ad = numbers(fmt.get("AD", ""), int)
                assert len(ad) == len(alts)+1 and all(x is not None for x in ad)
                refn, altn, adn = ad[0], ad[ai+1], sum(ad)
                af, tlod = val(fmt.get("AF", ""), ai), val(info.get("TLOD", ""), ai)
                assert af is None or 0 <= af <= 1
                sb = info.get("AS_SB_TABLE", "").split("|")
                strands = numbers(sb[ai+1], int) if len(sb) > ai+1 else []
                alt_fwd, alt_rev = strands if len(strands) == 2 else (None, None)
                flags = []
                if altn < 5: flags.append("alt_reads_lt5")
                if af is not None and af < .1: flags.append("caller_AF_lt0.10")
                if tlod is not None and tlod < 10: flags.append("TLOD_lt10")
                if altn >= 5 and alt_fwd is not None and min(alt_fwd, alt_rev) == 0: flags.append("alternate_support_one_strand")
                if info.get("STR") == "true": flags.append("tandem_repeat_context")
                if len(v[3]) != len(alts[ai]) and altn < 10: flags.append("indel_alt_reads_lt10")
                original = candidate.get("HGVSp_Short", "")
                row = dict(candidate_id="|".join(k), cell_line=cell, patient_id=models[cell]["patient_id"],
                    histotype=models[cell]["histotype_code"], wes_passage=models[cell]["wes_passage"], gene=r["Hugo_Symbol"],
                    panel_scope="core36" if r["Hugo_Symbol"] in CORE else "existing_driver_extension",
                    genome_build="GRCh38", maf_start_1based=r["Start_Position"], maf_end_1based=r["End_Position"],
                    maf_ref=r["Reference_Allele"], maf_alt=r["Tumor_Seq_Allele2"], source_protein_label=original,
                    source_protein_label_canonical_flag=candidate.get("hgvsp_canonical", ""),
                    source_protein_label_reconstructed_flag=candidate.get("hgvsp_reconstructed", ""),
                    source_transcript=r["Transcript_ID"], source_consequence=r["Consequence"], source_variant_class=r["Variant_Classification"],
                    source_existing_variation=r.get("Existing_variation", ""), source_CLIN_SIG=r.get("CLIN_SIG", ""),
                    original_research_tier=tiers.get(k, {}).get("tier", "not in selected driver table"),
                    original_tier_rationale=tiers.get(k, {}).get("rationale", ""),
                    chromosome=v[0], vcf_pos_1based=int(v[1]), vcf_ref=v[3], vcf_alt=alts[ai], vcf_alt_index_1based=ai+1,
                    vcf_number_alt_alleles=len(alts), vcf_FILTER=v[6], maf_FILTER=r["FILTER"],
                    minimal_edit_pos_0based=maf_edit(r)[0], minimal_edit_deleted=maf_edit(r)[1], minimal_edit_inserted=maf_edit(r)[2],
                    allele_match="exact minimal edit and archived MAF/VCF row linkage; not reference left-normalized",
                    source_vcf_sample=sample, vcf_GT=fmt.get("GT", ""), AD_raw=fmt.get("AD", ""),
                    AD_ref=refn, AD_alt=altn, AD_sum_all_alleles=adn, AD_alt_fraction_all=altn/adn if adn else None,
                    AD_alt_fraction_ref_plus_selected=altn/(refn+altn) if refn+altn else None,
                    caller_AF=af, FORMAT_DP=val(fmt.get("DP", ""), 0, int),
                    FAD_raw=fmt.get("FAD", ""), F1R2_ref=val(fmt.get("F1R2", ""), 0, int),
                    F1R2_alt=val(fmt.get("F1R2", ""), ai+1, int), F2R1_ref=val(fmt.get("F2R1", ""), 0, int),
                    F2R1_alt=val(fmt.get("F2R1", ""), ai+1, int), alt_forward=alt_fwd, alt_reverse=alt_rev,
                    phase_PGT=fmt.get("PGT", ""), phase_PID=fmt.get("PID", ""), phase_PS=fmt.get("PS", ""),
                    TLOD=tlod, **{"INFO_"+name: info.get(name, "") for name in INFO_KEEP if name != "TLOD"},
                    vcf_repeat_flag=info.get("STR") == "true", review_flags=";".join(flags),
                    nearby_selected_candidates="", same_phase_group_selected_candidates="", bam_review_priority="",
                    clinical_significance="not assigned by this extraction", variant_origin="unresolved by tumour-only calling",
                    biallelic_status="not established", current_stock_orthogonal_confirmation="not assessed",
                    source_vcf=m["source_vcf"], source_vcf_sha256=m["vcf_sha256"], source_maf=m["source_maf"], source_maf_sha256=m["maf_sha256"],
                    source_record_1based=count, source_record_sha256=hashlib.sha256(line.encode()).hexdigest(),
                    vep_version=versions[cell]["vep_version"], vep_assembly=versions[cell]["vep_assembly"],
                    source_clinvar_release=versions[cell]["clinvar_release"])
                evidence.append(row)
                for annotation in info.get("CSQ", "").split(","):
                    ann = dict(zip(csq_fields, annotation.split("|")))
                    if ann.get("SYMBOL") == r["Hugo_Symbol"]:
                        csq_rows.append(dict(candidate_id=row["candidate_id"], **ann))
                matched.add(k)
                kept += 1
        files.append(dict(cell_line=cell, source_vcf=m["source_vcf"], source_vcf_sha256=m["vcf_sha256"],
            source_vcf_hash_scope="compressed file bytes" if vp.suffix == ".gz" else "uncompressed file bytes",
            source_maf=m["source_maf"], source_maf_sha256=m["maf_sha256"], vcf_maf_records=count,
            selected_candidates=kept, source_hashes_match=True, row_counts_match=True))
        total_records += count
    assert total_records == 582474 and matched == set(selected_by_key)
    prior_path = source("reports/clinical_classification_2026-09-06/hrr_candidate_vcf_evidence.csv")
    prior_checks = 0
    for old in read(prior_path):
        found = [r for r in evidence if (r["cell_line"], r["chromosome"], r["vcf_pos_1based"], r["vcf_ref"], r["vcf_alt"]) ==
                 (old["cell_line"], old["chromosome"], int(old["vcf_pos"]), old["ref"], old["alt"])]
        assert len(found) == 1
        r = found[0]
        assert r["AD_raw"] == old["AD"] and abs(r["caller_AF"]-float(old["caller_AF"])) < 1e-10 and abs(r["TLOD"]-float(old["TLOD"])) < 1e-10
        prior_checks += 1
    for r in evidence:
        peers = [p for p in evidence if p is not r and p["cell_line"] == r["cell_line"] and p["chromosome"] == r["chromosome"]]
        near = [p["candidate_id"] for p in peers if abs(p["vcf_pos_1based"]-r["vcf_pos_1based"]) <= 10]
        phase = [p["candidate_id"] for p in peers if r["phase_PID"] not in ("", ".") and p["phase_PID"] == r["phase_PID"]]
        r["nearby_selected_candidates"] = ";".join(near)
        r["same_phase_group_selected_candidates"] = ";".join(phase)
        if near: r["review_flags"] += (";" if r["review_flags"] else "") + "nearby_selected_call_within10bp"
        if phase: r["review_flags"] += (";" if r["review_flags"] else "") + "shared_local_phase_group"
        flags = set(r["review_flags"].split(";"))
        substantive = flags - {"", "tandem_repeat_context"}
        r["bam_review_priority"] = "priority_1_complex_haplotype" if phase else "priority_2_read_support" if substantive else "priority_3_repeat_context_only" if flags != {""} else "no_heuristic_flag"
    evidence.sort(key=lambda r: (r["cell_line"], r["chromosome"], r["vcf_pos_1based"]))
    save(args.out / "variant_read_evidence.csv", evidence)
    save(args.out / "variant_transcript_consequences.csv", csq_rows)
    save(args.out / "source_manifest.csv", files)
    review = [r for r in evidence if r["review_flags"]]
    if review: save(args.out / "bam_review_candidates.csv", review)
    queue = sorted({(r["genome_build"], r["chromosome"], r["vcf_pos_1based"], r["vcf_ref"], r["vcf_alt"]) for r in evidence})
    with (args.out / "external_annotation_request.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["genome_build", "chromosome", "position_1based", "reference", "alternate"])
        writer.writerows(queue)
    with (args.out / "bam_review_regions.bed").open("w") as f:
        for r in review:
            f.write("\t".join(map(str, [r["chromosome"], max(0, r["vcf_pos_1based"]-101),
                r["vcf_pos_1based"]+len(r["vcf_ref"])+100, r["cell_line"]+"_"+r["gene"], r["bam_review_priority"]]))+"\n")
    compounds = []
    seen_phase = set()
    for r in evidence:
        if not r["same_phase_group_selected_candidates"]: continue
        pk = (r["cell_line"], r["chromosome"], r["phase_PID"])
        if pk in seen_phase: continue
        seen_phase.add(pk)
        group = [p for p in evidence if (p["cell_line"], p["chromosome"], p["phase_PID"]) == pk]
        compounds.append(dict(cell_line=r["cell_line"], gene=r["gene"], chromosome=r["chromosome"],
            source_PID=r["phase_PID"], source_PS=r["phase_PS"], source_PGTs=";".join(sorted({p["phase_PGT"] for p in group})),
            candidate_ids=";".join(p["candidate_id"] for p in group), record_count=len(group),
            net_length_change_if_same_haplotype=sum(len(p["minimal_edit_inserted"])-len(p["minimal_edit_deleted"]) for p in group),
            interpretation="Caller local-phase evidence only; BAM haplotype and combined transcript consequence require review; do not count as independent gene hits"))
    if compounds: save(args.out / "compound_haplotype_review.csv", compounds)
    assertions_path = args.out / "curated_assertions.csv"
    curated_count = 0
    if assertions_path.exists():
        inputs[str(assertions_path.relative_to(ROOT))] = sha(assertions_path)
        assertions = read(assertions_path)
        assertion_map = {(r["chromosome"], int(r["vcf_pos_1based"]), r["vcf_ref"], r["vcf_alt"]): r for r in assertions}
        assert len(assertions) == len(assertion_map), "Curation requires one combined evidence row per exact allele"
        assert set(assertion_map) <= {(r["chromosome"], r["vcf_pos_1based"], r["vcf_ref"], r["vcf_alt"]) for r in evidence}
        annotated = []
        for r in evidence:
            a = assertion_map.get((r["chromosome"], r["vcf_pos_1based"], r["vcf_ref"], r["vcf_alt"]))
            if a:
                assert a["gene"] == r["gene"]
                annotated.append({**r, **{"curation_"+k:v for k,v in a.items() if k not in ("chromosome", "vcf_pos_1based", "vcf_ref", "vcf_alt", "gene")}})
        save(args.out / "curated_model_evidence.csv", annotated)
        curated_count = len(assertions)
    summary = dict(script_sha256=sha(Path(__file__)), input_sha256=inputs, core_genes=CORE,
        additional_genes=args.additional_genes, retained_panel_candidates=len(evidence),
        core_candidates=sum(r["panel_scope"] == "core36" for r in evidence),
        unique_alleles=len({(r["chromosome"], r["vcf_pos_1based"], r["vcf_ref"], r["vcf_alt"]) for r in evidence}),
        source_models=len(files), source_records=total_records, exact_allele_matches=len(matched),
        source_hashes_all_match=True, source_row_counts_all_match=True,
        multiallelic_selected_records=sum(r["vcf_number_alt_alleles"] > 1 for r in evidence),
        reviewed_flagged_candidates=len(review), source_transcript_annotations=len(csq_rows),
        external_coordinate_queue_rows=len(queue), compound_phase_groups=len(compounds), curated_unique_alleles=curated_count,
        prior_independent_hrr_records_reconciled=prior_checks,
        all_selected_vcf_and_maf_PASS=all(r["vcf_FILTER"] == r["maf_FILTER"] == "PASS" for r in evidence),
        recovered_indel_AF_count=sum(len(r["vcf_ref"]) != len(r["vcf_alt"]) for r in evidence),
        review_priority_counts=dict(collections.Counter(r["bam_review_priority"] for r in evidence)),
        flag_counts=dict(collections.Counter(f for r in evidence for f in r["review_flags"].split(";") if f)),
        flags_are="BAM-review prioritisation heuristics, not validated filters or evidence of a false call",
        allele_representation="Prefix/suffix-trimmed exact event; no reference-based left alignment or HGVS normalization claimed",
        canonical_outputs_unchanged=all(sha(paths[p]) == inputs[p] for p in canonical_paths))
    assert summary["canonical_outputs_unchanged"]
    (args.out / "validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k not in ("input_sha256", "core_genes")}, indent=2))


if __name__ == "__main__":
    main()
