#!/usr/bin/env python3
"""Validate the received cluster return and prepare small, explicit pilot inputs.

Standard library only. No cluster access, analysis jobs or changes to source data.
Paths are observations reported by the cluster agent, not workstation observations.
"""
import csv
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
DOC = ROOT / "docs/cluster/molecular_extension_2026-09-06"
REC = ROOT / "docs/cluster/recovery/2026-09-06-molecular_extension"


def read(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    evidence_manifest = json.loads((REC / "received_files.json").read_text())
    for relative, expected in evidence_manifest["files"].items():
        assert digest(ROOT / relative) == expected, f"Received evidence changed: {relative}"
    inventory = read(REC / "input_inventory.tsv")
    models = read(DOC / "models.tsv")
    checks = {r["cell_line"]: r for r in read(REC / "target_only_cnr_check.tsv")}
    assert len(inventory) == 204 and len(models) == len(checks) == 23
    assert len({m["patient_id"] for m in models}) == 16
    assert set(checks) == {m["cell_line"] for m in models}
    roles = {
        "complete_md_cram": "duplicate-marked genome-wide CRAM",
        "published_recal_cram": "published recalibrated CRAM",
        "cnv_input_recal_cram": "Oct-24 recalibrated CRAM in work/",
        "cnv_input_bam": "manual-CNVkit input BAM",
        "target_only_cnr": "target-only CNR reconstructed",
        "corrected_cns": "corrected target-only CNS",
    }
    execution = []
    for model in models:
        line = model["cell_line"]
        row = {k: model[k] for k in ["cell_line", "patient_id", "wes_passage", "cnv_sample_id"]}
        for name, prefix in roles.items():
            matches = [r for r in inventory if r["cell_line"] == line and r["input_role"].startswith(prefix)]
            assert len(matches) == 1 and matches[0]["status"] == "verified", (line, prefix)
            item = matches[0]
            row[name] = item["verified_path"]
            if item["passage"]:
                assert item["passage"] == model["wes_passage"], (line, "passage")
            if name not in {"corrected_cns"}:
                assert model["cnv_sample_id"] in item["verified_path"], (line, name)
            if name == "target_only_cnr":
                assert item["sha256_if_small"] == model["target_only_cnr_sha256"]
            if name == "corrected_cns":
                assert item["verified_path"] == "repo:" + model["corrected_cns"]
                assert digest(ROOT / model["corrected_cns"]) == item["sha256_if_small"] == model["corrected_cns_sha256"]
        check = checks[line]
        assert check["sample"] == model["cnv_sample_id"]
        assert check["matches_repo_source_cnr"] == check["matches_repo_target_only"] == "True"
        assert check["cluster_cnr_sha256"] == model["source_cnr_sha256"]
        assert check["reconstructed_target_only_sha256"] == model["target_only_cnr_sha256"]
        assert int(check["n_rows"]) == 204706 and int(check["n_antitarget_dropped"]) == 40923
        row.update(target_only_cnr_sha256=model["target_only_cnr_sha256"],
                   corrected_cns_sha256=model["corrected_cns_sha256"],
                   msi_pilot=str(line in {"TOV21G", "TOV2414", "TOV3392D", "OV3331"}).lower(),
                   ascn_pilot=str(line in {"TOV81D", "OV2085", "TOV2929D", "OV3331", "OV1369-R2"}).lower(),
                   observation="cluster return 2026-09-06; recheck readability/index/RG before execution")
        execution.append(row)
    write_tsv(DOC / "execution_models.tsv", execution)
    normals = [{"reference_id": r["cell_line"], "source_bam": r["verified_path"],
                "source_bytes": r["bytes"], "processing": "bwa-mem2; not duplicate-marked; no BQSR",
                "role": "capture-matched unmatched reference; optional mapping-bias sensitivity after duplicate marking"}
               for r in inventory if r["input_role"].startswith("public CNV-reference exome sorted BAM")]
    assert len(normals) == 5
    write_tsv(DOC / "execution_normals.tsv", normals)

    variants = read(ROOT / "reports/molecular_extension_2026-09-06/variants/variant_read_evidence.csv", ",")
    flagged = read(ROOT / "reports/molecular_extension_2026-09-06/variants/bam_review_candidates.csv", ",")
    selected = {r["candidate_id"] for r in flagged}
    # The original 25-row BED omits the five strong CDK12 loss candidates.
    selected.update(r["candidate_id"] for r in variants if r["gene"] == "CDK12" and
                    r["source_variant_class"] in {"Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"})
    requests = []
    model_by_name = {m["cell_line"]: m for m in models}
    original_ids = {r["candidate_id"] for r in flagged}
    for v in variants:
        if v["candidate_id"] not in selected:
            continue
        pos = int(v["vcf_pos_1based"])
        requests.append({k: v[k] for k in ["candidate_id", "cell_line", "patient_id", "wes_passage", "gene",
            "chromosome", "vcf_pos_1based", "vcf_ref", "vcf_alt", "source_variant_class", "source_protein_label"]} |
            {"cnv_sample_id": model_by_name[v["cell_line"]]["cnv_sample_id"],
             "review_start_0based": max(0, pos - 1 - 100),
             "review_end_exclusive": pos - 1 + len(v["vcf_ref"]) + 100,
             "reason": v["bam_review_priority"] if v["candidate_id"] in original_ids else "CDK12_loss_candidate_completion"})
    assert len(requests) == len(selected) == 30
    write_tsv(DOC / "execution_variant_requests.tsv", requests)
    with (DOC / "execution_variant_regions.bed").open("w") as handle:
        for r in requests:
            handle.write("\t".join(map(str, [r["chromosome"], r["review_start_0based"],
                r["review_end_exclusive"], r["cell_line"] + "_" + r["gene"], r["reason"]])) + "\n")
    with (DOC / "execution_locus_regions.bed").open("w") as handle:
        loci = read(DOC / "locus_review_requests.tsv")
        for r in loci:
            handle.write("\t".join([r["chromosome"], r["review_start_0based"], r["review_end_exclusive"],
                r["cell_line"] + "_" + r["gene"], r["task_id"]]) + "\n")
    assert len(loci) == 25
    for kind, field in [("msi", "msi_pilot"), ("ascn", "ascn_pilot")]:
        (DOC / f"execution_{kind}_pilot_models.txt").write_text("".join(
            r["cnv_sample_id"] + "\n" for r in execution if r[field] == "true"))
    sources = [REC / "received_files.json", REC / "input_inventory.tsv", REC / "target_only_cnr_check.tsv",
               DOC / "models.tsv", DOC / "locus_review_requests.tsv",
               ROOT / "reports/molecular_extension_2026-09-06/variants/variant_read_evidence.csv",
               ROOT / "reports/molecular_extension_2026-09-06/variants/bam_review_candidates.csv"]
    result = dict(received_files_verified=len(evidence_manifest["files"]), inventory_rows=len(inventory),
                  wes_models=23, patients=16, cnr_manifest_matches=23, cns_local_hash_matches=23,
                  reference_normals=5, variant_requests=30, variant_models=len({r["cell_line"] for r in requests}),
                  locus_requests=25, msi_pilot_models=4, ascn_pilot_models=5,
                  cluster_execution=False, source_hashes={str(p.relative_to(ROOT)): digest(p) for p in sources},
                  output_hashes={str(p.relative_to(ROOT)): digest(p) for p in sorted(DOC.glob("execution_*"))
                                 if p.suffix in {".tsv", ".bed", ".txt"}})
    (DOC / "execution_input_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if not k.endswith("hashes")}, indent=2))


if __name__ == "__main__":
    main()
