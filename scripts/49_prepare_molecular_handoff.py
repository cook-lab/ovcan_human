#!/usr/bin/env python3
"""Prepare small model/locus checklists and an overlay for a public repo clone.

No cluster access, external annotation, job submission or raw-data packaging.
Run after the molecular-extension report modules have completed.
"""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

ROOT = Path(os.environ.get("OVCAN_PROJ", Path.cwd())).resolve()
DOC = ROOT / "docs/cluster/molecular_extension_2026-09-06"
REPORT = ROOT / "reports/molecular_extension_2026-09-06"


def read(path):
    with (ROOT / path).open(newline="") as f:
        return list(csv.DictReader(f))


def tsv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    DOC.mkdir(parents=True, exist_ok=True)
    models = {r["cell_line"]: r for r in read("metadata/resource_models.csv")}
    hints = {r["cell_line"]: r for r in read("reports/audit_2026-09-05/wes_cluster_models.csv")}
    received_inventory = ROOT / "docs/cluster/recovery/2026-09-06-molecular_extension/input_inventory.tsv"
    observed = set()
    if received_inventory.exists():
        with received_inventory.open(newline="") as handle:
            observed = {r["cell_line"] for r in csv.DictReader(handle, delimiter="\t")
                        if r["input_role"].startswith("duplicate-marked genome-wide CRAM") and r["status"] == "verified"}
    rows = []
    for r in read("output/wes_cnv_target_only/manifest.csv"):
        m = models[r["cell_line"]]
        rows.append(dict(cell_line=r["cell_line"], patient_id=m["patient_id"], histotype=m["histotype_code"],
            wes_passage=m["wes_passage"], rna_passage=m["rna_passage"], cnv_sample_id=r["sample_id"],
            historical_bam_hint=hints[r["cell_line"]]["historical_bam_path"],
            current_alignment_status=("cluster-observed 2026-09-06; see execution_models.tsv; recheck before running"
                                      if r["cell_line"] in observed else "unverified: resolve on cluster"),source_cnr=r["source_cnr"],
            source_cnr_sha256=r["source_cnr_sha256"],target_only_cnr_sha256=r["target_cnr_sha256"],
            target_only_cnr_scratch_path=r["target_cnr_scratch_path"],corrected_cns=r["cns_path"],
            corrected_cns_sha256=r["cns_sha256"]))
    assert len(rows) == 23 and len({r["patient_id"] for r in rows}) == 16
    tsv(DOC/"models.tsv", rows, list(rows[0]))
    loci = {r["gene"]:r for r in read("reports/molecular_extension_2026-09-06/copy_number/gene_loci.csv")}
    evidence = {(r["cell_line"],r["gene"]):r for r in read("reports/molecular_extension_2026-09-06/copy_number/gene_model_evidence.csv")}
    requests = {
        "AKT2":("OV3331", "same-model DNA/RNA/protein concordance; resolve amplicon and absolute CN"),
        "CCNE1":("TOV2929D OV2085 OV3291 TOV3133D TOV2881EP TOV2835EP", "relative gain/bin-segment sensitivity; absolute CN and local boundaries"),
        "CCND1":("OV3133-R OV3133-R2 TOV3133G TOV3133D", "within-patient CCND1 versus CCNE1 patterns; three targets"),
        "NF1":("TOV2835EP TOV3121D TOV3121EP", "deep/partial exon depletion; residual reads, mapping and breakpoint evidence"),
        "CDKN2A":("OV1369-R2 TOV1369 TOV2414", "zero-depth/partial locus; distinguish p16 and p14ARF exons"),
        "ERBB2":("TOV3392D TOV2835EP", "expression versus focal-bin and broad-segment discrepancy"),
        "BRCA2":("TOV81D", "pathogenic allele plus retained wild-type historical evidence; present allele-specific status"),
        "CDK12":("OV3133-R OV3133-R2 TOV3133D TOV3133G OV3291", "strong loss-candidate read support; allele-specific status and splice evidence")
    }
    out = []
    for gene,(names,purpose) in requests.items():
        l=loci[gene]
        for model in names.split():
            e=evidence[(model,gene)]
            out.append(dict(task_id="MEX03" if gene not in {"BRCA2","CDK12"} else "MEX08",
                cell_line=model,patient_id=models[model]["patient_id"],gene=gene,genome_build="GRCh38",
                chromosome=l["chromosome"],gene_start_0based=l["start"],gene_end_exclusive=l["end"],
                review_start_0based=max(0,int(l["start"])-1000),review_end_exclusive=int(l["end"])+1000,
                purpose=purpose,gene_segment_log2=e["gene_segment_log2"],
                positive_bin_median_log2=e["target_bin_median_log2"],target_bins=e["target_bins"],
                zero_depth_target_bins=e["zero_depth_target_bins"],annotation_source=l["source_url"]))
    tsv(DOC/"locus_review_requests.tsv",out,list(out[0]))
    tsv(DOC/"input_inventory.template.tsv",[],["task_id","cell_line","input_role","status","verified_path",
        "symlink_target","bytes","sha256_if_small","sample_run_id","genome_reference","passage",
        "software_version","source_command_or_log","evidence_notes"])
    tsv(DOC/"task_status.template.tsv",[{"task_id":f"MEX{i:02d}","status":"not_searched"} for i in range(1,12)],
        ["task_id","status","models_assessed","evidence_files","missing_inputs","proposed_command_file",
         "estimated_resources","operator_authorization","qc_result","next_action"])
    if received_inventory.exists():
        subprocess.run([sys.executable, str(ROOT / "scripts/51_prepare_cluster_execution_inputs.py")], check=True)
    # The archive overlays an existing repository clone; existing public files
    # remain there, while all locally added report modules/docs/scripts travel.
    files=set()
    for base in [ROOT/"docs/cluster/molecular_extension_2026-09-06",REPORT,
                 ROOT/"reports/clinical_classification_2026-09-06"]:
        for p in base.rglob("*"):
            if not p.is_file() or p.is_symlink() or "cache" in p.relative_to(base).parts or p.name=="bundle_validation.json":
                continue
            if p.suffix.lower() in {".md",".csv",".tsv",".json",".bed",".py",".R".lower(),".xml",".txt"}:
                files.add(p)
    for p in (ROOT/"scripts").iterdir():
        if p.name[:2].isdigit() and 40<=int(p.name[:2])<=51 and p.is_file():
            files.add(p)
    if received_inventory.exists():
        recovered = ROOT / "docs/cluster/recovery/2026-09-06-molecular_extension"
        received_manifest = json.loads((recovered / "received_files.json").read_text())
        files.update(ROOT / relative for relative in received_manifest["files"])
        files.update([recovered / "received_files.json", recovered / "LOCAL_REVIEW.md",
                      recovered / "proposed_commands/README.md",
                      ROOT / "docs/cluster/recovery/2026-09-05/README.md",
                      ROOT / "docs/cluster/recovery/2026-09-06/FOLLOWUP.md"])
    files.update(ROOT/p for p in ["docs/cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md",
        "README.md","docs/PROJECT_STATUS.md","docs/REPRODUCIBILITY.md","AGENTS.md",".gitignore",
        "scripts/09_wes_hrd.R","scripts/22_wes_signature_refit.R","output/wes_hrd_feasibility.md"])
    for p in files:
        assert p.is_relative_to(ROOT) and not p.is_symlink() and p.stat().st_size<20_000_000,p
        assert "data/cluster_wes_retrieval" not in str(p.relative_to(ROOT)),p
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    dest=ROOT/"output/cluster_handoff"
    dest.mkdir(parents=True,exist_ok=True)
    target=dest/"ovcan_molecular_extension_tasks_2026-09-06.tar.gz"
    note=("OvCAN molecular-extension overlay, prepared 6 September 2026.\n\n"
          "Extract into a working clone of cook-lab/ovcan_human (base ec8c47b or later).\n"
          "Entry point: docs/cluster/molecular_extension_2026-09-06/EXECUTION_PLAN.md\n"
          "This overlay contains small reports, candidate loci, templates and scripts.\n"
          "It contains no raw alignments, sequence files or whole-exome VCFs.\n"
          "Model-linked processed findings remain research annotations. No jobs run automatically.\n"
          "Incoming proposed_commands are historical and unexecuted; revise them per RECIPES.md before use.\n"
          "The author-approved 53-allele NCBI query is complete. Reuse the saved ClinVar snapshot.\n"
          "PDF previews are available on the workstation and are omitted from this task overlay.\n"
          "molecular_handoff_sha256.json verifies every included repository file.\n")
    with tarfile.open(target,"w:gz") as tar:
        for p in sorted(files):
            tar.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
        for name,data in [("MOLECULAR_HANDOFF_README.txt",note.encode()),
                          ("molecular_handoff_sha256.json",(json.dumps(manifest,indent=2)+"\n").encode())]:
            info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644
            tar.addfile(info,io.BytesIO(data))
    with tarfile.open(target,"r:gz") as tar:
        for member in tar:
            assert member.isfile() and not Path(member.name).is_absolute() and ".." not in Path(member.name).parts
            if member.name in manifest:
                assert hashlib.sha256(tar.extractfile(member).read()).hexdigest()==manifest[member.name]
    (DOC/"bundle_validation.json").write_text(json.dumps(dict(
        archive=str(target.relative_to(ROOT)),archive_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        archive_bytes=target.stat().st_size,repository_files=len(files),models=len(rows),
        patients=len({r["patient_id"] for r in rows}),locus_requests=len(out),
        safe_members_and_hashes_verified=True,base_commit="ec8c47b",no_cluster_execution=True),indent=2)+"\n")
    print(json.dumps({"archive":str(target),"files":len(files),"bytes":target.stat().st_size,
                      "models":len(rows),"locus_requests":len(out)},indent=2))


if __name__=="__main__":
    main()
