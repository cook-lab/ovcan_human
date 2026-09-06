#!/usr/bin/env python3
"""Audit recovered WES provenance without executing recovered code or changing calls.

Uses only Python's standard library. Full validation requires the ignored recovery
bundle and original archived inputs; the resulting CSV/JSON files are reviewable
without those bulk inputs. Hashes explicitly distinguish file bytes from gzip's
uncompressed stream. The recovery report is evidence to check, not an authority.
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

csv.field_size_limit(10_000_000)
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT / "data/cluster_wes_retrieval/2026-09-05/ovcan_human_wes_handoff_2026-09-05"
SUBPOPS = ("AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "SAS")


def read_table(path, delimiter=","):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def file_hash(path, uncompressed=False):
    opener = gzip.open if uncompressed and path.suffix == ".gz" else open
    digest = hashlib.sha256()
    with opener(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_table(path, rows):
    assert rows, path
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def bed_rows(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) >= 3 and fields[1].isdigit():
                assert int(fields[2]) > int(fields[1]) >= 0
                rows.append(fields)
    return rows


def union_bp(rows):
    intervals = collections.defaultdict(list)
    for row in rows:
        intervals[row[0]].append((int(row[1]), int(row[2])))
    total = 0
    for chromosome in intervals.values():
        end = -1
        for start, stop in sorted(chromosome):
            total += max(0, stop - max(start, end))
            end = max(end, stop)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    assert bundle.is_dir(), f"Restore the recovery bundle: {bundle}"
    out = ROOT / "output"
    report = ROOT / "reports/wes_completion_2026-09-05"
    out.mkdir(exist_ok=True)
    report.mkdir(parents=True, exist_ok=True)
    evidence = set()

    def ev(relative):
        path = bundle / relative
        assert path.is_file(), path
        evidence.add(path)
        return path

    inventory_path = ev("recovery/file_inventory.tsv")
    inventory = read_table(inventory_path, "\t")
    inventory_by_path = {r["copied_relative_path"].removeprefix("wes/"): r for r in inventory
                         if r["copied_relative_path"]}
    declared_hashes = {r["cell_line"]: r for r in read_table(ev("recovery/vcf_hash_check.tsv"), "\t")}
    inputs = read_table(ROOT / "output/wes_input_manifest.csv")
    baseline = {r["cell_line"]: r for r in read_table(ROOT / "reports/audit_2026-09-05/wes_cluster_models.csv")}
    assert len(inputs) == len(baseline) == len(declared_hashes) == 23
    assert len({r["patient_id"] for r in baseline.values()}) == 16
    samplesheet = read_table(ev("scratch_archive_extract/samplesheet.csv"))
    assert len(samplesheet) == 23
    samples = {r["sample"]: r for r in samplesheet}
    assert set(samples) == {r["cnv_sample_id"] for r in baseline.values()}
    assert all(r["patient"] == r["sample"] and r["status"] == "1" for r in samplesheet)

    trace_path = ev("sarek_3.5.1_run/pipeline_info/execution_trace_2025-10-23_15-35-13.txt")
    trace = read_table(trace_path, "\t")
    assert len(trace) == 38541 and all(r["status"] == "COMPLETED" and r["exit"] == "0" for r in trace)
    vep_tasks = {re.search(r"\((.*?)\)", r["name"])[1]: r for r in trace
                 if r["name"].split(" (")[0].endswith(":ENSEMBLVEP_VEP")}
    assert len(vep_tasks) == 23
    history = [line.split("\t", 6) for line in ev("sarek_3.5.1_run/launch/nextflow_history.tsv").read_text().splitlines()]
    success = [r for r in history if r[3] == "OK"]
    assert len(success) == 1 and success[0][2] == "angry_allen"
    params = json.loads(ev("sarek_3.5.1_run/pipeline_info/params_2025-10-23_15-35-28.json").read_text())
    ev("sarek_3.5.1_run/pipeline_info/nf_core_sarek_software_mqc_versions.yml")
    ev("sarek_3.5.1_run/launch/run_sarek.sh")
    ev("sarek_3.5.1_run/launch/obcf-graham.cfg")
    fai_path = ev("scratch_archive_extract/references/Homo_sapiens_assembly38.fasta.fai")
    reference_lengths = {r[0]: int(r[1]) for line in fai_path.read_text().splitlines()
                         for r in [line.split("\t")]}

    model_rows = []
    for record in inputs:
        cell = record["cell_line"]
        sample = baseline[cell]["cnv_sample_id"]
        vcf = ROOT / record["source_vcf"]
        maf = ROOT / record["source_maf"]
        file_sha = file_hash(vcf)
        stream_sha = file_hash(vcf, uncompressed=True)
        assert file_sha == record["vcf_sha256"]
        assert file_hash(maf) == record["maf_sha256"]
        opener = gzip.open if vcf.suffix == ".gz" else open
        with opener(vcf, "rt") as handle:
            lines = handle.readlines()
        header = next(line for line in lines if line.startswith("##VEP="))
        vcf_contigs = {match[1]: int(match[2]) for line in lines if line.startswith("##contig=")
                       for match in [re.search(r"ID=([^,]+),length=(\d+)", line)] if match}
        assert vcf_contigs and all(reference_lengths.get(chrom) == length for chrom, length in vcf_contigs.items())
        work_hash = re.search(r"/work/([^/]+/[^/]+)/113_GRCh38", header)[1]
        short_hash = work_hash.split("/")[0] + "/" + work_hash.split("/")[1][:6]
        assert short_hash == vep_tasks[sample]["hash"]
        assert 'time="2025-10-23 ' in header and '##VEP="v113.0"' in header
        variants = [line.rstrip().split("\t") for line in lines if not line.startswith("#")]
        with maf.open() as handle:
            maf_rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
        assert len(variants) == len(maf_rows)
        changes = collections.Counter()
        for variant, row in zip(variants, maf_rows):
            # The converter preserves record order; normalized indels can move coordinates.
            assert variant[0] == row["Chromosome"]
            changes[(variant[6], row["FILTER"])] += 1
            expected_common = any(float((row.get(f"gnomADe_{pop}_AF") or "0").split("/")[0]) > .0004
                                  for pop in SUBPOPS)
            assert ("common_variant" in row["FILTER"]) == expected_common
        other_changes = sum(count for (old, new), count in changes.items()
                            if old != new and new != ("" if old in ("PASS", ".") else old + ";") + "common_variant")
        assert other_changes == 0
        stats_path = ev(f"sarek_3.5.1_run/reports/bcftools/mutect2/{sample}/{sample}.mutect2.filtered.bcftools_stats.txt")
        before_vep = int(next(line.split("\t")[3] for line in stats_path.read_text().splitlines()
                             if line.startswith("SN\t0\tnumber of records:")))
        assert before_vep == len(variants)
        declared = declared_hashes[cell]["cluster_uncompressed_sha256"]
        scope = "uncompressed stream" if declared == stream_sha else "compressed file bytes" if declared == file_sha else "neither"
        assert scope != "neither"
        model_rows.append(dict(cell_line=cell, sample_id=sample, patient_id=baseline[cell]["patient_id"],
            source_vcf=record["source_vcf"], source_maf=record["source_maf"],
            vcf_file_sha256=file_sha, vcf_uncompressed_sha256=stream_sha,
            prior_manifest_hash_scope="compressed file bytes" if vcf.suffix == ".gz" else "uncompressed file bytes",
            recovery_declared_cluster_hash=declared, recovery_hash_matching_local_scope=scope,
            recovery_uncompressed_label_correct=declared == stream_sha,
            cluster_bytes_rehashed_locally=False,
            vep_task_hash_in_vcf=work_hash, vep_task_hash_in_successful_trace=vep_tasks[sample]["hash"],
            vep_task_trace_match=True, filtered_vcf_records_before_vep=before_vep,
            annotated_vcf_records=len(variants), maf_records=len(maf_rows),
            vcf_pass=sum(v[6] == "PASS" for v in variants),
            pass_retagged_common_variant=changes[("PASS", "common_variant")],
            maf_pass=sum(m["FILTER"] == "PASS" for m in maf_rows),
            other_filter_changes=other_changes, common_variant_subpop_cutoff=.0004))
    assert sum(r["annotated_vcf_records"] for r in model_rows) == 582474
    assert sum(r["vcf_pass"] for r in model_rows) == 19816
    assert sum(r["pass_retagged_common_variant"] for r in model_rows) == 3735
    assert sum(r["maf_pass"] for r in model_rows) == 16081
    assert sum(not r["recovery_uncompressed_label_correct"] for r in model_rows) == 1
    write_table(out / "wes_recovered_provenance_models.csv", model_rows)

    cnv_root = ROOT / "judy_archive/data/cnvkit wes - new"
    cnv_matches = []
    for recovered in sorted((bundle / "cnvkit_0.9.10_manual").glob("*_new/*.cns")):
        evidence.add(recovered)
        digest = file_hash(recovered)
        matches = sorted(cnv_root.rglob(recovered.name))
        assert matches and all(file_hash(f) == digest for f in matches)
        cnv_matches.append(dict(recovered_path=str(recovered.relative_to(bundle)),
            basename=recovered.name, sha256=digest, archived_copy_count=len(matches),
            all_archived_copies_identical=True, archived_paths=";".join(str(f.relative_to(ROOT)) for f in matches)))
    assert len(cnv_matches) == 69
    write_table(out / "wes_recovered_provenance_cnv_files.csv", cnv_matches)
    pooled_references = sorted(cnv_root.rglob("reference.cnn"))
    reference_sha = {file_hash(f) for f in pooled_references}
    assert len(pooled_references) == 35 and len(reference_sha) == 1
    pooled_sha = next(iter(reference_sha))
    write_table(out / "wes_recovered_provenance_cnv_reference.csv", [dict(
        role="empirical pooled CNV reference", archived_copies=len(pooled_references),
        bytes=pooled_references[0].stat().st_size, sha256=pooled_sha,
        source_path=str(pooled_references[0].relative_to(ROOT)),
        all_archived_copies_identical=True,
        per_normal_coverage_in_returned_bundle=False)])
    reference_rows = read_table(pooled_references[0], "\t")
    reference_kept_keys = set()
    reference_exclusions = []
    for group in ("target", "antitarget"):
        counts = collections.Counter()
        for row in reference_rows:
            if (row["gene"] == "Antitarget") != (group == "antitarget"):
                continue
            values = {k: float(row[k] or "nan") for k in ("log2", "depth", "spread", "gc")}
            flags = dict(log2_outside_pm5=abs(values["log2"]) > 5,
                         spread_gt1=values["spread"] > 1,
                         zero_depth=values["depth"] == 0,
                         gc_outside_0_3_to_0_7=values["gc"] < .3 or values["gc"] > .7)
            counts["total_bins"] += 1
            counts["excluded_union"] += any(flags.values())
            counts.update({k: int(v) for k, v in flags.items()})
            if not any(flags.values()):
                reference_kept_keys.add((row["chromosome"], row["start"], row["end"]))
        reference_exclusions.append(dict(bin_group=group, **counts,
            criteria_overlap=True, interpretation="Reconstructed default-rule union; exact retained coordinate set checked against every CNR"))
    write_table(out / "wes_recovered_provenance_cnv_reference_filters.csv", reference_exclusions)
    del reference_rows

    off_target = []
    for cell, baseline_row in baseline.items():
        sample = baseline_row["cnv_sample_id"]
        def cnv(suffix):
            path = sorted(cnv_root.rglob(sample + "." + suffix))[0]
            return path, read_table(path, "\t")
        anti_path, anti = cnv("antitargetcoverage.cnn")
        target_path, target = cnv("targetcoverage.cnn")
        cnr_path, cnr = cnv("cnr")
        cns_path, cns = cnv("cns")
        assert {(r["chromosome"], r["start"], r["end"]) for r in cnr} == reference_kept_keys
        n_positive = sum(float(r["depth"]) > 0 for r in cnr)
        n_probes = sum(int(r["probes"]) for r in cns)
        assert n_positive == n_probes
        off_target.append(dict(cell_line=cell, sample_id=sample, antitarget_bins=len(anti),
            original_target_bins=len(target), original_target_depth_zero=sum(float(r["depth"]) == 0 for r in target),
            target_bins_in_cnr=sum(r["gene"] != "Antitarget" for r in cnr),
            positive_target_bins_in_cnr=sum(r["gene"] != "Antitarget" and float(r["depth"]) > 0 for r in cnr),
            reference_rule_coordinate_set_matches_cnr=True,
            antitarget_depth_zero=sum(float(r["depth"]) == 0 for r in anti),
            antitarget_depth_positive=sum(float(r["depth"]) > 0 for r in anti),
            positive_antitargets_in_cnr=sum(r["gene"] == "Antitarget" and float(r["depth"]) > 0 for r in cnr),
            positive_depth_bins_in_cnr=n_positive, total_segment_probes=n_probes,
            antitarget_source=str(anti_path.relative_to(ROOT)), antitarget_sha256=file_hash(anti_path),
            cnr_source=str(cnr_path.relative_to(ROOT)), cnr_sha256=file_hash(cnr_path),
            target_source=str(target_path.relative_to(ROOT)), target_sha256=file_hash(target_path),
            cns_source=str(cns_path.relative_to(ROOT)), cns_sha256=file_hash(cns_path)))
    assert all(r["antitarget_depth_zero"] == 42705 and r["antitarget_bins"] == 42709 for r in off_target)
    write_table(out / "wes_recovered_provenance_cnv_support.csv", off_target)

    bed_files = ["SeqCap_EZ_Exome_v3.targets.bed", "hglft_genome_2576bc_8f0.bed",
                 "hglft_genome_2576bc_8f0.sorted.bed", "intervals_sorted.bed",
                 "intervals_sorted.target.bed", "intervals.antitarget.bed"]
    beds, bed_stats = {}, []
    for name in bed_files:
        path = ev("scratch_archive_extract/" + name)
        rows = beds[name] = bed_rows(path)
        bed_stats.append(dict(file=name, sha256=file_hash(path), intervals=len(rows),
            summed_interval_bp=sum(int(r[2])-int(r[1]) for r in rows), union_bp=union_bp(rows),
            coordinate_convention="0-based half-open", path=str(path.relative_to(ROOT))))
    original = {f"{r[0]}:{int(r[1])+1}-{r[2]}" for r in beds[bed_files[0]]}
    lifted = {r[3] for r in beds[bed_files[1]]}
    ignored_rows = bed_rows(ev("scratch_archive_extract/hglft_genome_ignored.txt"))
    ignored = {f"{r[0]}:{int(r[1])+1}-{r[2]}" for r in ignored_rows}
    assert len(original) == 242232 and len(lifted) == 242215
    assert original - lifted == ignored and len(ignored) == 17 and not lifted - original
    for name in bed_files[2:4]:
        assert collections.Counter(tuple(r[:3]) for r in beds[name]) == collections.Counter(tuple(r[:3]) for r in beds[bed_files[1]])
    assert file_hash(ev("scratch_archive_extract/intervals_sorted.target.bed")) == file_hash(out / "wes_cnvkit_target_intervals.bed")
    vendor_path = ev("mcgill_r004741_provenance/SeqCap_EZ_Exome_v3.targets.bed")
    vendor_rows = bed_rows(vendor_path)
    # The scratch copy adds UCSC chromosome prefixes; it is not byte-identical.
    assert [("chr" + r[0], r[1], r[2]) for r in vendor_rows] == [tuple(r[:3]) for r in beds[bed_files[0]]]
    bed_stats.insert(0, dict(file="provider_SeqCap_EZ_Exome_v3.targets.bed", sha256=file_hash(vendor_path),
        intervals=len(vendor_rows), summed_interval_bp=sum(int(r[2])-int(r[1]) for r in vendor_rows),
        union_bp=union_bp(vendor_rows), coordinate_convention="0-based half-open; bare chromosome names",
        path=str(vendor_path.relative_to(ROOT))))
    write_table(out / "wes_recovered_provenance_intervals.csv", bed_stats)

    fastp = read_table(ev("sarek_3.5.1_run/multiqc/multiqc_data/multiqc_fastp.txt"), "\t")
    assert len(fastp) == 23
    for row in fastp:
        assert "--disable_adapter_trimming" in row["command"] and "--length_required 15" in row["command"]
        assert "--disable_quality_filtering" not in row["command"] and "--disable_length_filtering" not in row["command"]
        summary = ast.literal_eval(row["summary"])
        assert summary["fastp_version"] == "0.23.4"
        assert summary["sequencing"] == "paired end (100 cycles + 100 cycles)"
    # Pin the exact selected executed commands, without sourcing or executing them.
    for task in sorted((bundle / "sarek_3.5.1_run/task_evidence").iterdir()):
        for name in (".command.sh", ".exitcode", "versions.yml"):
            if (task / name).exists():
                evidence.add(task / name)
        assert (task / ".exitcode").read_text().strip() == "0"
    for name in ("commands.txt", "commands1.txt", "run_cnvkit.sh"):
        ev("cnvkit_0.9.10_manual/" + name)
    for path in (bundle / "cnvkit_0.9.10_manual/normal_samples").iterdir():
        if path.is_file():
            evidence.add(path)
    for name in ("Homo_sapiens_assembly38.dict", "Homo_sapiens_assembly38.fasta.fai"):
        ev("scratch_archive_extract/references/" + name)
    fingerprints = []
    for path in sorted(evidence):
        relative = str(path.relative_to(bundle))
        digest = file_hash(path)
        old = inventory_by_path.get(relative)
        match = digest == old["checksum"] if old and old["checksum_algorithm"] == "sha256" else None
        assert match is not False, f"Recovery evidence hash mismatch: {relative}"
        fingerprints.append(dict(bundle_path=relative, bytes=path.stat().st_size, sha256=digest,
            recovery_inventory_record=old["record_id"] if old else "",
            recovery_inventory_hash_match=match if match is not None else "not listed"))
    write_table(out / "wes_recovered_provenance_evidence.csv", fingerprints)

    result = dict(
        bundle=str(bundle.relative_to(ROOT)), evidence_files=len(fingerprints),
        models=23, patients=16, successful_run=success[0][2], run_session=success[0][5],
        run_start_history=success[0][0], history_duration=success[0][1],
        workflow="nf-core/sarek", workflow_declared_release="3.5.1", nextflow_version="24.10.2",
        workflow_history_script_hash=success[0][4], verified_git_commit=None,
        successful_trace_tasks=len(trace), all_trace_tasks_exit_zero=True,
        annotated_vcf_headers_linked_to_successful_vep_tasks=23,
        archived_vcf_file_hashes_match_manifest=23,
        recovery_cluster_hashes_match_local_uncompressed_stream=22,
        recovery_hash_scope_exception="TOV3121D: declared uncompressed hash equals compressed local file bytes",
        recovered_cnv_files_match_archived_sources=len(cnv_matches),
        pooled_cnv_reference_archived_copies=35, pooled_cnv_reference_sha256=pooled_sha,
        reference_fai_contigs=len(reference_lengths), all_vcf_contig_lengths_match_reference_fai=True,
        annotated_records=582474, caller_pass_records=19816,
        converter_common_variant_retags_among_pass=3735, maf_pass_records=16081,
        converter_common_variant_rule="Any supplied gnomADe AFR/AMR/ASJ/EAS/FIN/NFE/SAS AF >0.0004",
        converter_rule_verified_rows=582474,
        original_22_maf_converter_version=None, recovered_tov3121d_converter_version="1.6.22",
        vep_actual_version="113.0", vep_cache_version="113", vep_params_declared_version=params["vep_version"],
        vep_filter_common_configured="exclude co-located 1KG_ALL AF>0.01; no record-count reduction observed",
        pre_vs_post_vep_record_count_matches=23,
        fastp_adapter_trimming_disabled=True, fastp_quality_filtering_disabled=False,
        fastp_length_filtering_disabled=False, fastp_minimum_length=15,
        nominal_fastp_read_cycles="100+100; TOV81D mean read length is shorter in the provider-derived input",
        vendor_targets=242232, lifted_source_intervals=242215, unmapped_source_intervals=17,
        lifted_output_intervals=242421, calling_interval_sum_bp=63709951, calling_interval_union_bp=63514049,
        cnv_target_bins=290475, cnv_antitarget_bins=42709,
        cnv_target_bins_excluded_by_reference_rules=85769, cnv_target_bins_after_reference_filter=204706,
        zero_antitarget_bins_per_model=42705, positive_antitarget_bins_per_model=4,
        cnv_support_interpretation="Effectively target-supported relative CNV; segment spans interpolate between measured bins",
        verification_scope="Returned small evidence and original workstation VCF/MAF/CNV files; cluster alignments/VCF bytes not supplied",
        limitations=[
            "Cluster VCFs and alignments are not included; cluster VCF identity uses a reported checksum plus independent archived-header/trace linkage.",
            "Retained recalibrated-CRAM overwrite dates, quickcheck and identical-reexecution claims are recovery-reported, not independently verifiable from absent CRAM headers/bytes.",
            "Original FASTA/dataset/container bytes, exact workflow Git commit/local edits and liftOver chain/command version are not recovered here.",
            "The manual CNV command scripts and identical output files support provenance, but successful manual run logs and normal bwa-mem2 execution/version are not included.",
            "No matched normals establish variant somatic origin; relative target-derived CNV does not establish absolute ploidy, LOH, HRD or MSI."
        ])
    # No analysis/call tables or cohort metadata are modified by this audit.
    (out / "wes_recovered_provenance.json").write_text(json.dumps(result, indent=2) + "\n")
    (report / "provenance_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("models", "patients", "successful_run", "evidence_files", "annotated_records", "caller_pass_records", "maf_pass_records", "recovered_cnv_files_match_archived_sources")}, indent=2))


if __name__ == "__main__":
    main()
