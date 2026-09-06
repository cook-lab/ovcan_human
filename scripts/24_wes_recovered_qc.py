#!/usr/bin/env python3
"""Reconcile recovered WES QC reports without reading alignments or changing analyses.

Python standard library only. Run from any directory. The default evidence bundle
is ignored by Git; --source accepts the same archived internal layout elsewhere.
Outputs are new QC summaries, a source/checksum inventory and a validation report.
No external tools, network requests or pipeline commands are executed.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = PROJECT / "data/cluster_wes_retrieval/2026-09-05/ovcan_human_wes_handoff_2026-09-05"
NORMALS = ["SRR4039087", "SRR4039088", "SRR4039089", "SRR4039096", "SRR4039097"]
THRESHOLDS = (1, 10, 20, 30, 50, 100)
MOSDEPTH_SOURCE = "https://github.com/brentp/mosdepth/blob/v0.3.8/mosdepth.nim"


def table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def require(condition, description):
    if not condition:
        raise ValueError(description)


def near(left, right, tolerance=1e-8):
    return math.isclose(float(left), float(right), rel_tol=0, abs_tol=tolerance)


def samtools_stats(path):
    values = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] == "SN":
            name, text = fields[1].rstrip(":"), fields[2]
            values[name] = float(text) if any(x in text.lower() for x in (".", "e")) else int(text)
    require(values, f"No samtools summary numbers: {path}")
    return values


def picard_metrics(path):
    lines = path.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith("LIBRARY\t")]
    require(len(starts) == 1, f"Expected one Picard metric header: {path}")
    start = starts[0]
    rows = []
    for line in lines[start + 1:]:
        if not line or line.startswith("#"):
            break
        rows.append(dict(zip(lines[start].split("\t"), line.split("\t"))))
    require(len(rows) == 1, f"Expected one library per model: {path}")
    return rows[0]


def bed_accounting(path):
    intervals = defaultdict(list)
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        chrom, start, end = line.split("\t")[:3]
        start, end = int(start), int(end)
        require(0 <= start < end, f"Invalid BED interval: {line}")
        intervals[chrom].append((start, end))
    per_contig_sum = {chrom: sum(b - a for a, b in vals) for chrom, vals in intervals.items()}
    union_bp = 0
    for vals in intervals.values():
        left = right = None
        for start, end in sorted(vals):
            if left is None:
                left, right = start, end
            elif start <= right:
                right = max(right, end)
            else:
                union_bp += right - left
                left, right = start, end
        union_bp += right - left
    return dict(interval_count=sum(map(len, intervals.values())),
                summed_bp=sum(per_contig_sum.values()), union_bp=union_bp,
                overlap_counted_extra_bp=sum(per_contig_sum.values()) - union_bp,
                per_contig_summed_bp=per_contig_sum)


def distribution(path):
    total = {}
    cdf_sums = defaultdict(float)
    for line in path.read_text().splitlines():
        chrom, threshold, fraction = line.split("\t")
        # Reproduce the archived MultiQC per-contig values from its printed CDF
        # rows. Missing high-depth levels and rounded fractions make these a
        # different quantity from the authoritative integer depth-sum / length.
        cdf_sums[chrom] += float(fraction)
        if chrom == "total":
            require(re.fullmatch(r"\d+\.\d{2}", fraction), f"Unexpected distribution precision: {path}")
            require(int(threshold) not in total, f"Duplicate total threshold: {path}")
            total[int(threshold)] = float(fraction)
    require(all(t in total for t in range(151)), f"Missing 0–150x thresholds: {path}")
    require(total[0] == 1, f"Coverage at zero must equal one: {path}")
    ordered = [total[t] for t in sorted(total)]
    require(all(0 <= x <= 1 for x in ordered), f"Invalid coverage fraction: {path}")
    require(all(a >= b for a, b in zip(ordered, ordered[1:])), f"Nonmonotone distribution: {path}")
    return total, {chrom: value - 1 for chrom, value in cdf_sums.items()}


def write_csv(path, rows):
    require(rows, f"Refusing empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=PROJECT / "output")
    parser.add_argument("--report-dir", type=Path, default=PROJECT / "reports/wes_completion_2026-09-05")
    args = parser.parse_args()
    source = args.source.resolve()
    run = source / "sarek_3.5.1_run"
    mqc = run / "multiqc/multiqc_data"
    csv.field_size_limit(100_000_000)
    inventory = {row["copied_relative_path"].removeprefix("wes/"): row
                 for row in table(source / "recovery/file_inventory.tsv", "\t") if row["copied_relative_path"]}
    sources, source_ids, checks, metric_defs, values, metric_sources = [], {}, [], {}, defaultdict(list), defaultdict(set)

    def register(path, role, model="", absent_ok=False):
        path = Path(path)
        key = str(path.resolve())
        if key in source_ids:
            return source_ids[key]
        relative = path.relative_to(source).as_posix() if path.is_relative_to(source) else ""
        original = inventory.get(relative, {})
        exists = path.is_file()
        require(exists or absent_ok, f"Required source absent: {path}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if exists else ""
        declared = original.get("checksum", "") if original.get("checksum_algorithm") == "sha256" else ""
        if exists and declared:
            require(digest == declared, f"Recovered inventory hash mismatch: {relative}")
        sid = f"WQC{len(sources) + 1:04d}"
        source_ids[key] = sid
        try:
            project_path = path.relative_to(PROJECT).as_posix()
        except ValueError:
            project_path = str(path)
        sources.append(dict(source_id=sid, role=role, model_or_reference=model,
                            project_path=project_path, bundle_relative_path=relative,
                            sha256=digest, bytes=path.stat().st_size if exists else "",
                            original_cluster_path=original.get("original_path", ""),
                            inventory_sha256=declared,
                            inventory_hash_match=(digest == declared) if exists and declared else "",
                            availability="present" if exists else "not_in_local_bundle"))
        return sid

    def add(row, column, value, unit, denominator, stage, source_id, notes=""):
        row[column] = value
        definition = dict(unit=unit, denominator=denominator, alignment_stage=stage, notes=notes)
        require(column not in metric_defs or metric_defs[column] == definition, f"Inconsistent metric definition: {column}")
        metric_defs[column] = definition
        values[column].append((row["cell_line"], value))
        metric_sources[column].add(source_id)

    # Independently reconcile canonical model, patient and passage aliases.
    family_path = PROJECT / "metadata/line_family_map.csv"
    alias_path = PROJECT / "reports/audit_2026-09-05/wes_cluster_models.csv"
    manifest_path = PROJECT / "output/wes_input_manifest.csv"
    family = {row["cell_line"]: row for row in table(family_path)}
    aliases = table(alias_path)
    require(len(aliases) == 23 and len({x["cell_line"] for x in aliases}) == 23, "Expected 23 unique baseline models")
    require({x["cell_line"] for x in aliases} == {x["cell_line"] for x in table(manifest_path)}, "WES alias/manifest mismatch")
    sample_sheet_path = source / "scratch_archive_extract/samplesheet.csv"
    sample_sheet = {row["sample"]: row for row in table(sample_sheet_path)}
    require({x["cnv_sample_id"] for x in aliases} == set(sample_sheet), "Sarek sample set differs from baseline")
    for path, role in ((family_path, "canonical patient mapping"), (alias_path, "audited WES model/sample/passages"),
                       (manifest_path, "audited 23-model variant input manifest"), (sample_sheet_path, "Sarek input samplesheet"),
                       (source / "recovery/file_inventory.tsv", "cluster recovery checksum inventory")):
        register(path, role)

    bed = source / "scratch_archive_extract/intervals_sorted.bed"
    bed_stats = bed_accounting(bed)
    bed_id = register(bed, "mosdepth --by BED; overlaps retained")
    require(bed_stats["interval_count"] == 242421 and bed_stats["summed_bp"] == 63709951, "Unexpected recovered target BED")
    target_denom = f'{bed_stats["summed_bp"]} summed BED bases; overlapping intervals counted repeatedly'
    cnv_bed = source / "scratch_archive_extract/intervals_sorted.target.bed"
    cnv_bed_stats = bed_accounting(cnv_bed)
    register(cnv_bed, "CNVkit target bins, distinct from mosdepth BED")

    # Executed options, version and pinned implementation determine depth semantics.
    mos_command = next(run.glob("task_evidence/oct23_MOSDEPTH*/.command.sh"))
    mos_version = mos_command.with_name("versions.yml")
    command = mos_command.read_text()
    require("--by intervals_sorted.bed" in command and "mosdepth: 0.3.8" in mos_version.read_text(), "Unexpected mosdepth invocation/version")
    require(not re.search(r"(?:^|\s)(?:-F|--flag|-Q|--mapq|-x|--fast-mode|--use-median)(?:\s|=)", command), "Depth options changed; reassess semantics")
    register(mos_command, "executed mosdepth command; representative md task")
    register(mos_version, "executed mosdepth version")
    bqsr_command = next(run.glob("task_evidence/oct23_GATK4_APPLYBQSR*/.command.sh"))
    require("--intervals" in bqsr_command.read_text(), "Expected interval-restricted ApplyBQSR evidence")
    register(bqsr_command, "executed ApplyBQSR command showing interval restriction")

    fast_path, raw_path = mqc / "multiqc_fastp.txt", mqc / "multiqc_data.json"
    fast_rows = table(fast_path, "\t")
    fast = {row["Sample"].rsplit("-", 1)[0]: row for row in fast_rows}
    require(len(fast_rows) == 23 and set(fast) == set(sample_sheet), "Expected one fastp lane record per baseline model")
    fast_id = register(fast_path, "MultiQC-preserved raw fastp JSON fields; original JSONs omitted from bundle")
    raw_id = register(raw_path, "MultiQC JSON cross-check of raw fastp records and report date")
    raw = json.loads(raw_path.read_text())
    raw_fast = raw["report_saved_raw_data"]["multiqc_fastp"]
    mqc_source_path = mqc / "multiqc_sources.txt"
    register(mqc_source_path, "MultiQC original report source/sample associations")
    mqc_sources = table(mqc_source_path, "\t")
    fast_origin = {x["Sample Name"]: x["Source"] for x in mqc_sources if x["Module"].startswith("FastP")}
    mqc_sam = {r["Sample"]: r for r in table(mqc / "multiqc_samtools_stats.txt", "\t")}
    mqc_pic = {r["Sample"]: r for r in table(mqc / "multiqc_picard_dups.txt", "\t")}
    mqc_cov = {r["Sample"]: r for r in table(mqc / "mosdepth_cumcov_dist.txt", "\t")}
    mqc_chrom = {r["Sample"]: r for r in table(mqc / "mosdepth_perchrom.txt", "\t")}
    for filename in ("multiqc_samtools_stats.txt", "multiqc_picard_dups.txt", "mosdepth_cumcov_dist.txt", "mosdepth_perchrom.txt"):
        register(mqc / filename, "contemporaneous MultiQC cross-check")

    models, profiles = [], []
    comparisons = defaultdict(int)
    for alias in sorted(aliases, key=lambda x: x["cell_line"]):
        cell, sample = alias["cell_line"], alias["cnv_sample_id"]
        fam = family[cell]
        require(sample == cell.replace("-", "_") + "_" + alias["recorded_wes_passage"], f"Ambiguous sample alias: {cell}")
        require(alias["patient_id"] == fam["patient_id"] and alias["histotype"] == fam["subtype"], f"Patient/histotype mismatch: {cell}")
        row = dict(cell_line=cell, patient_id=fam["patient_id"], histotype=fam["subtype"], wes_sample_id=sample,
                   wes_passage=alias["recorded_wes_passage"], input_lane=sample_sheet[sample]["lane"],
                   target_bed_source_id=bed_id, target_summed_bp=bed_stats["summed_bp"], target_union_bp=bed_stats["union_bp"])
        fr = fast[sample]
        require(fr["Sample"] == sample + "-" + row["input_lane"], f"fastp lane mismatch: {sample}")
        row.update(fastp_source_id=fast_id, fastp_record_key=fr["Sample"], fastp_original_report=fast_origin[fr["Sample"]])
        parsed = {key: ast.literal_eval(fr[key]) for key in ("summary", "filtering_result", "duplication")}
        for key, val in parsed.items():
            require(val == raw_fast[fr["Sample"]][key], f"MultiQC TSV/JSON mismatch: {sample}/{key}")
            comparisons["fastp_tsv_json_fields"] += 1
        require(fr["command"] == raw_fast[fr["Sample"]]["command"], f"fastp command mismatch: {sample}")
        require("--disable_adapter_trimming" in fr["command"], f"Unexpected fastp adapter options: {sample}")
        before, after = parsed["summary"]["before_filtering"], parsed["summary"]["after_filtering"]
        filtering = parsed["filtering_result"]
        require(before["total_reads"] % 2 == after["total_reads"] % 2 == 0, f"Unpaired fastp count: {sample}")
        require(sum(filtering.values()) == before["total_reads"] and filtering["passed_filter_reads"] == after["total_reads"], f"fastp filtering accounting: {sample}")
        for label, entry in (("before", before), ("after", after)):
            add(row, f"fastp_reads_{label}", entry["total_reads"], "read ends (R1+R2)", "both reads from each paired fragment", "fastp", fast_id)
            add(row, f"fastp_pairs_{label}", entry["total_reads"] // 2, "paired fragments", "R1+R2 count divided by two; paired-end balance verified", "fastp", fast_id)
        add(row, "fastp_read_retention_fraction", after["total_reads"] / before["total_reads"], "fraction", "fastp reads before filtering (R1+R2)", "fastp", fast_id)
        add(row, "fastp_bases_after", after["total_bases"], "bases", "both mates after fastp filtering", "fastp", fast_id)
        add(row, "fastp_q30_bases_after", after["q30_bases"], "bases", "bases after filtering with Phred quality >=30", "fastp", fast_id)
        q30 = after["q30_bases"] / after["total_bases"]
        require(near(q30, after["q30_rate"], 0.00000051), f"fastp Q30 fraction: {sample}")
        add(row, "fastp_q30_fraction_after", q30, "fraction", "all bases after fastp filtering", "fastp", fast_id, "Recomputed from integer base counts; source q30_rate is rounded to six decimals.")
        add(row, "fastp_gc_fraction_after", after["gc_content"], "fraction", "all bases after fastp filtering", "fastp", fast_id, "Source precision six decimals.")
        for mate in (1, 2):
            add(row, f"fastp_read{mate}_mean_length_after", after[f"read{mate}_mean_length"], "bp", f"read{mate} after fastp filtering", "fastp", fast_id)

        stage_stats = {}
        for stage in ("md", "recal"):
            stats_path = run / f"reports/samtools/{sample}/{sample}.{stage}.cram.stats"
            stats_id = register(stats_path, f"samtools {stage} primary report", cell)
            st = samtools_stats(stats_path)
            stage_stats[stage] = st
            row[f"{stage}_samtools_source_id"] = stats_id
            denom = f"{stage} primary read ends, excluding secondary and supplementary alignments"
            require(st["reads mapped"] + st["reads unmapped"] == st["raw total sequences"], f"samtools mapped accounting: {sample}/{stage}")
            require(st["1st fragments"] + st["last fragments"] == st["raw total sequences"], f"samtools read-end accounting: {sample}/{stage}")
            for name, column in (("raw total sequences", "primary_reads"), ("reads mapped", "mapped_reads"),
                                 ("reads unmapped", "unmapped_reads"), ("reads properly paired", "properly_paired_reads"),
                                 ("reads MQ0", "mq0_reads"), ("reads duplicated", "duplicate_reads"),
                                 ("supplementary alignments", "supplementary_alignments")):
                count_denom = f"all {stage} alignment records; supplementary records excluded from primary-read fractions" if "alignments" in column else denom
                add(row, f"{stage}_{column}", st[name], "alignments" if "alignments" in column else "read ends", count_denom, stage, stats_id)
            add(row, f"{stage}_mapped_fraction", st["reads mapped"] / st["raw total sequences"], "fraction", denom, stage, stats_id)
            add(row, f"{stage}_properly_paired_fraction", st["reads properly paired"] / st["raw total sequences"], "fraction", denom, stage, stats_id)
            add(row, f"{stage}_mq0_fraction_of_mapped", st["reads MQ0"] / st["reads mapped"], "fraction", f"{stage} mapped primary read ends", stage, stats_id)
            if stage == "md":
                require(st["raw total sequences"] == after["total_reads"] and st["1st fragments"] == st["last fragments"] == after["total_reads"] // 2, f"fastp/samtools count mismatch: {sample}")
                for name, val in st.items():
                    key = name.replace(" ", "_")
                    if key in mqc_sam[sample + ".md"]:
                        require(near(val, mqc_sam[sample + ".md"][key], 1e-6), f"samtools/MultiQC mismatch: {sample}/{name}")
                        comparisons["samtools_summary_fields"] += 1

            depth_path = run / f"reports/mosdepth/{sample}/{sample}.{stage}.mosdepth.summary.txt"
            dist_path = run / f"reports/mosdepth/{sample}/{sample}.{stage}.mosdepth.region.dist.txt"
            depth_id = register(depth_path, f"mosdepth {stage} integer depth-sum and target length", cell)
            dist_id = register(dist_path, f"mosdepth {stage} per-base target cumulative fractions", cell)
            row[f"{stage}_mosdepth_summary_source_id"], row[f"{stage}_mosdepth_distribution_source_id"] = depth_id, dist_id
            summary = {x["chrom"]: x for x in table(depth_path, "\t")}
            target = summary["total_region"]
            if stage == "md":
                require(int(target["length"]) == bed_stats["summed_bp"], f"BED/summary denominator mismatch: {sample}/{stage}")
            included_length = 0
            for chrom, length in bed_stats["per_contig_summed_bp"].items():
                observed_length = int(summary.get(chrom + "_region", {}).get("length", 0))
                require(observed_length == length or (stage == "recal" and observed_length == 0), f"Contig target denominator mismatch: {sample}/{stage}/{chrom}")
                included_length += observed_length
            depth_sum, length = int(target["bases"]), int(target["length"])
            require(included_length == length, f"Per-contig lengths do not sum to total: {sample}/{stage}")
            row[f"mosdepth_{stage}_target_denominator_bp"] = length
            row[f"mosdepth_{stage}_target_omitted_no_read_contig_bp"] = bed_stats["summed_bp"] - length
            stage_target_denom = target_denom if stage == "md" else "sample-specific mosdepth_recal_target_denominator_bp; overlapping intervals repeated; no-read contigs omitted by mosdepth 0.3.8"
            require(near(depth_sum / length, target["mean"], 0.00501), f"mosdepth mean not depth-sum/length: {sample}/{stage}")
            add(row, f"mosdepth_{stage}_target_depth_sum", depth_sum, "read-base depth sum", stage_target_denom, stage, depth_id)
            add(row, f"mosdepth_{stage}_mean_target_depth_x", depth_sum / length, "x", stage_target_denom, stage, depth_id,
                "Recomputed from integer sum/length; duplicates excluded (FLAG 1796), MAPQ>=0, overlapping mates corrected; source mean printed to two decimals.")
            dist, rounded_means = distribution(dist_path)
            for threshold in THRESHOLDS:
                add(row, f"mosdepth_{stage}_fraction_target_ge_{threshold}x", dist[threshold], "fraction", stage_target_denom, stage, dist_id,
                    "Cumulative per-base BED coverage, rounded to two decimals; exact base count unavailable, do not multiply into a claimed exact count.")
                if stage == "md" and mqc_cov[sample + ".md"].get(str(threshold)):
                    require(near(100 * dist[threshold], mqc_cov[sample + ".md"][str(threshold)], 1e-7), f"mosdepth/MultiQC coverage mismatch: {sample}/{threshold}")
                    comparisons["mosdepth_coverage_points"] += 1
            if stage == "md":
                for chrom, val in mqc_chrom[sample + ".md"].items():
                    if chrom != "Sample" and val:
                        require(near(rounded_means[chrom], val, 1e-7), f"mosdepth/MultiQC CDF-derived contig mean mismatch: {sample}/{chrom}")
                        comparisons["mosdepth_cdf_derived_contig_means"] += 1
            for threshold in range(151):
                profiles.append(dict(cell_line=cell, patient_id=fam["patient_id"], wes_sample_id=sample,
                                     stage=stage, depth_x=threshold,
                                     fraction_ge_depth=dist[threshold], fraction_decimal_places=2,
                                     target_summed_bp=length, source_id=dist_id))

        add(row, "recal_primary_read_retention_fraction", stage_stats["recal"]["raw total sequences"] / stage_stats["md"]["raw total sequences"],
            "fraction", "genome-wide md primary read ends", "recal_vs_md", row["recal_samtools_source_id"],
            "Interval restriction after per-shard ApplyBQSR; not a filtering-loss or genome-wide on-target-rate estimate.")
        add(row, "recal_primary_reads_not_retained", stage_stats["md"]["raw total sequences"] - stage_stats["recal"]["raw total sequences"],
            "read ends", "genome-wide md primary read ends minus interval-restricted recal primary read ends", "recal_vs_md", row["recal_samtools_source_id"])
        for derived in ("recal_primary_read_retention_fraction", "recal_primary_reads_not_retained"):
            metric_sources[derived].add(row["md_samtools_source_id"])

        pic_path = run / f"reports/markduplicates/{sample}/{sample}.md.cram.metrics"
        pic_id = register(pic_path, "Picard/GATK MarkDuplicates primary library metrics", cell)
        pic = picard_metrics(pic_path)
        require(pic["LIBRARY"] == sample, f"Picard library alias mismatch: {sample}")
        row["picard_source_id"] = pic_id
        for name, val in pic.items():
            if name == "LIBRARY":
                require(val == mqc_pic[sample + ".md"][name], f"Picard/MultiQC library mismatch: {sample}")
            else:
                require(near(val, mqc_pic[sample + ".md"][name], 1e-8), f"Picard/MultiQC mismatch: {sample}/{name}")
                comparisons["picard_fields"] += 1
        examined = int(pic["UNPAIRED_READS_EXAMINED"]) + 2 * int(pic["READ_PAIRS_EXAMINED"])
        duplicates = int(pic["UNPAIRED_READ_DUPLICATES"]) + 2 * int(pic["READ_PAIR_DUPLICATES"])
        require(examined == stage_stats["md"]["reads mapped"], f"Picard/samtools mapped denominator mismatch: {sample}")
        require(duplicates == stage_stats["md"]["reads duplicated"], f"Picard/samtools duplicate count mismatch: {sample}")
        require(near(duplicates / examined, pic["PERCENT_DUPLICATION"], 0.00000051), f"Picard duplicate fraction formula: {sample}")
        add(row, "picard_examined_mapped_read_ends", examined, "read ends", "UNPAIRED_READS_EXAMINED + 2*READ_PAIRS_EXAMINED", "md", pic_id)
        add(row, "picard_duplicate_read_ends", duplicates, "read ends", "UNPAIRED_READ_DUPLICATES + 2*READ_PAIR_DUPLICATES", "md", pic_id)
        add(row, "picard_duplication_fraction", float(pic["PERCENT_DUPLICATION"]), "fraction", "Picard examined mapped read ends", "md", pic_id,
            "Despite its name PERCENT_DUPLICATION is a 0–1 fraction, reported to six decimals. Duplicate reads were marked, not removed from the CRAM.")
        con_path = run / f"variant_calling/mutect2/{sample}/{sample}.mutect2.contamination.table"
        con_id = register(con_path, "GATK CalculateContamination primary table", cell)
        contamination = table(con_path, "\t")
        require(len(contamination) == 1 and contamination[0]["sample"] == f"{sample}_{sample}", f"Contamination sample mapping: {sample}")
        row["contamination_source_id"] = con_id
        for field, name in (("contamination", "gatk_contamination_fraction"), ("error", "gatk_contamination_reported_error")):
            add(row, name, float(contamination[0][field]), "fraction", "GATK model using GetPileupSummaries at population-resource sites in calling intervals", "tumour_only_calling", con_id,
                "Caller estimate/error field; not evidence of stock authentication or mycoplasma status.")
        models.append(row)

    require(len(models) == 23 and len({r["patient_id"] for r in models}) == 16, "Expected 23 models from 16 patients")
    checks.extend(["23 unique baseline sample/library/contamination identifiers map to 16 canonical patients",
                   "23 paired fastp read counts agree exactly with genome-wide samtools primary counts and balanced R1/R2",
                   "23 Picard fractions and count formulae agree with primary samtools counts and contemporaneous MultiQC",
                   "23 md target denominators match the full summed BED; all 46 stage denominators match included per-contig lengths",
                   "46 cumulative distributions are bounded, monotonic and complete at 0–150x",
                   "Available original-report SHA-256 values match the recovered inventory"])
    summaries = []
    for column, definition in metric_defs.items():
        vals = [value for _, value in values[column]]
        summaries.append(dict(cohort="baseline_models", metric=column, n=len(vals), n_primary_verified=len(vals),
                              minimum=min(vals), median=statistics.median(vals), maximum=max(vals),
                              minimum_models=";".join(cell for cell, value in values[column] if value == min(vals)),
                              maximum_models=";".join(cell for cell, value in values[column] if value == max(vals)),
                              status="validated_recovered_reports", source_ids=";".join(sorted(metric_sources[column])), **definition))

    # The small handoff omits normal coverage CNNs. Preserve that limitation instead
    # of upgrading the recovery generator's arithmetic into primary validation.
    recovered_qc_path = source / "recovery/qc_metrics.tsv"
    recovered_qc_id = register(recovered_qc_path, "cluster-produced recovery table; secondary evidence only for omitted normal CNNs")
    recovered_qc = table(recovered_qc_path, "\t")
    normal_rows = []
    for normal in NORMALS:
        reported = [r for r in recovered_qc if r["model_or_reference"] == normal and r["metric"] == "cnvkit_target_mean_depth_length_weighted"]
        require(len(reported) == 1, f"Expected one reported normal target mean: {normal}")
        report = reported[0]
        cnn = source / report["source_path"].removeprefix("wes/")
        cnn_id = register(cnn, "CNV-reference targetcoverage CNN required for independent validation", normal, absent_ok=True)
        value, status = float(report["value"]), "reported_only_primary_missing"
        if cnn.is_file():
            rows = table(cnn, "\t")
            length = sum(int(r["end"]) - int(r["start"]) for r in rows)
            require(len(rows) == cnv_bed_stats["interval_count"] and length == cnv_bed_stats["summed_bp"], f"CNV-reference bin denominator: {normal}")
            value = sum((int(r["end"]) - int(r["start"])) * float(r["depth"]) for r in rows) / length
            require(near(value, report["value"], 0.00501), f"CNV-reference reported mean mismatch: {normal}")
            status = "validated_primary_coverage_cnn"
        normal_rows.append(dict(reference=normal, target_mean_depth_x=value, status=status,
                                coverage_source_id=cnn_id, reported_summary_source_id=recovered_qc_id,
                                denominator=report["denominator"]))
    normals_verified = sum(r["status"] == "validated_primary_coverage_cnn" for r in normal_rows)
    normal_values = [r["target_mean_depth_x"] for r in normal_rows]
    summaries.append(dict(cohort="cnv_reference_exomes", metric="cnvkit_normal_mean_target_depth_x", n=5,
                          n_primary_verified=normals_verified, minimum=min(normal_values), median=statistics.median(normal_values), maximum=max(normal_values),
                          minimum_models=";".join(r["reference"] for r in normal_rows if r["target_mean_depth_x"] == min(normal_values)),
                          maximum_models=";".join(r["reference"] for r in normal_rows if r["target_mean_depth_x"] == max(normal_values)),
                          status="validated_primary_coverage_cnn" if normals_verified == 5 else "reported_only_or_partially_verified",
                          source_ids=";".join([recovered_qc_id] + [r["coverage_source_id"] for r in normal_rows]), unit="x",
                          denominator=f'{cnv_bed_stats["summed_bp"]} bp in {cnv_bed_stats["interval_count"]} CNVkit target bins',
                          alignment_stage="independently_aligned_CNV_reference_exomes",
                          notes="Secondary recovery-table values are retained only with explicit verification status; do not mix with mosdepth model depth."))

    report = dict(schema_version=1, source_bundle=str(source), baseline_models=23, baseline_patients=16,
                  multiqc_creation_date=raw["config_creation_date"], multiqc_version=raw["config_version"],
                  target_bed=bed_stats, cnv_target_bed=cnv_bed_stats, metric_definitions=metric_defs,
                  checks=checks, crosscheck_counts=dict(comparisons), sources=len(sources),
                  source_files_hash_checked=sum(r["inventory_hash_match"] is True for r in sources),
                  reference_normals=normal_rows, reference_normals_primary_verified=normals_verified,
                  method_references=[dict(url=MOSDEPTH_SOURCE, version="0.3.8", lines="219–254, 355–413, 612–627, 737–754",
                                          facts="FLAG 1796 excludes unmapped, secondary, QC-failed and duplicate records; MAPQ default0; overlapping mates corrected without fast-mode; BED distribution uses per-base depths; precision defaults2")],
                  limitations=["Coverage fractions have two decimal places; exact covered-base counts cannot be recovered from these summaries.",
                               "Target denominator sums overlapping BED intervals; a distinct union-footprint coverage calculation was not performed.",
                               "Recal summaries omit 507–3843 target bases on no-read contigs; each actual denominator is retained. All md summaries cover the full summed BED length.",
                               "Archived MultiQC per-contig values equal sums of printed cumulative fractions excluding the zero-depth level; rounded/sparse high-depth levels make these differ from the primary integer-sum/length means used here.",
                               "Original fastp JSON files are omitted from this bundle; their raw fields persist consistently in both MultiQC TSV and JSON.",
                               "Primary metrics files may reflect later identical-workflow reruns; md numerical QC was reconciled against the Oct 23 MultiQC report, not assumed byte-identical to calling alignments.",
                               "No stock-authentication inference follows from the contamination estimates.",
                               "Reference-normal target means remain secondary reports if their coverage CNNs are absent; no normal alignment QC is present."],
                  corrections=["Recovery generator used the wrong Picard path and omitted all 23 model duplication records; nested primary reports now parsed.",
                               "Recovery prose saying mosdepth counted duplicate reads is incorrect for the executed default FLAG 1796.",
                               "Recal read loss is interval restriction, not evidence of failed filtering; restricted-BAM target-read percentages are not genome-wide on-target rates."],
                  model_findings=[dict(cell_line="TOV2929D", finding="Lowest mean target depth and >=10/20/30x coverage in cohort; inspect locus coverage before interpreting negative calls; no automatic exclusion justified by these summaries."),
                                  dict(cell_line="OV3291", finding="Largest contamination estimate remains below 0.005; model-based estimate alone does not establish contamination failure."),
                                  dict(cell_line="TOV2881EP", finding="Lowest post-fastp Q30 fraction; report distribution without inferring a failed library."),
                                  dict(cell_line="TOV81D", finding="fastp-before count 130101168 read ends = 65050584 pairs; not an exact match to provider 65085803 reads/pairs proposed by acquisition review.")])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "wes_qc_model_summary.csv", models)
    write_csv(args.output_dir / "wes_qc_metric_summary.csv", summaries)
    write_csv(args.output_dir / "wes_qc_sources.csv", sources)
    write_csv(args.output_dir / "wes_qc_coverage_profile.csv", profiles)
    (args.report_dir / "qc_validation.json").write_text(json.dumps(report, indent=2) + "\n")

    primary = {r["metric"]: r for r in summaries if r["cohort"] == "baseline_models"}
    selected = ["fastp_reads_after", "fastp_pairs_after", "fastp_q30_fraction_after", "fastp_read_retention_fraction", "md_mapped_fraction", "md_properly_paired_fraction", "picard_duplication_fraction", "mosdepth_md_mean_target_depth_x", "mosdepth_md_fraction_target_ge_10x", "mosdepth_md_fraction_target_ge_20x", "mosdepth_md_fraction_target_ge_30x", "gatk_contamination_fraction", "recal_primary_read_retention_fraction", "md_mq0_fraction_of_mapped"]
    lines = ["# Recovered WES QC validation", "", "Validated all **23 baseline models from 16 patients** against the recovered primary reports and contemporaneous MultiQC records. No alignments were read, analyses rerun, or existing scientific result tables modified.", "", "## Quantitative summary", "", "All rows below have n=23. Fractions are on 0–1 scales. Full machine precision and per-model extrema are in `output/wes_qc_metric_summary.csv`; integer numerators remain in `output/wes_qc_model_summary.csv`.", "", "| Metric | Minimum | Median | Maximum | Unit |", "| --- | ---: | ---: | ---: | --- |"]
    for name in selected:
        r = primary[name]
        fmt = lambda value: f"{value:,}" if isinstance(value, int) else f"{value:.9g}"
        lines.append(f'| `{name}` | {fmt(r["minimum"])} | {fmt(r["median"])} | {fmt(r["maximum"])} | {r["unit"]} |')
    lines += ["", "## Denominators and source semantics", "",
              "- fastp counts are R1+R2 read ends. Division by two is supported by equal first/last read-end counts in all 23 genome-wide samtools reports. Raw fastp summary/filtering/duplication fields agree between the MultiQC TSV and JSON. Original standalone fastp JSONs were not included in the handoff.",
              "- Mapping and proper-pair fractions use genome-wide duplicate-marked (`md`) primary read ends, excluding secondary/supplementary alignments from that denominator. Picard duplication is a fraction despite the field name `PERCENT_DUPLICATION`; its numerator and denominator were independently rebuilt from paired/unpaired counts and reconciled to samtools.",
              f'- mosdepth 0.3.8 ran `--by intervals_sorted.bed` with no flag, mapping-quality or fast-mode override. Its [version-pinned implementation]({MOSDEPTH_SOURCE}#L219-L254) excludes FLAG 1796 (unmapped, secondary, QC-failed and duplicate alignments), uses MAPQ>=0 and corrects overlapping mates. Supplementary alignments are not excluded by 1796. This corrects the handoff prose claiming duplicate-inclusive depth.',
              f'- The [BED-mode distribution code]({MOSDEPTH_SOURCE}#L612-L627) accumulates per-base depth, rather than the means used for numeric windows. The md denominator is **{bed_stats["summed_bp"]:,} summed BED bases across {bed_stats["interval_count"]:,} intervals**, not the **{bed_stats["union_bp"]:,} bp union**; {bed_stats["overlap_counted_extra_bp"]:,} bp are counted again because of overlaps. All 23 md summaries include the full summed length. Recal summaries omit 507–3,843 bp on no-read contigs; their sample-specific lengths are retained. Every stage denominator was reconciled to its included per-contig lengths.',
              "- Means were reconstructed from integer depth sums divided by integer target lengths. The distribution files retain only two decimals: coverage fractions are approximate to 0.01, and exact threshold-covered base counts cannot be reconstructed. No spurious exact counts were calculated from those rounded fractions.",
              "- Archived MultiQC per-contig values were reproduced as sums of the printed cumulative fractions excluding the zero-depth level. Because the depth distribution is rounded and omits unchanged high-depth levels, these plotted values differ from the authoritative integer-sum/length means. The latter are used in the new model summary and manuscript.",
              "- Recalibrated CRAM QC describes alignments retained after interval-sharded ApplyBQSR. Its lower primary-read count is not a quality-filter failure. A target-read percentage calculated from those restricted alignments is not a genome-wide on-target fraction. Later report/CRAM rewrites are distinguished from the Oct 23 calling run; agreement with contemporaneous MultiQC establishes numerical QC consistency, not alignment byte identity.",
              "- GATK contamination tables map uniquely through their repeated patient/sample IDs for all 23 models. The estimate and reported error field are retained without calling the error a confidence interval or inferring stock authenticity.", "", "## Source checks and corrections", "",
              f'- Checked {report["source_files_hash_checked"]} present source-file hashes against the recovery inventory. All matched. `output/wes_qc_sources.csv` contains local SHA-256 values, cluster paths, inventory hashes and explicit missing-file status.',
              f'- Cross-checks: {dict(comparisons)}. Every compared value passed. The metric definitions, checks and source IDs are also recorded in `qc_validation.json`.',
              "- The recovery generator searched the wrong Picard path, omitting all 23 model duplication records. This builder reads the actual nested sample directories and supplies the missing values.",
              "", "## Five CNV reference exomes", "",
              f'Primary coverage CNNs independently validated: **{normals_verified}/5**. The recovery table reports length-weighted means of {min(normal_values):.2f}–{max(normal_values):.2f}x (median {statistics.median(normal_values):.2f}x), using {cnv_bed_stats["summed_bp"]:,} bp in {cnv_bed_stats["interval_count"]:,} CNVkit bins. These are a distinct method/denominator from model mosdepth. Missing primary CNNs prevent independently endorsing their reported arithmetic; missing paths and inventory hashes are preserved. No normal alignment, duplication or contamination QC was supplied.',
              "", "## Actionable interpretation", "",
              "TOV2929D is consistently the lowest-depth/lowest-coverage model (69.7151x mean; rounded fractions 0.91/0.79/0.66 at 10/20/30x). Locus-level depth should be checked before interpreting an absent variant in this model. These aggregate summaries alone do not justify exclusion or demonstrate a biological cause. All contamination estimates are below 0.005; the largest is OV3291. A reference-bundle mapping fraction includes MAPQ0 alignments (12.26–15.95% of mapped primary read ends), so high alignment rate is not synonymous with uniquely usable target coverage.",
              "", "## Figure suggestion and reproduction", "",
              "Use two compact panels: all 23 per-model mean target depths, and rounded cumulative target-coverage curves (or the 10/20/30x distributions), with a consistent model order. Mean-depth and coverage filters/denominators belong in the external legend. Source-ready curves for both stages are in `output/wes_qc_coverage_profile.csv`; the primary manuscript should use `stage=md`.",
              "", "```sh", "python3 scripts/24_wes_recovered_qc.py", "```", "", "The script uses only the Python standard library and validates source hashes, aliases, denominators, filtering accounting and cross-report agreement before writing its new outputs. It does not modify metadata, existing analysis outputs, figures, manuscript or release records."]
    (args.report_dir / "qc_validation.md").write_text("\n".join(lines) + "\n")
    print(f"PASS: {len(models)} models, {len(profiles)} coverage points, {len(summaries)} metric summaries, {len(sources)} source records; normals independently verified {normals_verified}/5")


if __name__ == "__main__":
    main()
