#!/usr/bin/env python3
"""Reproduce the bounded OV1369-R2 centering sensitivity without changing calls.

Standard-library Python performs span arithmetic; base R and matrixStats verify
the exact median convention used by script08. OVCAN_DATA can restore source CNR
and CNS files from an external archived data directory.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MODEL = "OV1369-R2"
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
csv.field_size_limit(10_000_000)
R_CHECK = r'''
args <- commandArgs(trailingOnly=TRUE)
old <- read.delim(args[1], stringsAsFactors=FALSE)
new <- read.delim(args[2], stringsAsFactors=FALSE)
bins <- read.delim(args[3], stringsAsFactors=FALSE)
auto <- paste0("chr", 1:22)
old <- old[old$chromosome %in% auto, ]
new <- new[new$chromosome %in% auto, ]
bins <- bins[bins$chromosome %in% auto & bins$gene != "Antitarget" & bins$depth > 0, ]
values <- c(n=nrow(bins),
  bin_unweighted=median(bins$log2),
  bin_CNR_weighted=matrixStats::weightedMedian(bins$log2, bins$weight),
  archived_segment_center=matrixStats::weightedMedian(old$log2, old$probes),
  target_only_segment_center=matrixStats::weightedMedian(new$log2, new$probes))
for (key in names(values)) cat(key, sprintf("%.17g", values[[key]]), sep="\t", fill=TRUE)
'''


def read(path, delimiter=","):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_path(path):
    p = Path(path)
    prefix = Path("judy_archive/data")
    if os.environ.get("OVCAN_DATA") and p.is_relative_to(prefix):
        return Path(os.environ["OVCAN_DATA"]) / p.relative_to(prefix)
    return ROOT / p


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def totals(segments, center):
    denominator = sum(int(r["end"]) - int(r["start"]) for r in segments)
    gained = sum(int(r["end"]) - int(r["start"]) for r in segments if float(r["log2"]) - center > 0.20)
    lost = sum(int(r["end"]) - int(r["start"]) for r in segments if float(r["log2"]) - center < -0.20)
    return denominator, gained, lost


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", type=Path, default=shutil.which("Rscript"))
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reports/wes_cnv_coverage_2026-09-06")
    args = parser.parse_args()
    if not args.rscript:
        parser.error("Rscript with matrixStats is required")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    baseline_path = ROOT / "reports/wes_cnv_coverage_2026-09-06/provenance_archived_centering_baseline.csv"
    baseline = next(r for r in read(baseline_path) if r["cell_line"] == MODEL)
    original = next(r for r in read(ROOT / "output/wes_recovered_provenance_cnv_support.csv") if r["cell_line"] == MODEL)
    corrected = next(r for r in read(ROOT / "output/wes_cnv_target_only/manifest.csv") if r["cell_line"] == MODEL)
    old_path, cnr_path = source_path(original["cns_source"]), source_path(original["cnr_source"])
    new_path = ROOT / corrected["cns_path"]
    assert sha(old_path) == original["cns_sha256"] == baseline["archived_cns_sha256"]
    assert sha(cnr_path) == original["cnr_sha256"] == corrected["source_cnr_sha256"]
    assert sha(new_path) == corrected["cns_sha256"]
    output = subprocess.check_output([str(args.rscript), "--vanilla", "-e", R_CHECK,
                                      str(old_path), str(new_path), str(cnr_path)], text=True)
    values = {line.split()[0]: float(line.split()[1]) for line in output.strip().splitlines()}
    assert all(math.isfinite(x) for x in values.values())
    old_center, new_center = values["archived_segment_center"], values["target_only_segment_center"]
    assert abs(old_center - float(baseline["archived_autosome_center"])) < 1e-12
    assert int(values["n"]) == 195544
    old = [r for r in read(old_path, "\t") if r["chromosome"] in AUTOSOMES]
    new = [r for r in read(new_path, "\t") if r["chromosome"] in AUTOSOMES]
    checks = []
    for name, segments in (("archived", old), ("target_only", new)):
        for center_name, center in (("archived", old_center), ("target_only", new_center)):
            denominator, gained, lost = totals(segments, center)
            checks.append(dict(segment_source=name, center_source=center_name, center=center,
                n_autosomal_segments=len(segments), autosomal_span_bp=denominator,
                gained_bp=gained, lost_bp=lost, fga=(gained + lost) / denominator))
    crossings = []
    for r in new:
        if float(r["log2"]) - old_center >= -0.20 and float(r["log2"]) - new_center < -0.20:
            crossings.append({k: r[k] for k in ("chromosome", "start", "end", "log2", "probes")} |
                dict(span_bp=int(r["end"]) - int(r["start"]),
                     log2_at_archived_center=float(r["log2"]) - old_center,
                     log2_at_new_center=float(r["log2"]) - new_center))
    scenarios = [("archived_segment_probe_weighted_median", old_center),
                 ("target_only_segment_probe_weighted_median", new_center),
                 ("positive_target_bin_unweighted_median", values["bin_unweighted"]),
                 ("positive_target_bin_CNR_weighted_median", values["bin_CNR_weighted"]),
                 ("no_additional_centering", 0)]
    alternatives = []
    for label, center in scenarios:
        denominator, gained, lost = totals(new, center)
        alternatives.append(dict(center_definition=label, center=center,
            fga_on_target_only_segments=(gained + lost) / denominator,
            frac_gain=gained / denominator, frac_loss=lost / denominator,
            autosomal_span_bp=denominator))
    write_csv(out / "ov1369_centering_2x2.csv", checks)
    write_csv(out / "ov1369_centering_threshold_crossings.csv", crossings)
    write_csv(out / "ov1369_alternative_centers.csv", alternatives)
    summary = dict(script="scripts/29c_cnv_centering_sensitivity.py", model=MODEL,
        archived_center=old_center, target_only_center=new_center, center_difference=new_center-old_center,
        independent_fga_target_only=checks[-1]["fga"], loss_threshold_crossings=len(crossings),
        crossing_bp=sum(r["span_bp"] for r in crossings),
        autosomal_span_difference_bp=checks[-1]["autosomal_span_bp"]-checks[0]["autosomal_span_bp"],
        sources=[dict(path=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p), sha256=sha(p))
                 for p in (old_path, new_path, cnr_path, baseline_path)],
        formula="FGA=sum(end-start for autosomal segments where abs(log2-center)>0.20)/sum(end-start); retain segment sums as in script08, not a base-resolution truth claim",
        bin_center_check=dict(positive_autosomal_target_bins=int(values["n"]),
            unweighted_median=values["bin_unweighted"], CNR_weighted_median=values["bin_CNR_weighted"],
            calculation="Base R median and matrixStats::weightedMedian(log2, weight), after gene != Antitarget, depth > 0, chromosome chr1..chr22"))
    (out / "ov1369_centering_diagnostic.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
