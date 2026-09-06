#!/usr/bin/env python3
"""Read-only provenance/antitarget checks for the September CNV coverage return.

No historical command is executed and no biological output is replaced. Raw
bin-weighted means using FIXED existing segment probe membership are diagnostic;
they do not rerun CBS, its smoothing or its breakpoint selection. Membership is
reconstructed from ordered positive-depth CNR rows and CNS probe counts and is
verified against every archived CNS log2 mean to its printed precision.
"""
import argparse
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(10_000_000)
ROOT = Path(__file__).resolve().parents[1]


def rows(path, delimiter="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rel(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write_csv(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def span(row):
    return row["chromosome"], int(row["start"]), int(row["end"])


def bounds(values):
    return {"min": min(values), "median": statistics.median(values), "max": max(values)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=ROOT / "data/cluster_wes_retrieval/2026-09-06/ovcan_human_wes_cnv_coverage_2026-09-06")
    parser.add_argument("--previous-bundle", type=Path, default=ROOT / "data/cluster_wes_retrieval/2026-09-05/ovcan_human_wes_handoff_2026-09-05")
    parser.add_argument("--report-dir", type=Path, default=ROOT / "reports/wes_cnv_coverage_2026-09-06")
    args = parser.parse_args()
    bundle, previous, out = args.bundle.resolve(), args.previous_bundle.resolve(), args.report_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    all_new = sorted(p for p in bundle.rglob("*") if p.is_file())
    assert len(all_new) == 109
    old_by_name = defaultdict(list)
    for p in (ROOT / "judy_archive/data/cnvkit wes - new").rglob("*"):
        if p.is_file():
            old_by_name[p.name].append(p)
    inventory = rows(previous / "recovery/file_inventory.tsv")
    identities = []
    for p in all_new:
        kind = None
        comparators = []
        if p.name in ("commands.txt", "commands1.txt", "command_diagram.txt", "run_cnvkit.sh"):
            kind = "previous_returned_command"
            other = previous / "cnvkit_0.9.10_manual" / p.name
            assert other.exists(), other
            comparators = [(rel(other), sha(other), "current file bytes")]
        elif p.suffix == ".cnn" or p.suffix == ".bed":
            if p.name.startswith("SRR"):
                kind = "normal_CNN_historical_inventory"
                comparators = [(r["original_path"], r["checksum"], r["checksum_evidence"]) for r in inventory
                               if Path(r["original_path"]).name == p.name]
                assert len(comparators) == 23, (p.name, len(comparators))
            else:
                kind = "archived_reference" if p.name == "reference.cnn" else "archived_target_BED" if p.suffix == ".bed" else "archived_model_CNN"
                comparators = [(rel(other), sha(other), "current file bytes") for other in old_by_name[p.name]]
                assert comparators, p
        if kind:
            checksum = sha(p)
            assert all(c[1] == checksum for c in comparators), p
            identities.append(dict(kind=kind, returned_path=rel(p), bytes=p.stat().st_size, sha256=checksum,
                                   comparator_count=len(comparators), all_comparators_match=True,
                                   comparator_paths=" | ".join(c[0] for c in comparators),
                                   comparator_checksum_evidence=" | ".join(sorted(set(c[2] for c in comparators)))))
    assert Counter(r["kind"] for r in identities) == {"previous_returned_command": 4, "archived_model_CNN": 46,
                                                      "normal_CNN_historical_inventory": 10, "archived_reference": 1, "archived_target_BED": 1}
    write_csv(out / "provenance_file_identity.csv", identities)

    bed = bundle / "normals/intervals_sorted.target.bed"
    targets = defaultdict(list)
    for line in bed.read_text().splitlines():
        if line and not line.startswith(("#", "track", "browser")):
            chrom, start, end = line.split("\t")[:3]
            targets[chrom].append((int(start), int(end)))
    support = rows(ROOT / "output/wes_recovered_provenance_cnv_support.csv", ",")
    assert len(support) == 23
    # Immutable 23-row snapshot of archived centers, pinned before correction.
    # Do not use canonical segments here: subsequent target-only correction is
    # expected to replace those outputs. The baseline links its source checksum.
    baselines = {r["cell_line"]: r for r in rows(out / "provenance_archived_centering_baseline.csv", ",")}
    assert set(baselines) == {r["cell_line"] for r in support}
    context, model_summaries, segment_diagnostics = [], [], []
    segment_mean_errors = []
    positive_coords = None
    for model in support:
        cnr_path, cns_path = ROOT / model["cnr_source"], ROOT / model["cns_source"]
        assert sha(cnr_path) == model["cnr_sha256"]
        assert sha(cns_path) == model["cns_sha256"]
        assert baselines[model["cell_line"]]["archived_cns_sha256"] == model["cns_sha256"]
        cnr = rows(cnr_path)
        cns = rows(cns_path)
        positive = [r for r in cnr if float(r["depth"]) > 0]
        antitargets = [r for r in positive if r["gene"] == "Antitarget"]
        target_bins = [r for r in positive if r["gene"] != "Antitarget"]
        coords = set(map(span, antitargets))
        assert len(antitargets) == 4
        if positive_coords is None:
            positive_coords = coords
        assert coords == positive_coords
        target_weights = [float(r["weight"]) for r in target_bins]
        model_summaries.append(dict(cell_line=model["cell_line"], sample_id=model["sample_id"],
            positive_target_bins=len(target_bins), positive_antitarget_bins=4,
            target_weight_min=min(target_weights), target_weight_median=statistics.median(target_weights),
            target_weight_max=max(target_weights), antitarget_weight_min=min(float(r["weight"]) for r in antitargets),
            antitarget_weight_max=max(float(r["weight"]) for r in antitargets)))
        # Probe counts partition ordered positive-depth CNR records. This avoids
        # assigning an overlapping antitarget to two segments by coordinates.
        # Mean reconstruction certifies this membership without rerunning CBS.
        relevant = [s for s in cns if any(span(a)[0] == span(s)[0] and span(a)[1] < span(s)[2] and span(a)[2] > span(s)[1] for a in antitargets)]
        metrics = {}
        cursor = 0
        for s in cns:
            chrom, start, end = span(s)
            n_probes = int(s["probes"])
            bins = positive[cursor:cursor + n_probes]
            cursor += n_probes
            assert len(bins) == n_probes and all(r["chromosome"] == chrom for r in bins)
            all_weight = sum(float(r["weight"]) for r in bins)
            all_mean = sum(float(r["weight"]) * float(r["log2"]) for r in bins) / all_weight
            mean_error = abs(all_mean - float(s["log2"]))
            segment_mean_errors.append(mean_error)
            assert mean_error < 0.0000051, (model["cell_line"], span(s), mean_error)
            if s not in relevant:
                continue
            tbins = [r for r in bins if r["gene"] != "Antitarget"]
            abins = [r for r in bins if r["gene"] == "Antitarget"]
            tw = sum(float(r["weight"]) for r in tbins)
            aw = sum(float(r["weight"]) for r in abins)
            target_mean = sum(float(r["weight"]) * float(r["log2"]) for r in tbins) / tw
            combined_mean = sum(float(r["weight"]) * float(r["log2"]) for r in bins) / (tw + aw)
            diag = dict(cell_line=model["cell_line"], sample_id=model["sample_id"], chromosome=chrom,
                segment_start=start, segment_end=end, segment_log2=float(s["log2"]), segment_probes=int(s["probes"]),
                segment_recorded_weight=float(s["weight"]), contributing_target_bins=len(tbins), contributing_antitarget_bins=len(abins),
                contributing_target_weight=tw, contributing_antitarget_weight=aw,
                antitarget_fraction_contributing_weight=aw / (tw + aw),
                target_only_raw_weighted_mean=target_mean, all_bins_raw_weighted_mean=combined_mean,
                raw_weighted_mean_change_removing_antitargets=target_mean - combined_mean,
                archived_segment_minus_target_only_raw_mean=float(s["log2"]) - target_mean)
            center = float(baselines[model["cell_line"]]["archived_autosome_center"])
            archived_centered = float(s["log2"]) - center
            corrected_at_same_center = target_mean - center
            old_gain, new_gain = archived_centered > 0.20, corrected_at_same_center > 0.20
            old_loss, new_loss = archived_centered < -0.20, corrected_at_same_center < -0.20
            diag.update(archived_autosome_centered_log2=archived_centered,
                        target_only_log2_at_unchanged_autosome_center=corrected_at_same_center,
                        autosomal_gain_loss_label_changes_at_fixed_center=bool(chrom in {f"chr{i}" for i in range(1, 23)} and (old_gain != new_gain or old_loss != new_loss)))
            metrics[(chrom, start, end)] = (diag, set(map(span, abins)))
            segment_diagnostics.append(diag)
        assert cursor == len(positive)
        for a in antitargets:
            chrom, start, end = span(a)
            intersections = [max(0, min(end, t1) - max(start, t0)) for t0, t1 in targets[chrom]]
            n_overlap, bp_overlap = sum(x > 0 for x in intersections), sum(intersections)
            assert bp_overlap > 0
            for s in relevant:
                sc, ss, se = span(s)
                if chrom == sc and start < se and end > ss:
                    context.append(dict(cell_line=model["cell_line"], sample_id=model["sample_id"], chromosome=chrom,
                        antitarget_start=start, antitarget_end=end, antitarget_length=end-start,
                        overlapping_target_bins=n_overlap, target_overlap_bp=bp_overlap,
                        depth=float(a["depth"]), cnr_log2=float(a["log2"]), cnr_weight=float(a["weight"]),
                        segment_start=ss, segment_end=se, segment_log2=float(s["log2"]), segment_probes=int(s["probes"]),
                        segment_recorded_weight=float(s["weight"]),
                        bin_weight_fraction_of_segment_recorded_weight=float(a["weight"]) / float(s["weight"]),
                        bin_is_contributor_to_segment=span(a) in metrics[(sc, ss, se)][1],
                        cnr_source=model["cnr_source"], cns_source=model["cns_source"]))
    assert len(context) == 94
    write_csv(out / "provenance_antitarget_context.csv", context)
    write_csv(out / "provenance_segment_influence.csv", segment_diagnostics)
    write_csv(out / "provenance_model_bin_weights.csv", model_summaries)
    # Do not mistake the bundled current Slurm wrapper for a successful run log.
    commands1 = next(p for p in all_new if p.name == "commands1.txt")
    wrapper = next(p for p in all_new if p.name == "run_cnvkit.sh")
    assert "--array=1-27" in wrapper.read_text()
    assert 'commands1.txt' in wrapper.read_text()
    assert len(commands1.read_text().splitlines()) == 6
    execution_files = [rel(p) for p in all_new if p.suffix in (".log", ".out", ".err") or p.name.startswith(".command")]
    assert not execution_files
    unique_context = {(r["cell_line"], r["chromosome"], r["antitarget_start"]): r for r in context}
    assert len(unique_context) == 92
    summary = dict(script="scripts/27_cnv_coverage_provenance.py", bundle=rel(bundle),
        input_payload_files=108, input_total_files=109,
        identity_categories=dict(Counter(r["kind"] for r in identities)),
        all_compared_identities_match=True, normal_historical_inventory_links=230,
        commands_unchanged=True, new_execution_logs_or_version_records=False,
        current_wrapper_array="1-27", current_wrapper_selected_command_lines=6,
        models=23, distinct_positive_antitarget_coordinates=4, positive_antitarget_records=92,
        positive_antitarget_records_overlapping_targets=92, antitarget_segment_context_rows=len(context),
        distinct_overlapping_segments=len(segment_diagnostics),
        distinct_contributed_segments=sum(r["contributing_antitarget_bins"] > 0 for r in segment_diagnostics),
        all_CNS_weighted_means_reconstructed=len(segment_mean_errors),
        max_absolute_CNS_mean_reconstruction_error=max(segment_mean_errors),
        autosomal_segment_labels_changed_at_fixed_center=sum(r["autosomal_gain_loss_label_changes_at_fixed_center"] for r in segment_diagnostics),
        cnr_antitarget_log2=bounds([r["cnr_log2"] for r in unique_context.values()]),
        cnr_antitarget_weight=bounds([r["cnr_weight"] for r in unique_context.values()]),
        target_bin_weight_medians=bounds([r["target_weight_median"] for r in model_summaries]),
        affected_segment_log2=bounds([r["segment_log2"] for r in segment_diagnostics]),
        affected_segment_probes=bounds([r["segment_probes"] for r in segment_diagnostics]),
        antitarget_fraction_contributing_weight=bounds([r["antitarget_fraction_contributing_weight"] for r in segment_diagnostics]),
        raw_weighted_mean_change_removing_antitargets=bounds([r["raw_weighted_mean_change_removing_antitargets"] for r in segment_diagnostics]),
        archived_segment_minus_target_only_raw_mean=bounds([r["archived_segment_minus_target_only_raw_mean"] for r in segment_diagnostics]),
        scientific_outputs_changed=False,
        diagnostic_limit="Fixed-membership means partition ordered positive-depth CNR rows by CNS probe counts and reproduce all archived means within printed precision. Removing antitarget contributions here does not rerun CBS or establish alternative breakpoints/copy-number calls.")
    (out / "provenance_assessment.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
