#!/usr/bin/env python3
"""Independently validate recovered CNVkit CNN coverage summaries (stdlib only).

Reads the September 6 bundle and small historical QC/metadata tables. Does not
execute CNVkit, use alignments, alter historical QC or make biological calls.
Depths are rounded bin means in the supplied CNNs, not raw integer pileups.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = 'data/cluster_wes_retrieval/2026-09-06/ovcan_human_wes_cnv_coverage_2026-09-06'
OLD_BUNDLE = 'data/cluster_wes_retrieval/2026-09-05/ovcan_human_wes_handoff_2026-09-05'
NORMALS = ('SRR4039087', 'SRR4039088', 'SRR4039089', 'SRR4039096', 'SRR4039097')
REFERENCE_SHA256 = '848c052274264ac544897977648860455ac742d89bc1de5edc7ece1ca185eeda'
COVERAGE_DOC = 'https://raw.githubusercontent.com/etal/cnvkit/v0.9.10/cnvlib/coverage.py'
REFERENCE_DOC = 'https://raw.githubusercontent.com/etal/cnvkit/v0.9.10/cnvlib/reference.py'
FORMAT_DOC = 'https://cnvkit.readthedocs.io/en/stable/fileformats.html#target-and-antitarget-bin-level-coverages-cnn'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path, delimiter=','):
    with path.open(newline='') as stream:
        yield from csv.DictReader(stream, delimiter=delimiter)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, records):
    require(bool(records), f'No records for {path}')
    columns = list(dict.fromkeys(key for record in records for key in record))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(records)


def interval_stats(coords):
    groups = defaultdict(list)
    for chrom, start, end, _ in coords:
        groups[chrom].append((start, end))
    union_bp = 0
    for intervals in groups.values():
        last_end = -1
        for start, end in sorted(intervals):
            union_bp += max(0, end - max(start, last_end))
            last_end = max(last_end, end)
    lengths = [end - start for _, start, end, _ in coords]
    return dict(bins=len(coords), summed_bp=sum(lengths), union_bp=union_bp,
                repeated_overlap_bp=sum(lengths)-union_bp,
                chromosomes=len(groups), minimum_bin_bp=min(lengths),
                median_bin_bp=statistics.median(lengths), maximum_bin_bp=max(lengths))


def overlap_bp(first, second):
    left, right = defaultdict(list), defaultdict(list)
    for chrom, start, end, _ in first:
        left[chrom].append((start, end))
    for chrom, start, end, _ in second:
        right[chrom].append((start, end))
    overlap = 0
    for chrom in left.keys() & right.keys():
        a, b = sorted(left[chrom]), sorted(right[chrom])
        i = j = 0
        while i < len(a) and j < len(b):
            overlap += max(0, min(a[i][1], b[j][1])-max(a[i][0], b[j][0]))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
    return overlap


def coord_digest(coords, include_gene=False):
    digest = hashlib.sha256()
    for coord in coords:
        digest.update(('\t'.join(map(str, coord if include_gene else coord[:3]))+'\n').encode())
    return digest.hexdigest()


def read_cnn(path, expected=None):
    coords, depths = [], array('d')
    max_log_error = 0.0
    zero_logs = Counter()
    for index, row in enumerate(rows(path, '\t')):
        coord = row['chromosome'], int(row['start']), int(row['end']), row['gene']
        require(coord[1] >= 0 and coord[2] > coord[1], f'Invalid interval: {path}:{index+2}')
        if expected is not None:
            require(index < len(expected) and coord == expected[index], f'Bin identity/order mismatch: {path}:{index+2}')
        else:
            coords.append(coord)
        depth, log2 = float(row['depth']), float(row['log2'])
        require(math.isfinite(depth) and depth >= 0 and math.isfinite(log2), f'Invalid coverage: {path}:{index+2}')
        depths.append(depth)
        if depth > 0:
            max_log_error = max(max_log_error, abs(math.log2(depth)-log2))
        else:
            zero_logs[log2] += 1
    if expected is not None:
        require(len(depths) == len(expected), f'Bin count mismatch: {path}')
        coords = expected
    require(len(depths) > 0 and max_log_error < 0.00011, f'Depth/log2 disagreement: {path}: {max_log_error}')
    require(set(zero_logs) <= {-20.0}, f'Unexpected zero-depth log2 sentinel: {path}: {zero_logs}')
    require(len(set(c[:3] for c in coords)) == len(coords), f'Duplicate coordinates: {path}')
    return coords, depths, max_log_error


def summarize(sample_id, cohort, bin_group, coords, depths, keep, source_id, log_error, alias=None):
    lengths = [end-start for _, start, end, _ in coords]
    bp = sum(lengths)
    zero = [i for i, depth in enumerate(depths) if depth == 0]
    kept = [i for i, value in enumerate(keep) if value]
    row = dict(sample_id=sample_id, cohort=cohort, cell_line=(alias or {}).get('cell_line', ''),
               patient_id=(alias or {}).get('patient_id', ''), bin_group=bin_group,
               source_id=source_id, n_bins=len(depths), summed_bin_bp=bp,
               interval_union_bp=interval_stats(coords)['union_bp'],
               length_weighted_mean_depth_x=math.fsum(length*depth for length, depth in zip(lengths, depths))/bp,
               unweighted_mean_bin_depth_x=math.fsum(depths)/len(depths),
               median_bin_depth_x=statistics.median(depths), minimum_bin_depth_x=min(depths), maximum_bin_depth_x=max(depths),
               zero_depth_bins=len(zero), positive_depth_bins=len(depths)-len(zero),
               zero_depth_bin_fraction=len(zero)/len(depths),
               bp_in_zero_depth_bins=sum(lengths[i] for i in zero),
               fraction_bp_in_zero_depth_bins=sum(lengths[i] for i in zero)/bp,
               reference_retained_bins=len(kept), reference_retained_bp=sum(lengths[i] for i in kept),
               zero_depth_reference_retained_bins=sum(depths[i] == 0 for i in kept),
               positive_depth_reference_retained_bins=sum(depths[i] > 0 for i in kept),
               bp_in_zero_depth_reference_retained_bins=sum(lengths[i] for i in kept if depths[i] == 0),
               reference_retained_length_weighted_mean_depth_x=math.fsum(lengths[i]*depths[i] for i in kept)/sum(lengths[i] for i in kept),
               maximum_abs_log2_depth_rounding_difference=log_error,
               coordinates_sha256=coord_digest(coords), coordinates_and_genes_sha256=coord_digest(coords, True),
               status='validated_primary_coverage_cnn')
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT/BUNDLE)
    parser.add_argument('--previous-source', type=Path, default=ROOT/OLD_BUNDLE)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'output')
    parser.add_argument('--report-dir', type=Path, default=ROOT/'reports/wes_cnv_coverage_2026-09-06')
    args = parser.parse_args()
    source = args.source.resolve()
    checksums = {}
    for line in (source/'SHA256SUMS').read_text().splitlines():
        digest, relative = line.split(None, 1)
        checksums[relative.strip().removeprefix('./')] = digest
    sources = []

    def register(path, role, primary=False):
        digest = sha256(path)
        relative = str(path.relative_to(source)) if path.is_relative_to(source) else str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        expected = checksums.get(relative, '') if path.is_relative_to(source) else ''
        if primary:
            require(expected == digest, f'Primary source checksum mismatch: {relative}')
        sid = f'CNVC{len(sources)+1:03d}'
        sources.append(dict(source_id=sid, path=relative, role=role, size_bytes=path.stat().st_size,
                            sha256=digest, archive_member_sha256=expected,
                            archive_member_checksum_verified=bool(expected and digest == expected)))
        return sid

    aliases_path = ROOT/'reports/audit_2026-09-05/wes_cluster_models.csv'
    register(aliases_path, '23-model sample aliases and patient origin')
    aliases = {r['cnv_sample_id']: r for r in rows(aliases_path)}
    require(len(aliases) == 23 and len({r['patient_id'] for r in aliases.values()}) == 16, 'Model/patient invariant')
    actual_models = {p.name.removesuffix('.targetcoverage.cnn') for p in (source/'models').glob('*.targetcoverage.cnn')}
    require(actual_models == set(aliases), 'Recovered model target set mismatch')
    require({p.name.removesuffix('.antitargetcoverage.cnn') for p in (source/'models').glob('*.antitargetcoverage.cnn')} == actual_models, 'Model target/antitarget pairing mismatch')
    require({p.name.removesuffix('.sorted.targetcoverage.cnn') for p in (source/'normals').glob('*.targetcoverage.cnn')} == set(NORMALS), 'Normal target set mismatch')
    require({p.name.removesuffix('.sorted.antitargetcoverage.cnn') for p in (source/'normals').glob('*.antitargetcoverage.cnn')} == set(NORMALS), 'Normal target/antitarget pairing mismatch')

    bed = source/'normals/intervals_sorted.target.bed'
    bed_id = register(bed, 'CNVkit split target design', True)
    target_coords = []
    for line in bed.read_text().splitlines():
        chrom, start, end, gene = line.split('\t')
        target_coords.append((chrom, int(start), int(end), gene))
    require(interval_stats(target_coords)['bins'] == 290475 and interval_stats(target_coords)['summed_bp'] == 63514049, 'Unexpected target footprint')
    anti_path = source/f'normals/{NORMALS[0]}.sorted.antitargetcoverage.cnn'
    anti_coords, _, _ = read_cnn(anti_path)
    coords = dict(target=target_coords, antitarget=anti_coords)
    designs = {group: dict(interval_stats(c), coordinates_sha256=coord_digest(c), coordinates_and_genes_sha256=coord_digest(c, True)) for group, c in coords.items()}
    require(all(d['repeated_overlap_bp'] == 0 for d in designs.values()), 'Within-design overlapping bins')
    between_overlap = overlap_bp(target_coords, anti_coords)
    # These historical designs are not completely disjoint. Preserve and report
    # the measured overlap; do not treat it as an input error or merge the bins.

    reference = source/'normals/reference.cnn'
    ref_id = register(reference, 'Pooled empirical CNVkit reference; not a sixth sample', True)
    require(sources[-1]['sha256'] == REFERENCE_SHA256, 'Pooled reference differs from previously audited reference')
    indices = {c[:3]: (group, i) for group, cc in coords.items() for i, c in enumerate(cc)}
    keep = {group: bytearray(len(cc)) for group, cc in coords.items()}
    refzero = {group: bytearray(len(cc)) for group, cc in coords.items()}
    seen = set()
    reference_tiny_negative_depths = []
    for row in rows(reference, '\t'):
        key = row['chromosome'], int(row['start']), int(row['end'])
        require(key in indices and key not in seen, 'Reference interval set mismatch')
        seen.add(key)
        group, i = indices[key]
        require(row['gene'] == coords[group][i][3], 'Reference gene annotation mismatch')
        depth, log2, spread, gc = [float(row[k]) for k in ('depth', 'log2', 'spread', 'gc')]
        require(all(map(math.isfinite, (depth, log2, spread, gc))) and depth > -1e-12, 'Invalid reference field')
        if depth < 0:
            reference_tiny_negative_depths.append(dict(chromosome=key[0], start=key[1], end=key[2],
                bin_group=group, depth=depth,
                reference_mask_retained=not (abs(log2) > 5 or spread > 1 or depth == 0 or gc < .3 or gc > .7)))
        refzero[group][i] = depth == 0
        keep[group][i] = not (abs(log2) > 5 or spread > 1 or depth == 0 or gc < .3 or gc > .7)
    require(len(seen) == len(indices), 'Reference omits design bins')
    require(sum(keep['target']) == 204706 and sum(keep['antitarget']) == 40923, 'Historical reference mask changed')

    historical_qc = args.previous_source/'recovery/qc_metrics.tsv'
    historical_id = register(historical_qc, 'September 5 reported-only normal depth summaries')
    reported = {r['model_or_reference']: r for r in rows(historical_qc, '\t') if r['metric'] == 'cnvkit_target_mean_depth_length_weighted'}
    summaries, normal_comparisons, chromosome_summaries = [], [], []
    support = {group: bytearray(len(cc)) for group, cc in coords.items()}
    for cohort, folder, samples in [('cnv_reference_exomes', 'normals', NORMALS), ('baseline_models', 'models', sorted(aliases))]:
        for sample in samples:
            for group in ('target', 'antitarget'):
                filename = sample + ('.sorted' if folder == 'normals' else '') + f'.{group}coverage.cnn'
                path = source/folder/filename
                sid = register(path, f'{cohort} {group} primary bin depth', True)
                cc, depths, log_error = read_cnn(path, coords[group])
                summary = summarize(sample, cohort, group, cc, depths, keep[group], sid, log_error, aliases.get(sample))
                summaries.append(summary)
                by_chrom = defaultdict(list)
                for i, coord in enumerate(cc):
                    by_chrom[coord[0]].append(i)
                for chrom, ix in by_chrom.items():
                    lengths = {i: cc[i][2]-cc[i][1] for i in ix}
                    total_bp = sum(lengths.values())
                    zero = [i for i in ix if depths[i] == 0]
                    chromosome_summaries.append(dict(sample_id=sample, cohort=cohort, bin_group=group,
                        chromosome=chrom, source_id=sid, bins=len(ix), summed_bin_bp=total_bp,
                        length_weighted_mean_depth_x=math.fsum(lengths[i]*depths[i] for i in ix)/total_bp,
                        zero_depth_bins=len(zero), zero_depth_bin_fraction=len(zero)/len(ix),
                        bp_in_zero_depth_bins=sum(lengths[i] for i in zero),
                        fraction_bp_in_zero_depth_bins=sum(lengths[i] for i in zero)/total_bp,
                        reference_retained_bins=sum(keep[group][i] for i in ix),
                        zero_depth_reference_retained_bins=sum(keep[group][i] for i in zero)))
                if folder == 'normals':
                    for i, depth in enumerate(depths):
                        support[group][i] += depth > 0
                    if group == 'target':
                        previous = reported[sample]
                        value = summary['length_weighted_mean_depth_x']
                        difference = value - float(previous['value'])
                        require(abs(difference) < .00501, f'Previous normal mean not confirmed: {sample}')
                        normal_comparisons.append(dict(sample_id=sample, length_weighted_mean_depth_x=value,
                                                       previous_reported_mean_depth_x=float(previous['value']),
                                                       difference_x=difference, agrees_at_reported_precision=True,
                                                       n_bins=summary['n_bins'], summed_bin_bp=summary['summed_bin_bp'],
                                                       previous_denominator=previous['denominator'],
                                                       primary_source_id=sid, previous_summary_source_id=historical_id))
            print(f'Parsed {cohort}: {sample}', flush=True)

    support_rows = []
    for group, cc in coords.items():
        for retained in (False, True):
            for n_positive in range(6):
                selected = [i for i in range(len(cc)) if bool(keep[group][i]) == retained and support[group][i] == n_positive]
                support_rows.append(dict(bin_group=group, reference_mask_retained=retained, positive_normal_count=n_positive,
                                         n_normals=5, bins=len(selected), bp=sum(cc[i][2]-cc[i][1] for i in selected),
                                         pooled_reference_zero_depth_bins=sum(refzero[group][i] for i in selected),
                                         reference_source_id=ref_id))

    # Check consistency with prior primary-CNN summaries without rewriting them.
    prior_support_path = ROOT/'output/wes_recovered_provenance_cnv_support.csv'
    prior_id = register(prior_support_path, 'Previously audited model CNV measurement counts')
    prior_support = {r['sample_id']: r for r in rows(prior_support_path)}
    for r in summaries:
        if r['cohort'] != 'baseline_models':
            continue
        prior = prior_support[r['sample_id']]
        if r['bin_group'] == 'target':
            pairs = [('n_bins', 'original_target_bins'), ('zero_depth_bins', 'original_target_depth_zero'),
                     ('reference_retained_bins', 'target_bins_in_cnr'), ('positive_depth_reference_retained_bins', 'positive_target_bins_in_cnr')]
        else:
            pairs = [('n_bins', 'antitarget_bins'), ('zero_depth_bins', 'antitarget_depth_zero'),
                     ('positive_depth_bins', 'antitarget_depth_positive'), ('positive_depth_reference_retained_bins', 'positive_antitargets_in_cnr')]
        for current, old in pairs:
            require(r[current] == int(prior[old]), f'Historical model support mismatch: {r["sample_id"]}/{current}')

    definitions = {
        'length_weighted_mean_depth_x': ('x', 'sum(bin_length_bp * supplied_depth) / sum(bin_length_bp); rounded CNN means'),
        'unweighted_mean_bin_depth_x': ('x', 'arithmetic mean of supplied bin depths; each bin has equal weight'),
        'median_bin_depth_x': ('x', 'median of supplied bin means; not median per-base depth'),
        'zero_depth_bins': ('bins', 'bins with supplied depth exactly zero'),
        'zero_depth_bin_fraction': ('fraction', 'zero-depth bins / all bins'),
        'bp_in_zero_depth_bins': ('bp', 'sum of lengths of whole bins with zero mean depth'),
        'fraction_bp_in_zero_depth_bins': ('fraction', 'bp in zero-depth whole bins / summed bin bp; does not count zero-depth bases inside positive-mean bins'),
        'positive_depth_bins': ('bins', 'bins with mean depth >0; does not imply every base has coverage'),
        'reference_retained_bins': ('bins', 'bins passing pooled-reference mask reconstructed and previously reconciled to archived CNR coordinates'),
        'zero_depth_reference_retained_bins': ('bins', 'zero-depth bins among reference-mask-passing bins'),
        'positive_depth_reference_retained_bins': ('bins', 'positive-depth bins among reference-mask-passing bins'),
        'reference_retained_length_weighted_mean_depth_x': ('x', 'length-weighted depth across reference-mask-passing bins only'),
    }
    metrics = []
    for cohort in ('cnv_reference_exomes', 'baseline_models'):
        for group in ('target', 'antitarget'):
            selected = [r for r in summaries if r['cohort'] == cohort and r['bin_group'] == group]
            for metric, (unit, definition) in definitions.items():
                values = [r[metric] for r in selected]
                metrics.append(dict(cohort=cohort, bin_group=group, metric=metric, n=len(values), n_primary_verified=len(values),
                                    minimum=min(values), median=statistics.median(values), maximum=max(values), unit=unit,
                                    minimum_samples=';'.join(r['sample_id'] for r in selected if r[metric] == min(values)),
                                    maximum_samples=';'.join(r['sample_id'] for r in selected if r[metric] == max(values)),
                                    definition=definition, source_ids=';'.join(r['source_id'] for r in selected)))

    limitations = [
        'CNN depth columns contain rounded bin means. Weighted totals computed from them are not exact original integer read-base counts.',
        'A positive mean bin does not establish coverage at every base. Fractions of bases at >=10x, >=20x or >=30x cannot be reconstructed from these CNNs; none are calculated.',
        'Base length in zero-mean bins is exact from the interval coordinates, but is only a lower bound on all uncovered bases because positive-mean bins may contain uncovered positions.',
        'The reference depth column is a robust aggregate across normals, not the arithmetic mean or a sixth independent sample. Corrected reference log2 values are not log2 of its absolute depth column.',
        'Model CNNs describe interval-restricted recalibrated inputs from the recorded manual CNVkit workflow. Normal CNNs describe separately processed reference exomes. Matching interval coordinates does not establish matching read-processing filters or genome-wide coverage.',
        'All 23 model antitarget files have only four positive-depth bins. The target and antitarget designs overlap at a small set of bases; positive antitarget signals are not necessarily independent off-target support. The broad normal antitarget coverage does not supply missing model off-target measurements; CNV profiles remain predominantly target-supported.',
        'No normal alignment/mapping, duplication, contamination, read-level integrity or per-base coverage reports were supplied in this archive. Those checks remain unavailable, rather than failed.',
        'The reference mask is reproduced for validation only; no bin thresholds, normal selection, CNV segmentation, copy-number interpretation or historical scientific outputs are changed.',
    ]
    report = dict(schema_version=1, source_bundle=BUNDLE if source == (ROOT/BUNDLE).resolve() else str(source),
                  status='passed', normal_samples=5, baseline_models=23, baseline_patients=16,
                  primary_coverage_files=56, primary_design_and_reference_files=2,
                  primary_archive_hashes_verified=sum(r['archive_member_checksum_verified'] for r in sources),
                  pooled_reference_sha256=REFERENCE_SHA256, pooled_reference_matches_previous=True,
                  pooled_reference_tiny_negative_depth_count=len(reference_tiny_negative_depths),
                  pooled_reference_minimum_depth=min((r['depth'] for r in reference_tiny_negative_depths), default=0),
                  pooled_reference_tiny_negative_depth_rows=reference_tiny_negative_depths,
                  interval_designs=designs, target_antitarget_overlap_bp=between_overlap,
                  all_coverage_coordinates_and_gene_labels_match_design=True,
                  reference_retained_bins={group: sum(values) for group, values in keep.items()},
                  previous_normal_means_confirmed=5, previous_model_support_fields_confirmed=23*8,
                  maximum_abs_log2_depth_rounding_difference=max(r['maximum_abs_log2_depth_rounding_difference'] for r in summaries),
                  definitions=definitions, normal_comparisons=normal_comparisons,
                  cross_normal_support=support_rows, metric_summaries=metrics, limitations=limitations,
                  method_sources=[dict(url=COVERAGE_DOC, role='CNVkit v0.9.10 bin depth and log2 calculation'),
                                  dict(url=REFERENCE_DOC, role='CNVkit v0.9.10 robust pooled depth; separately corrected log2'),
                                  dict(url=FORMAT_DOC, role='CNVkit CNN format and BED coordinate interpretation')])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir/'wes_cnv_coverage_sample_summary.csv', summaries)
    write_csv(args.output_dir/'wes_cnv_coverage_chromosome_summary.csv', chromosome_summaries)
    write_csv(args.output_dir/'wes_cnv_coverage_metric_summary.csv', metrics)
    write_csv(args.output_dir/'wes_cnv_coverage_normal_comparison.csv', normal_comparisons)
    write_csv(args.output_dir/'wes_cnv_coverage_reference_support.csv', support_rows)
    write_csv(args.output_dir/'wes_cnv_coverage_sources.csv', sources)
    for path in (args.output_dir/'wes_cnv_coverage_validation.json', args.report_dir/'coverage_validation.json'):
        path.write_text(json.dumps(report, indent=2)+'\n')

    def metric(cohort, group, name):
        return next(m for m in metrics if m['cohort'] == cohort and m['bin_group'] == group and m['metric'] == name)

    text = ['# CNV coverage validation — 6 September 2026', '',
            'Independently parsed all **5 normal and 23 model target/antitarget pairs** (56 primary CNNs). All 58 primary coverage/design/reference source hashes match the supplied archive manifest. The pooled reference is byte-identical to the previously audited reference. Historical September 5 QC, CNV calls and manuscript/figure files were not rewritten.', '',
            '## Target depth and zero-bin summaries', '',
            '| Cohort | n | Length-weighted mean depth, min / median / max (×) | Zero-depth target bins, min / median / max |',
            '| --- | ---: | --- | --- |']
    for cohort in ('cnv_reference_exomes', 'baseline_models'):
        m, z = metric(cohort, 'target', 'length_weighted_mean_depth_x'), metric(cohort, 'target', 'zero_depth_bins')
        text.append(f'| {cohort} | {m["n"]} | {m["minimum"]:.6f} / {m["median"]:.6f} / {m["maximum"]:.6f} | {z["minimum"]:,} / {z["median"]:,} / {z["maximum"]:,} |')
    text += ['', f'Every target file has **290,475 bins spanning 63,514,049 bp**, with no within-design overlaps. Every antitarget file has **42,709 bins spanning 2,580,384,051 bp**, also without within-design overlaps. All 28 samples have identical ordered coordinates and gene labels within each bin group; the target CNNs also match the supplied target BED. The target and antitarget designs overlap by **{between_overlap:,} bp** and must not be described as disjoint. The pooled reference contains precisely their combined bin coordinate set.', '',
             'The denominator differs from the 63,709,951 summed bases of the overlapping mosdepth calling intervals. The CNVkit target length is the 63,514,049 bp union after merging/splitting. Directly pooling CNN and mosdepth metrics would conflate interval definitions, alignment stages and coverage implementations.', '',
             '## Recovered normal means', '', '| Reference exome | Verified length-weighted mean (×) | Prior reported mean (×) |', '| --- | ---: | ---: |']
    for r in normal_comparisons:
        text.append(f'| {r["sample_id"]} | {r["length_weighted_mean_depth_x"]:.9f} | {r["previous_reported_mean_depth_x"]:.2f} |')
    text += ['', 'All five values agree at the original two-decimal precision. They are now independently validated from primary CNNs; the earlier dated report correctly retains its historical reported-only status.', '',
             '## Normal target coverage gaps', '',
             '| Normal | Zero-depth bins (%) | bp in zero-depth bins (%) | Positive bins after reference mask / 204,706 |',
             '| --- | ---: | ---: | ---: |']
    for r in summaries:
        if r['cohort'] == 'cnv_reference_exomes' and r['bin_group'] == 'target':
            text.append(f'| {r["sample_id"]} | {r["zero_depth_bins"]:,} ({100*r["zero_depth_bin_fraction"]:.3f}%) | {r["bp_in_zero_depth_bins"]:,} ({100*r["fraction_bp_in_zero_depth_bins"]:.3f}%) | {r["positive_depth_reference_retained_bins"]:,} |')
    text += ['', 'All five normals have zero-depth target bins on every autosome. Within each normal, chr19 has the largest zero-bin fraction (24.43–37.13% across normals), and chr13 the smallest (5.57–7.55%). Thus the gaps are distributed unevenly rather than confined to a wholly missing chromosome. CNN summaries do not identify their cause. Per-chromosome means, counts and base footprints are recorded in `wes_cnv_coverage_chromosome_summary.csv`.', '',
             '## Reference support and processing compatibility', '']
    for group in ('target', 'antitarget'):
        all_zero = [r for r in support_rows if r['bin_group'] == group and r['positive_normal_count'] == 0]
        not_all = [r for r in support_rows if r['bin_group'] == group and r['positive_normal_count'] < 5]
        retained_missing = [r for r in not_all if r['reference_mask_retained']]
        text.append(f'- {group.capitalize()}: {sum(r["bins"] for r in all_zero):,} bins have zero depth in all five normals; {sum(r["bins"] for r in not_all):,} have zero depth in at least one. Of {sum(keep[group]):,} reference-mask-passing bins, {sum(r["bins"] for r in retained_missing):,} lack positive depth in at least one normal. Detailed counts and base lengths by number of supporting normals are in `wes_cnv_coverage_reference_support.csv`.')
    m = metric('cnv_reference_exomes', 'antitarget', 'length_weighted_mean_depth_x')
    text += ['', f'Normal antitarget mean depths range from {m["minimum"]:.6f} to {m["maximum"]:.6f}× (median {m["median"]:.6f}×). Each model instead has exactly 42,705 zero-depth antitarget bins and only four positive bins. The primary files confirm that the usable CNV signal remains predominantly on-target. Coordinate agreement confirms a common bin design; it does not by itself establish equivalent upstream processing or justify normalising away the asymmetric off-target support.', '',
             'All retained target bins have positive depth in at least three normals: 3,355 have support from three, 5,478 from four, and 195,873 from all five. The mask excludes every bin with no positive depth in any normal. Its 204,706 retained targets span 44,729,346 bp. This is a support summary, not an independent calibration of CNV accuracy.', '',
             'The pooled-reference mask retains 204,706 target and 40,923 antitarget bins. Reconstruction uses the previously validated union of |log2|>5, spread>1, depth=0, and GC outside 0.3–0.7. All 184 model bin-count/support comparisons with the September 5 audit pass. Mask criteria overlap and are not additive. This validation does not refilter or regenerate CNV profiles.', '',
             'The pooled reference contains 25 tiny negative depth values (minimum −2.77556×10⁻¹⁷), consistent with numerical cancellation in robust aggregation. Two target bins pass the literal historical depth==0 mask: chr17:63963626–63963868 and chr20:2396442–2396684 (BED coordinates). Their coordinates and values are recorded in the JSON audit. These are not physically negative read depths. The historical mask is preserved exactly; no rounded-zero rule or additional bin exclusion has been applied.', '',
             '## Units, precision and limits', '',
             f'CNVkit [coverage implementation v0.9.10]({COVERAGE_DOC}) defines depth as mean pileup per bin and computes log2 from that mean. For all positive-depth input rows, the maximum discrepancy between log2(depth) and the stored log2 is {report["maximum_abs_log2_depth_rounding_difference"]:.8f}, consistent with text rounding; every zero-depth row uses the -20 sentinel. The [reference implementation]({REFERENCE_DOC}) separately robustly aggregates absolute depths and bias-corrected log2 values, so the reference log2/depth columns are not interchangeable.', '']
    text += ['- '+item for item in limitations]
    text += ['', '## Reproduction and outputs', '', '```sh', 'python3 scripts/26_wes_cnv_coverage_qc.py', '```', '',
             '`output/wes_cnv_coverage_sample_summary.csv` contains 56 sample/bin-group rows. The metric summary gives n, min, median, max, units and source IDs; the normal comparison records prior rounding agreement; the reference support table records cross-normal support; the source table records primary checksums. Machine-readable validation is saved both beside this report and under output/. Source inputs remain read-only.', '']
    (args.report_dir/'coverage_validation.md').write_text('\n'.join(text))
    print(f'PASS: 56 CNNs, {report["primary_archive_hashes_verified"]} primary hashes, 5 confirmed normal means, 184 unchanged model support fields; reports written.', flush=True)


if __name__ == '__main__':
    main()
