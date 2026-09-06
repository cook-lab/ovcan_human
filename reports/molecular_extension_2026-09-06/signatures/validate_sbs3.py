#!/usr/bin/env python3
"""Validate standalone SBS3 results without R or access to private raw files."""
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

DEST = Path(__file__).resolve().parent
ROOT = DEST.parents[2]


def rows(name):
    with (DEST / name).open() as handle:
        return list(csv.DictReader(handle))


def close(a, b, tol=1e-10):
    assert math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol), (a, b)


def quantile(values, p):
    x = sorted(values)
    i = (len(x) - 1) * p
    lo, hi = math.floor(i), math.ceil(i)
    return x[lo] + (x[hi] - x[lo]) * (i - lo)


def main():
    manifest = json.loads((DEST / 'run_manifest.json').read_text())
    n_boot = manifest['n_boot']
    summary, boot = rows('sbs3_sensitivity_summary.csv'), rows('sbs3_bootstrap_replicates.csv')
    assert len(summary) == 276
    assert len(boot) == 184 * n_boot
    assert len({r['cell_line'] for r in summary}) == 23
    assert len({r['patient_id'] for r in summary}) == 16
    grouped = defaultdict(list)
    for r in boot:
        key = tuple(r[k] for k in ('cell_line', 'substrate', 'reference'))
        grouped[key].append(r)
        close(r['sbs3_fraction'], float(r['sbs3_count']) / float(r['total_fitted']))
    for s in summary:
        key = tuple(s[k] for k in ('cell_line', 'substrate', 'reference'))
        close(s['sbs3_fraction'], float(s['sbs3_count']) / float(s['total_fitted']))
        assert 0 <= float(s['sbs3_fraction']) <= 1
        assert float(s['nnls_without_sbs3_delta_sse_per_n_squared']) >= -1e-12
        if int(s['n_boots']) == 0:
            assert s['substrate'] == 'rare' and key not in grouped
            continue
        records = grouped[key]
        assert len(records) == n_boot
        assert sorted(int(r['replicate']) for r in records) == list(range(1, n_boot + 1))
        values = [float(r['sbs3_fraction']) for r in records]
        close(s['boot_fraction_lo95'], quantile(values, .025))
        close(s['boot_fraction_hi95'], quantile(values, .975))
        close(s['boot_fraction_median'], statistics.median(values))
        close(s['boot_selected_fraction'], sum(v > 0 for v in values) / n_boot)
    exposures = rows('all_signature_point_exposures.csv')
    assert len(exposures) == 23 * 3 * 2 * (60 + 22)
    fraction_sums = defaultdict(float)
    for r in exposures:
        fraction_sums[tuple(r[k] for k in ('cell_line', 'substrate', 'reference'))] += float(r['fraction'])
    assert len(fraction_sums) == 276
    for value in fraction_sums.values():
        close(value, 1)
    for s in rows('source_manifest.csv'):
        assert hashlib.sha256((ROOT / s['path']).read_bytes()).hexdigest() == s['sha256']
    assert hashlib.sha256((ROOT / 'scripts/44_wes_sbs3_sensitivity.R').read_bytes()).hexdigest() == manifest['script_sha256']
    # Existing sparse canonical point exposures: absent SBS3 rows are zero.
    with (ROOT / 'output/wes_signature_refit_exposures.csv').open() as handle:
        old = {(r['cell_line'], 'genome_' + r['reference_set']): float(r['rel_contribution'])
               for r in csv.DictReader(handle) if r['signature'] == 'SBS3'}
    legacy_errors = []
    for s in summary:
        if s['substrate'] == 'baseline' and s['opportunity_model'] == 'genome':
            previous = old.get((s['cell_line'], s['reference']), 0)
            close(previous, s['sbs3_fraction'])
            legacy_errors.append(abs(previous - float(s['sbs3_fraction'])))
    assert len(legacy_errors) == 46
    # Patient means are descriptive means across available models, not replicates.
    patients = rows('patient_descriptive_sensitivity.csv')
    assert len(patients) == 16 * 3 * 4
    for p in patients:
        models = [s for s in summary if all(s[k] == p[k] for k in ('patient_id', 'substrate', 'reference'))]
        assert len(models) == int(p['n_models'])
        close(p['mean_model_sbs3_fraction'], statistics.mean(float(s['sbs3_fraction']) for s in models))
    diagnostic_summary = []
    for substrate in ('baseline', 'rare', 'rare_read_supported'):
        for reference in ('genome_full', 'genome_restricted', 'target_full', 'target_restricted'):
            group = [s for s in summary if s['substrate'] == substrate and s['reference'] == reference]
            diagnostic_summary.append(dict(substrate=substrate, reference=reference, n_models=len(group),
                n_patients=len({s['patient_id'] for s in group}),
                n_models_nonzero_point_sbs3=sum(float(s['sbs3_fraction']) > 0 for s in group),
                n_models_boot_lower_gt_zero=sum(float(s['boot_fraction_lo95']) > 0 for s in group) if substrate != 'rare' else '',
                maximum_selection_fraction=max(float(s['boot_selected_fraction']) for s in group) if substrate != 'rare' else '',
                maximum_nnls_cosine_loss_without_sbs3=max(float(s['nnls_with_sbs3_cosine']) - float(s['nnls_without_sbs3_cosine']) for s in group)))
    with (DEST / 'dictionary_summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(diagnostic_summary[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(diagnostic_summary)
    out = dict(status='passed', models=23, patients=16, point_fits=276,
               bootstrap_fits=184, bootstrap_replicates=len(boot),
               canonical_genome_point_fits_reproduced=46,
               maximum_absolute_legacy_fraction_difference=max(legacy_errors),
               checks=['bootstrap percentile ranges and fractions', 'explicit complete reference grid',
                       'source and script hashes', 'canonical point-fit reproduction',
                       'nested NNLS residual inequality', 'patient descriptive means'])
    (DEST / 'validation.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
