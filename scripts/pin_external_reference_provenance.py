#!/usr/bin/env python3
"""Pin existing external reference bytes; never downloads or mutates the inputs.

The Figshare metadata was fetched from the exact public /versions/1 endpoint on
2026-09-05. Historical snapshot dates below come from project documentation, not
filesystem timestamps. Cellosaurus entry versions are not database releases.
"""
from pathlib import Path
import csv
import hashlib
import json
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / 'output/external'
REF = ROOT / 'data/reference/depmap_24Q4_figshare_v1_metadata.json'
META_VERIFIED_DATE = '2026-09-05'
meta = json.loads(REF.read_text())
assert meta['id'] == 27993248 and meta['version'] == 1
assert meta['title'] == 'DepMap 24Q4 Public'
files = {x['name']: x for x in meta['files']}
primary = ['Model.csv', 'OmicsExpressionProteinCodingGenesTPMLogp1.csv',
           'OmicsSomaticMutationsMatrixHotspot.csv', 'OmicsSomaticMutationsMatrixDamaging.csv']
derived = {
    'ovarian_ach.txt': 'Model.csv',
    'depmap_ovarian_models.csv': 'Model.csv',
    'depmap_expr_ovarian.csv': 'OmicsExpressionProteinCodingGenesTPMLogp1.csv',
    'depmap_expr_overlap5.csv': 'depmap_expr_ovarian.csv',
    'depmap_hotspot_overlap5.csv': 'OmicsSomaticMutationsMatrixHotspot.csv',
    'depmap_damaging_overlap5.csv': 'OmicsSomaticMutationsMatrixDamaging.csv',
}

def checksums(p):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            md5.update(chunk); sha.update(chunk)
    return {'bytes': p.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}

rows = []
for name in primary:
    f = EXT / name; ref = files[name]; h = checksums(f)
    assert h['md5'] == ref['computed_md5'] and h['bytes'] == ref['size'], name
    rows.append(dict(file=str(f.relative_to(ROOT)), resource='DepMap', kind='primary_download',
        **h, declared_release='DepMap Public 24Q4', verified_release=meta['title'],
        verified_release_version=meta['version'], doi=meta['doi'],
        source_url=ref['download_url'], source_file_id=ref['id'], source_md5=ref['computed_md5'],
        provenance_status='local bytes match public versioned Figshare metadata MD5 and size',
        metadata_verified_date=META_VERIFIED_DATE, documented_snapshot_date='2026-07-23',
        date_evidence='scripts/fetch_external_data.R:67-68; historical note, not HTTP retrieval record',
        record_versions=''))
for name, parent in derived.items():
    f = EXT / name
    rows.append(dict(file=str(f.relative_to(ROOT)), resource='DepMap', kind='derived_subset',
        **checksums(f), declared_release='DepMap Public 24Q4', verified_release='',
        verified_release_version='', doi=meta['doi'], source_url=parent, source_file_id='',
        source_md5='', provenance_status='derived file pinned locally; parent release independently verified',
        metadata_verified_date=META_VERIFIED_DATE, documented_snapshot_date='', date_evidence='',
        record_versions=''))
cellosaurus_details = []
for f in sorted((EXT / 'cellosaurus').glob('*.json')):
    d = json.loads(f.read_text())['Cellosaurus']
    records = []
    for entry in d.get('cell-line-list', []):
        acc = next((a['value'] for a in entry.get('accession-list', []) if a['type'] == 'primary'), '')
        records.append({'accession': acc, 'entry_version': entry.get('entry-version'),
                        'record_last_updated': entry.get('last-updated')})
    if f.name.startswith('search_'):
        url = 'https://api.cellosaurus.org/search/cell-line?q=' + quote(f.stem[7:], safe='') + '&format=json'
    else:
        url = 'https://api.cellosaurus.org/cell-line/' + f.stem.upper() + '?format=json'
    rows.append(dict(file=str(f.relative_to(ROOT)), resource='Cellosaurus', kind='cached_API_response',
        **checksums(f), declared_release='', verified_release='', verified_release_version='',
        doi='', source_url=url, source_file_id='', source_md5='',
        provenance_status='cached response pinned; database release and original HTTP retrieval timestamp absent',
        metadata_verified_date='', documented_snapshot_date='2026-07-23',
        date_evidence='reports/review_code.md:333; scripts/18_external_benchmarking.R dated 2026-07-23; documentary date only',
        record_versions=json.dumps(records, separators=(',', ':'))))
    cellosaurus_details.append({'file': str(f.relative_to(ROOT)), 'records': records})
assert len(cellosaurus_details) == 44
assert sum(Path(x['file']).name.startswith('search_') for x in cellosaurus_details) == 42
f = EXT / 'cellosaurus/names.txt'
rows.append(dict(file=str(f.relative_to(ROOT)), resource='Cellosaurus', kind='query_name_list',
    **checksums(f), declared_release='', verified_release='', verified_release_version='', doi='',
    source_url='', source_file_id='', source_md5='', provenance_status='local query list pinned',
    metadata_verified_date='', documented_snapshot_date='', date_evidence='', record_versions=''))

out = ROOT / 'output/external_reference_provenance.csv'
with out.open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
summary = {
    'audit_date': META_VERIFIED_DATE,
    'depmap': {'declared_release': 'DepMap Public 24Q4', 'verified_title': meta['title'],
        'verified_version': meta['version'], 'verified_doi': meta['doi'],
        'source_metadata_url': 'https://api.figshare.com/v2/articles/27993248/versions/1',
        'source_metadata_path': str(REF.relative_to(ROOT)), 'source_metadata_checksums': checksums(REF),
        'four_primary_files_match_public_md5_and_size': True,
        'historical_snapshot_date_documented': '2026-07-23',
        'historical_date_evidence': 'scripts/fetch_external_data.R:67-68',
        'historical_http_retrieval_record_available': False},
    'cellosaurus': {'cached_search_responses': 42, 'cached_accession_responses': 2,
        'historical_snapshot_date_documented': '2026-07-23',
        'historical_date_evidence': 'reports/review_code.md:333; dated analysis in scripts/18_external_benchmarking.R',
        'historical_http_retrieval_record_available': False, 'database_release': None,
        'limit': 'Individual entry versions and last-updated dates are not a database release or retrieval date.',
        'record_metadata': cellosaurus_details},
    'file_manifest': str(out.relative_to(ROOT)), 'files': rows,
}
(ROOT / 'output/external_reference_provenance.json').write_text(json.dumps(summary, indent=2) + '\n')
print(f'Pinned {len(rows)} files; all four primary DepMap files match verified Figshare v1 metadata.')
print('Cellosaurus database release not recoverable; July 23, 2026 is a documentary snapshot date only.')
