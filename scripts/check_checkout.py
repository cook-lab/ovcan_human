#!/usr/bin/env python3
"""Read-only validation of the committed handoff and processed release.

Python standard library only. --inputs additionally checks the local archived
input layout by size; it never searches the cluster or starts analysis.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import os
import sys

ROOT = Path(__file__).resolve().parents[1]

def rows(path, delimiter=','):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter=delimiter))

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',action='store_true',help='also check archived input paths/sizes using OVCAN_DATA or the default local archive')
    args=ap.parse_args()
    required=['README.md','AGENTS.md','CLAUDE.md','docs/PROJECT_STATUS.md',
      'docs/cluster/CLAUDE_TASK.md','docs/cluster/WES_RECOVERY.md','docs/REPRODUCIBILITY.md',
      'docs/data/archived_input_inventory.tsv','docs/cluster/evidence/manifest.csv',
      'reports/audit_2026-09-05/wes_cluster_path_hints.csv',
      'reports/audit_2026-09-05/wes_cluster_models.csv',
      'output/wes_input_manifest.csv','output/wes_vcf_header_provenance.csv',
      'output/wes_pipeline_parameters.csv','output/wes_cnvkit_target_intervals.bed',
      'docs/manuscript/v7/OvCAN_Scientific_Data_draft_v7.md','release/SHA256SUMS']
    missing=[p for p in required if not (ROOT/p).is_file()]
    if missing:
        print('FAIL missing committed handoff files:',*missing,sep='\n  ');return 1
    fam=rows(ROOT/'metadata/line_family_map.csv')
    models=rows(ROOT/'reports/audit_2026-09-05/wes_cluster_models.csv')
    hints=rows(ROOT/'reports/audit_2026-09-05/wes_cluster_path_hints.csv')
    checks=[('42 unique models',len(fam)==42 and len({r['cell_line'] for r in fam})==42),
            ('34 patients',len({r['patient_id'] for r in fam})==34),
            ('23 WES handoff models',len(models)==23 and len({r['cell_line'] for r in models})==23),
            ('23 historical BAM filenames',all(r['historical_bam_path'].endswith('.recal.bam') for r in models)),
            ('322 historical path hints',len(hints)==322)]
    for name,ok in checks:
        if not ok:print('FAIL',name);return 1
    release_count=0
    for line in (ROOT/'release/SHA256SUMS').read_text().splitlines():
        expected,relative=line.split(None,1); p=ROOT/'release'/relative.strip()
        if not p.is_file() or digest(p)!=expected:
            print('FAIL release checksum',relative);return 1
        release_count+=1
    for r in rows(ROOT/'docs/cluster/evidence/manifest.csv'):
        p=ROOT/r['evidence_path']
        if not p.is_file() or digest(p)!=r['sha256']:
            print('FAIL copied evidence checksum',r['evidence_path']);return 1
    print(f'PASS: 42 models / 34 patients; 23 WES aliases; 322 historical path hints; {release_count} release checksums; copied-command checksums.')
    print('This checks the committed handoff, not current cluster availability or full raw-input reproducibility.')
    if args.inputs:
        data=Path(os.environ.get('OVCAN_DATA',str(ROOT/'judy_archive/data'))).expanduser()
        absent=[]; changed=[]
        for r in rows(ROOT/'docs/data/archived_input_inventory.tsv','\t'):
            p=data/r['relative_path']
            if not p.is_file():absent.append(r['relative_path'])
            elif p.stat().st_size!=int(r['bytes']):changed.append(r['relative_path'])
        print(f'Archived input layout at {data}: {len(absent)} absent and {len(changed)} different-size files.')
        print('Size equality is not content or sample-identity verification; use source hashes where available.')
        for label,items in [('absent',absent),('different size',changed)]:
            for p in items[:10]:print(f'  {label}: {p}')
        return 2 if absent or changed else 0
    return 0

if __name__=='__main__':
    sys.exit(main())
