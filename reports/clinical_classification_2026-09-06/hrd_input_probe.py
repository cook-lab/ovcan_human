#!/usr/bin/env python3
"""Read-only WES eligibility/allele evidence probe; writes only this report's CSVs.

Run from any directory. OVCAN_PROJ/OVCAN_DATA support restored private inputs.
Counts describe retained records, not an unbiased germline SNP ascertainment set.
"""
import csv
import gzip
import hashlib
import os
from pathlib import Path
import re

ROOT = Path(os.environ.get('OVCAN_PROJ', Path(__file__).resolve().parents[2]))
DEST = ROOT / 'reports/clinical_classification_2026-09-06'
HRR = {'BRCA1', 'BRCA2', 'PALB2', 'RAD51C', 'RAD51D', 'BRIP1', 'BARD1', 'CDK12'}
FLAGS = ('genotype-germline-sites', 'genotype-pon-sites', 'interval-padding')


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def resolve(path):
    p = Path(path)
    if p.is_absolute():
        return p
    prefix = 'judy_archive/data/'
    if 'OVCAN_DATA' in os.environ and path.startswith(prefix):
        return Path(os.environ['OVCAN_DATA']) / path[len(prefix):]
    return ROOT / path


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def maf_key(chrom, pos, ref, alt):
    # Trim VCF shared anchor, retaining 1-based MAF coordinates.
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    if not ref:
        return (chrom, pos - 1, pos, '-', alt)
    return (chrom, pos, pos + len(ref) - 1, ref, alt or '-')


def write_csv(name, rows):
    with (DEST / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    mutations = read_csv(ROOT / 'output/wes_mutations_filtered.csv')
    wanted = {}
    for row in mutations:
        if row['Hugo_Symbol'] in HRR:
            key = (row['cell_line'], row['Chromosome'], int(row['Start_Position']),
                   int(row['End_Position']), row['Reference_Allele'], row['Tumor_Seq_Allele2'])
            assert key not in wanted
            wanted[key] = row
    summaries, evidence = [], []
    for source in read_csv(ROOT / 'output/wes_input_manifest.csv'):
        path, model = resolve(source['source_vcf']), source['cell_line']
        assert sha256(path) == source['vcf_sha256'], f'Input changed: {path}'
        stat = dict(cell_line=model, **{k: None for k in FLAGS}, annotated_records=0,
                    VCF_PASS=0, autosomal_biallelic_SNV=0, autosomal_SNV_with_AD_DP20=0,
                    AD_DP20_AF03_07=0, AD_DP20_POPAF_lt2=0, PASS_AD_DP20_POPAF_lt2=0,
                    PASS_autosomal_SNV=0, source_vcf=source['source_vcf'])
        op = gzip.open if path.suffix == '.gz' else open
        csq_fields = []
        with op(path, 'rt') as handle:
            for line in handle:
                if line.startswith('##GATKCommandLine=<ID=Mutect2,'):
                    for flag in FLAGS:
                        hit = re.search(r'--' + flag + r'\s+(\S+)', line)
                        stat[flag] = hit.group(1) if hit else 'not recorded'
                if line.startswith('##INFO=<ID=CSQ,'):
                    csq_fields = line.split('Format: ', 1)[1].split('"', 1)[0].split('|')
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                assert len(f) == 10, 'Expected one tumour sample'
                chrom, pos, _, ref, alt, _, filt, info, fmt, sample = f
                pos = int(pos)
                values = dict(zip(fmt.split(':'), sample.split(':')))
                tags = dict(item.split('=', 1) if '=' in item else (item, '') for item in info.split(';'))
                ad = [int(v) for v in values['AD'].split(',')] if values.get('AD', '.') != '.' else []
                raw_af = ad[1] / sum(ad) if len(ad) == 2 and sum(ad) else None
                stat['annotated_records'] += 1
                stat['VCF_PASS'] += filt == 'PASS'
                autosomal_snv = (chrom in {'chr' + str(i) for i in range(1, 23)}
                                  and len(ref) == len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT')
                if autosomal_snv:
                    stat['autosomal_biallelic_SNV'] += 1
                    stat['PASS_autosomal_SNV'] += filt == 'PASS'
                    if len(ad) == 2 and sum(ad) >= 20:
                        stat['autosomal_SNV_with_AD_DP20'] += 1
                        stat['AD_DP20_AF03_07'] += 0.3 <= raw_af <= 0.7
                        common = float(tags.get('POPAF', 'inf')) < 2
                        stat['AD_DP20_POPAF_lt2'] += common
                        stat['PASS_AD_DP20_POPAF_lt2'] += common and filt == 'PASS'
                key = (model,) + maf_key(chrom, pos, ref, alt)
                if key in wanted:
                    row = wanted[key]
                    annotations = [dict(zip(csq_fields, z.split('|'))) for z in tags.get('CSQ', '').split(',')]
                    annotation = next((z for z in annotations if z.get('Feature') == row['Transcript_ID']), annotations[0])
                    evidence.append(dict(cell_line=model, gene=row['Hugo_Symbol'], protein_label=row['HGVSp_Short'],
                        chromosome=chrom, vcf_pos=pos, ref=ref, alt=alt, filter=filt,
                        AD=values.get('AD'), AD_alt_fraction=raw_af, caller_AF=values.get('AF'),
                        TLOD=tags.get('TLOD'), HGVSc=annotation.get('HGVSc', ''),
                        HGVSp=annotation.get('HGVSp', ''), Existing_variation=annotation.get('Existing_variation', ''),
                        CLIN_SIG=annotation.get('CLIN_SIG', ''), source_vcf=source['source_vcf']))
        assert all(stat[k] is not None for k in FLAGS)
        summaries.append(stat)
    assert len(summaries) == 23
    assert len(evidence) == len(wanted) == 13
    write_csv('hrd_vcf_input_counts.csv', summaries)
    write_csv('hrr_candidate_vcf_evidence.csv', evidence)
    print(f'Verified {len(summaries)} source VCF hashes; matched {len(evidence)} HRR candidates.')


if __name__ == '__main__':
    main()
