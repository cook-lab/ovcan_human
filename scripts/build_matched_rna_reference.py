#!/usr/bin/env python3
"""Rebuild the exact transcript-to-gene map from the pinned Ensembl 93 cDNA FASTA.

The FASTA is an immutable reference input, not downloaded during analysis. All
versioned IDs and target lengths are independently checked in 01_rna_load_qc.R.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json

root = Path(__file__).resolve().parents[1]
ref = root / 'data/reference'
fasta = ref / 'Homo_sapiens.GRCh38.cdna.all.rel93.fa.gz'
expected_sha256 = '682a8cb58aff6084ad9053476c6bfc45d3afa20b32930c37807c37919a840f0f'
actual = hashlib.sha256(fasta.read_bytes()).hexdigest()
if actual != expected_sha256:
    raise ValueError(f'Pinned cDNA FASTA checksum mismatch: {actual}')
records = []
rec = None
with gzip.open(fasta, 'rt') as handle:
    for line in handle:
        if line.startswith('>'):
            if rec:
                records.append(rec)
            fields = line.strip().split()
            attrs = dict(field.split(':', 1) for field in fields[1:] if ':' in field)
            region = fields[2].split(':')[2]
            rec = {
                'ensembl_transcript_id': fields[0][1:].split('.')[0],
                'ensembl_gene_id': attrs['gene'].split('.')[0],
                'external_gene_name': attrs.get('gene_symbol', ''),
                'gene_biotype': attrs['gene_biotype'],
                'transcript_id_versioned': fields[0][1:],
                'gene_id_versioned': attrs['gene'],
                'seq_region': region,
                'primary_assembly': not region.startswith('CHR_'),
                'transcript_length': 0,
            }
        else:
            rec['transcript_length'] += len(line.strip())
    if rec:
        records.append(rec)
map_path = ref / 'tx2gene_ensembl_rel93.csv'
with map_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(records)
map_sha256 = hashlib.sha256(map_path.read_bytes()).hexdigest()
if map_sha256 != '57372b5a8933360f9fd89cfd394426cdb48250958f8908ce3092f0413ae2fb6b':
    raise ValueError(f'Transcript map changed unexpectedly: {map_sha256}')
provenance = {
    'reference': 'Ensembl release 93, Homo sapiens GRCh38 cDNA all',
    'original_source_url': 'https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz',
    'recovered_from': '/Users/dpcook/Analysis/test/data/reference/Homo_sapiens.GRCh38.cdna.all.rel93.fa.gz',
    'recovered_on': '2026-09-05',
    'source_verification': 'Official Ensembl release-93 CHECKSUMS retrieved 2026-09-05: BSD sum 19566, 65198 blocks for Homo_sapiens.GRCh38.cdna.all.fa.gz; recovered FASTA matches both. All versioned target IDs and lengths independently match every archived library. Original binary index sequence checksum remains unavailable.',
    'official_checksum_url': 'https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/CHECKSUMS',
    'official_bsd_sum': '19566 65198',
    'fasta_sha256': actual,
    'map_sha256': map_sha256,
    'map_md5': hashlib.md5(map_path.read_bytes()).hexdigest(),
    'map_builder': 'scripts/build_matched_rna_reference.py',
    'n_transcripts': len(records),
    'n_primary_assembly_transcripts': sum(r['primary_assembly'] for r in records),
    'n_genes': len({r['ensembl_gene_id'] for r in records}),
    'gene_policy': 'Retain all primary and alternative reference loci as annotated Ensembl gene rows. Symbol analyses sum TPM over gene rows with identical symbols.',
    'supersedes': 'Ensembl release-105 transcript map, which excluded 3529 index targets and 1.60-3.41% of library TPM.',
}
(ref / 'rna_reference_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(f'Wrote {len(records):,} exact transcript mappings; SHA256 {map_sha256}')
