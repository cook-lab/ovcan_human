# Data included in the public repository

The repository includes the current processed output tables and compact R objects, the complete validated processed release, metadata, the matched release-93 transcript map, the derived CNVkit target BED, source-command evidence and current manuscript/figures. It deliberately excludes the approximately 10 GB original input archive, duplicate result snapshots, large download caches and temporary rendering files. The local originals were preserved.

[archived_input_inventory.tsv](archived_input_inventory.tsv) records **906 files** available under the local `judy_archive/data/` at repository preparation, using paths relative to that directory and byte sizes. This is a local inventory, not a current cluster scan. It contains no new hashes; size equality alone does not verify content. The [WES input manifest](../../output/wes_input_manifest.csv) separately supplies original VCF/MAF SHA256 hashes, while [external-reference provenance](../../output/external_reference_provenance.csv) and [RNA-reference provenance](../../data/reference/rna_reference_provenance.json) identify reference inputs.

## Intentional omissions

| Omitted path | Recovery/rebuild route |
| --- | --- |
| `judy_archive/data/` | Restore from managed project storage, or point `OVCAN_DATA` to a directory with the same internal layout. Use the inventory and manifests to compare files. Includes RNA quantification outputs, proteomics workbooks and archived WES/CNV outputs; not a complete original sequencing-read deposit. |
| `output/rna_txi.rds` | Rebuilt by `scripts/01_rna_load_qc.R` from archived RNA abundance files and the matched release-93 map; processed RNA matrices are committed. |
| `output/external/OmicsExpressionProteinCodingGenesTPMLogp1.csv` and `OmicsSomaticMutationsMatrixDamaging.csv` | Obtain the exact DepMap 24Q4 Public version-1 files recorded in `output/external_reference_provenance.csv`; compare supplied MD5/size and recorded SHA256. Derived ovarian/overlap subsets and results are committed. |
| `data/reference/Homo_sapiens.GRCh38.cdna.all.rel93.fa.gz` | Exact source URL and SHA256 are in `data/reference/rna_reference_provenance.json`; the matched transcript map and reference-reconciliation results are committed. |
| `data/reference/tx2gene_ensembl_rel105.csv` | Historical mismatched map. It is intentionally excluded and must not replace the corrected release-93 map. |
| `reports/_*`, audit `before*`/`rna_before`, render PNG/PDF directories, `tmp/`, `logs/`, exploratory `figs/` | Historical or rebuildable duplicates and scratch output; current reports, decisions, current results and canonical figures are retained. Historical manifests may reference these excluded snapshots. |
| Earlier loose Word/PDF/PowerPoint drafts under `docs/` | Superseded working/meeting files; manuscripts v5–v7 and original source notebooks remain committed. |

The full `release/` package is tracked, including its 49 saved SHA256 checksums. New recovery work should not rewrite it until a reviewed data revision is deliberately made.

## Restore an input tree

From the clone root, point to the directory containing `rna_seq/`, `proteomics/`, `wes - old/` and `cnvkit wes - new/`:

```bash
export OVCAN_PROJ="$PWD"
export OVCAN_DATA="/actual/managed/path/to/archived/data"
python3 scripts/check_checkout.py --inputs
```

The example path is a placeholder, not a known cluster location. The optional check tests paths and sizes without changing files. Missing inputs in this check mean missing from that selected directory, not that files were never produced. Do not reconstruct historical folder names by guessing sample identity; use the canonical model/sample maps.

New cluster downloads belong in ignored `data/cluster_wes_retrieval/` or separate managed storage. Commit the curated inventory/report and small reviewed scripts/configs under `docs/cluster/recovery/`; keep sequencing and alignment payloads outside Git. `docs/cluster/evidence/manifest.csv` links copied original CNV commands to their local source paths and hashes.
