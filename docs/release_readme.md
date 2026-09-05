# OvCAN processed-data resource

The resource contains 42 ovarian cancer cell models from 34 patients. Start with `metadata/models.csv`: every row belongs to the resource, and no inclusion filter is required. The internal project retains its larger source inventory solely for provenance and bridge-sample accounting.

Use `file_catalog.csv` to locate records, `DATA_DICTIONARY.md` to interpret them, and `field_inventory.csv` to inspect individual column names and missing-value counts. `SHA256SUMS` records the content checksum of every packaged file. `validation.json` records the model counts and cross-file consistency checks performed when this package was built.

## Reuse

- Join records by `cell_line`. Columns carrying model names use the same identifiers. Model order can differ between RNA and protein matrices; align names before analysis.
- Related models share `patient_id`. The selected-model flag provides one model per patient for the complete resource. For a particular assay intersection, select a measured model within each patient; do not assume the global selection maximises every intersection.
- RNA TPM and estimated-count matrices have different gene filters. The reference map includes primary-assembly and alternative-locus transcripts present in the quantification reference. Use the reference annotations when restricting genomic context.
- Protein abundance is supplied on a log2 scale after normalisation using a pooled internal standard. The upstream scaling formula requires laboratory confirmation. Empty/`NA` measurements remain missing, including 70 identified features with no quantified abundance. Missingness is structured by TMT plex and does not indicate absent protein. `prot_qc.csv` rows align exactly to the protein matrix.
- Sequence alterations are tumour-only coding candidates. Evidence tiers support prioritisation; they do not establish somatic origin, loss of heterozygosity, pathogenicity, or drug response.
- Copy-number segments describe relative total copy ratio. Autosomal summary values are preferred. Arm calls intersect segments with reference arm boundaries and use the fraction of full annotated non-centromeric arm length meeting the declared threshold. Patient aggregation is recorded in each summary.
- Historical histotype groups remain visible. In particular, TOV112D has been reassigned to dedifferentiated ovarian carcinoma. Histotype and contributing centre are partly confounded.
- STR and mycoplasma documentation must be confirmed with contributing laboratories. A pending documentation field does not mean testing was not performed.

## Build and scope

From the project root, `bash scripts/run_all.sh --no-fetch` rebuilds processed analyses and figures using the existing source inputs and pinned reference data. Restore the R package environment recorded in `renv.lock` before a full scientific rerun. Then rebuild the local viewer and package:

```sh
python3 scripts/pin_external_reference_provenance.py  # pin/verify the existing external cache
python3 app/build_payload.py  # requires pandas and numpy
python3 app/build_single.py
python3 scripts/build_release.py  # standard library only
```

The `validation/` directory contains the principal consistency and sensitivity tables. The `exploratory/` directory contains target-opportunity signature sensitivity results with their limitations recorded in each table; these are not diagnostic classifications. A full analysis requires the scripts and source inputs in the accompanying project; this processed-data package alone is not a substitute for the raw data and code archive.

External-reference hashes identify the exact cached inputs. DepMap primary files match the official 24Q4 version-1 record. The Cellosaurus date is a historically documented snapshot date; original retrieval headers and the database release were not preserved. A future refresh must record its own retrieval metadata and must not be labelled as the original snapshot.

This directory packages processed data. Raw reads, raw proteomic files, accessions, licensing terms, complete experimental metadata, and repository reviewer access require separate finalisation; their presence is not implied by this package. Do not publish until those items and collaborator confirmations are complete. The analysis audit and its author action list are under `reports/audit_2026-09-05/` in the project.
