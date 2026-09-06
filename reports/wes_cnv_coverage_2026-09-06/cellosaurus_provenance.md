# Cellosaurus provenance recovery — 6 September 2026

**M06 is an analysis-team provenance item, not an author confirmation request.** The exact analysed snapshot is identifiable and internally reproducible. Its database-wide release number and original HTTP retrieval timestamp are not recoverable from the retained project records; these fields should be explicitly marked unavailable and the author-facing request closed.

## Evidence recovered

- **42 search-response JSON files, two direct accession-response JSON files, and one 42-name query list** remain under `output/external/cellosaurus/`. All 45 files independently match the byte counts, MD5 and SHA-256 values recorded by the September 5 audit in `output/external_reference_provenance.csv`.
- The 44 JSON bodies total **1,102,054 bytes**. They contain 37 record instances representing 30 distinct accessions. Their body keys are only `cell-line-list` and, when present, `publication-list`; no response timestamp or database-release object is present.
- Per-entry accession, creation date, entry version and last-update date survive. The distinct entry last-update dates are 10 April 2025 and 31 March 2026. Those are record metadata, not proof of a database release or download date.
- Thirty exact name matches and their reference STR marker counts agree with `output/cellosaurus_str_status.csv` and current `metadata/resource_models.csv`. Twelve named queries have empty cached results. This confirms which records were analysed without refreshing any external dataset.

The complete per-file checksums, documented API URLs, entry metadata and selected model accessions are in [cellosaurus_provenance.json](cellosaurus_provenance.json). A reproducible aggregate fingerprint for the 45-file sorted checksum manifest is `66c35d5eeed63fc4ec5dcde49e7f54718135283b9145eda909db0b63f5bb0f2e`; its precise construction is recorded in that JSON.

## What the dates establish

`reports/review_code.md:333` calls these the Cellosaurus JSON snapshot of 23 July, and `scripts/18_external_benchmarking.R:32` is dated 23 July 2026. This supports **23 July 2026 as a historically documented snapshot date**, without independently establishing the exact HTTP acquisition time.

The retrieval recipe `scripts/fetch_external_data.R` was added retrospectively: its header dated 24 July 2026 explicitly describes the original inputs as manually downloaded without a committed fetch script. Lines 189–207 implement the search/direct-accession endpoints and save response bodies with `download.file`; they do not save HTTP headers or acquisition times. Script 18 parses the saved files offline. Git history begins with the 5 September repository import and cannot supply the earlier acquisition timestamp. Filesystem timestamps were not used as historical download evidence.

The source endpoints are documented as `https://api.cellosaurus.org/search/cell-line?q=<name>&format=json` and `https://api.cellosaurus.org/cell-line/<accession>?format=json`. These are supported by the retained script and filenames, rather than a surviving original HTTP request log. Individual concrete URLs are recorded per file in the JSON audit.

A current API lookup would describe a new retrieval and cannot establish the missing original metadata. No current lookup, cache replacement, external-data refresh or scientific recomputation was performed.

## Recommended manuscript/checklist treatment

Suggested methods wording: “Cellosaurus information was taken from archived API responses documented as the 23 July 2026 snapshot. The retained records include accession and entry-version metadata, and the analysed response files are identified by checksums.” If a data-provenance field requires the database release or original HTTP timestamp, enter “unavailable in the retained snapshot” rather than guessing from entry dates.

Close M06 as: “Analysis-team provenance audit complete: cached Cellosaurus responses, source endpoints, entry versions and checksums are documented. Original database release and HTTP retrieval metadata were not retained.” This requires no additional author fact-finding.

M09 remains separate: Cellosaurus reference STR profiles do not authenticate the stocks profiled in this study, and current-stock STR/mycoplasma documentation still comes from the contributing laboratories. No accession, reference-status field or manuscript file was modified by this audit.
