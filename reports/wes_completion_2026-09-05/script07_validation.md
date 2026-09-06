# Script 07 guard validation

Ran `/usr/local/bin/Rscript scripts/07_wes_mutations.R` **once**, with exit status **0**. The log is `script07_validation.log`; before/after SHA-256 inventories and the comparison result are recorded beside this report.

The new assertion rejecting nonnumeric population-AF annotations in MAF PASS records passed. The existing retained-content guard also passed:

- 6,194 retained candidates from 23 models; median 206 candidates per model.
- Pinned R identity-column digest: `4febdd74482a09ef5818be7f6c7ffe2cd53f09e5e750bfcb4d86bacb1967e6b4`.
- Filtering cascade remained 582,474 annotated rows → 16,081 MAF PASS rows → 15,995 passing the additional population-AF filter → 6,194 coding nonsynonymous candidates.

| Protected files | Files present before run | Changed | Removed |
| --- | ---: | ---: | ---: |
| `output/wes*.csv` | 37 | 0 | 0 |
| All files under `metadata/` | 3 | 0 | 0 |
| Canonical manuscript figure PDFs/PNGs | 29 | 0 | 0 |

One new CSV, `output/wes_recovered_provenance_cnv_reference_filters.csv`, appeared during the observation window. Script 07 does not write this path; the provenance agent confirmed it was their concurrent Script 23 output, with no Script 07 dependency. All 37 pre-existing WES CSVs are byte-identical.

The script refreshed its exploratory oncoplot exports and normal runtime records. Canonical manuscript figures and the manuscript were not modified by this validation. The log contains 23 warnings whose individual text was not printed by the batch session, followed by successful VAF joining, and a ComplexHeatmap vectorisation notice. No assertion failed and no protected scientific-table checksum changed. The script was not rerun to expand these warnings.
