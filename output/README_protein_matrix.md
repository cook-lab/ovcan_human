# The protein layer: one matrix, one row annotation, one accounting table

Several different protein feature counts appear in this project, and all of them are correct —
they count different things. Rather than choose one and discard the rest, the deposit is
structured so that every count is derivable from two files.

## Start here

| File | Shape | What it is |
|---|---|---|
| `prot_abundance_matrix.csv` | 8,427 × 31 | **The canonical matrix.** Supplied log2 protein abundance; confirm the upstream pooled-standard normalization/scaling formula, one row per protein feature, one column per model. This is the file to reanalyse. |
| `prot_qc.csv` | 8,427 × 15 | **The row annotation.** One row per matrix row, in the same order, carrying peptide support, q-value, replicate CV, and the four presence flags below. |
| `prot_feature_accounting.csv` | 14 × 3 | **The reconciliation.** Every feature count quoted anywhere in this project, with its definition. |

`prot_matrix.rds` is the same matrix as `prot_abundance_matrix.csv` in R binary form, loaded by
`read_prot_matrix()` in `scripts/00_setup.R`.

## Choosing a feature set

Any subset a reanalysis needs is a filter on `prot_qc.csv`, not a different file:

| Want | Filter | n |
|---|---|---|
| Everything quantified in at least one model | `present_n_lines >= 1` | 8,357 |
| **Recommended default for cross-model comparison** | `pass_presence50` (`present_n_lines >= 16`) | 7,733 |
| Complete cases only, for PCA and variance decomposition | `complete_case` | 6,855 |
| One row per gene symbol | `symbol_representative` | 8,396 |

Two flags exist to prevent specific mistakes. `zero_plex` marks the 70 features present in the
search output but quantified in no model; they are retained so that the matrix row count matches
the search output, and they must be dropped before any per-protein statistic. `present_n_plex`
records how many of the five TMT plexes a feature was quantified in — missingness is structural
per plex, so a feature absent from a whole plex is absent from every model in that plex, and
set-conditional presence is not biology.

## Row naming

Rows are named by gene symbol. Where two accessions map to the same symbol, the representative is selected by most observed models, then most quantified peptides,
then lowest q-value, and finally accession in alphabetical order. It keeps the bare symbol; the others are named
`SYMBOL|UNIPROT` and carry `symbol_representative = FALSE` (31 rows). Three search-output rows
have a UniProt accession but no gene symbol and cannot enter a symbol-keyed analysis; they are
excluded from the matrix, which is why 8,430 search rows become 8,427 analysis rows.
