# Integration report — clean-room reproducibility test

**Date:** 2026-07-24
**Environment:** R 4.5.2, `Rscript`, aarch64-apple-darwin20, Bioconductor 3.21, `OVCAN_PROJ` default.
**Scope:** all of `scripts/` except the `3*` figure scripts. `docs/`, `reports/review_*.md`,
`reports/_snapshot_*/` and `output/external/` untouched.
**Companion:** `reports/integration_diff.md` (per-file comparison), `reports/integration_logs/`
(per-script stdout+stderr, `timings_pass1.tsv` / `timings_pass2.tsv`, the comparison and
verification scripts, and `pass1_08_wes_cnv_familymap_guard.log`).

---

## 1. Verdict

**Yes — a clean run from raw inputs reproduces the full output set.**

From an `output/` containing nothing but the protected `external/` directory, and with
`metadata/line_family_map.csv` deleted, 23 script invocations ran to completion with
**zero non-zero exits** and regenerated **all 122 baseline files**. Nothing is missing;
nothing is written by no script.

Of those 122 files: **109 identical**, **12** are `session_info*.txt` provenance records
that carry a run timestamp and therefore cannot be identical by construction, and
**1** — the bootstrap cache — legitimately grew because the baseline copy was written by a
superseded code path (§5). No file differs in a numeric value, a column name, or a row
count for any other reason. All **25** headline numbers reproduce exactly (§6).

Three qualifications, stated plainly:

1. **As certified, the run required network access.** `output/tx2gene_ensembl_rel105.csv`
   is a pinned reference *input*, but it used to live only in `output/` — the directory
   documented as "safe to regenerate". Emptying `output/` therefore made
   `01_rna_load_qc.R` fall back to a **live biomaRt query** against
   `dec2021.archive.ensembl.org`, in **both** passes. So the certified claim is precisely:
   *reproducible with network access to the Ensembl 105 archive*, **not** offline.
   The map came back **byte-identical both times** (md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`,
   266,615 transcripts, matching the hand-cached copy and the baseline snapshot), so
   nothing drifted — but the two claims are different and only the weaker one was earned.
   **Now fixed:** the map has been relocated to `data/reference/`, `01` reads it from
   there with an **asserted md5**, and the offline path is verified (§4 R6). A clean run
   from this commit forward needs **no network access**.
2. **The first pass was contaminated and had to be discarded.** The RNA repair agent was
   executing pipeline scripts into the same `output/` while my run was in flight — it has
   since confirmed this directly (its runs: `03, 04, 12, 13, 19, 17, 21`; last write
   15:29). The certified result above is a second, exclusive pass. Details in §2.
3. **`output/_preclean/` has been deleted**, as instructed, because every file was
   reproduced. `reports/_snapshot_postfix/` remains as the archive.

---

## 2. What went wrong operationally (and why the verdict is still sound)

This is not a code defect, but it determined how the test had to be run, so it belongs
in the record.

**Pass 1** (12 scripts, all exit 0) was contaminated. While `03_rna_de_signatures.R` was
running, files appeared in `output/` from scripts I had not launched —
`session_info_{12,13,19}.txt` at 15:20, `session_info_17.txt` at 15:24,
`session_info_21.txt` at 15:29 — and `metadata/line_family_map.csv` was recreated at
15:17:43 by another process hitting the `ensure_family_map()` guard. The RNA repair agent
was still executing its scripts. The outputs it produced were derived from my regenerated
upstream files and nothing looked wrong, but a pass in which 5 of 22 steps were executed
by an unlogged third party cannot be certified. I abandoned it as evidence, archived its
products and logs, and started again.

**Pass 2** is the certified run. Its driver (`reports/integration_logs/run_pass2.sh`)
refuses to start a script while any other `--file=scripts/` R process is alive, so pass 2
never shared `output/` with anything.

That guard fired four times on figure scripts from other agents (`30_fig1_overview.R`,
`33_supp_genomics.R` twice, `36_fig6_adc.R`) and waited them out — a total delay of
about 10 s. **Those were false positives, and I record the correction here because I
initially reported them as contamination.** Those agents run with `OVCAN_PROJ` pointed at
an isolated sandbox, so they read and write their own `output/`, `metadata/` and `docs/`,
never the real ones. The evidence: `00_setup.R` stamps the target root into
`output/session_info.txt`, and the sandbox copy names `project /private/tmp/.../figsandbox`
while the real copy names the real root; and all 28 files in the real
`docs/manuscript/figures/` are byte-identical to `reports/_snapshot_pre_revision/`, which
could not be true if a figure script had run against the real root. **The figure work is
not void and does not need rebuilding on my account.**

My stated mechanism for that claim — that figure scripts write `package_versions.csv` and
`session_info.txt` into the clean room — is only true when `OVCAN_PROJ` is *unset*, because
`00_setup.R:17-20` derives `OUT` from `PROJ`. It was set. The correct exclusivity test keys
on the **target root**, not the script name: treat a process as foreign only if its
`OVCAN_PROJ` is absent or equals the real project root (`ps -Eww -p <pid>`), or check
`lsof +D output/` directly. My name-based guard is over-broad — harmless here, but it would
stall indefinitely against a long sandboxed job.

Separately, and unrelated to the sandbox: none of `00`–`22` sources
`scripts/00b_figure_theme.R` or uses its palettes (`17` and `18` define their own local
`save_fig()`), so the concurrent revisions to that file could not have affected this run.

One self-inflicted incident worth recording: my first driver's exclusivity check matched
the string `"Rscript scripts/"`, which never matches, because `Rscript` execs R with
`--file=scripts/<name>`. It let `18_external_benchmarking.R` start alongside `22`. I
killed `18`, deleted the 5 partial CSVs it had written, fixed the pattern, and re-ran `18`
from scratch in pass 2. `output/external/` was verified intact throughout
(55 files, 652 MB, never opened for writing).

Both `metadata/samples.csv` and `output/external/` are unmodified. `reports/_snapshot_pre_revision/`
and `reports/_snapshot_postfix/` are unmodified.

---

## 3. Failure log

**No script failed in pass 2.** The defects below were found and fixed *before* or
*during* the run; each is a genuine cross-half or clean-checkout bug, not a symptom patch.

| # | Defect | `file:line` | Root cause | Fix |
|---|---|---|---|---|
| F1 | `ensure_family_map()` guard incomplete — `17`, `20`, `21` read `metadata/line_family_map.csv` bare | `17_variance_confounders.R:85`, `21_rna_sensitivity.R:56`, `20_supplement_table.R:46-49` | The guard was added to `07`/`08`/`16` only. `line_family_map.csv` is *generated* by `15`, so any entry point at `17`/`20`/`21` on a clean checkout reads a file no script has produced — **the same class of bug as the original `wes_mutations_filtered.csv` defect** | `17` and `21` now call `ensure_family_map()`. `20` deliberately sources nothing (pure join, must stay cheap), so it sources `15` directly under the same condition |
| F2 | Bootstrap cache validity checked keys but not coverage; an incomplete cache would be silently reused and would **drop rows** from the exposure table | `22_wes_signature_refit.R:202-213`, `:262-265` | Validity tested `n_boots`, `cache_max_delta`, `cache_n_ref_full`, `cache_n_ref_restricted` and the cell-line set. The short 1,250-row cache on disk passed all of them. The join then uses `replace_na(boot_selected_frac, 0)` in its final `filter()`, so an absent cache row is indistinguishable from a never-selected signature | Validity now also requires the complete `n_lines x n_signatures` grid per reference set; and a `stopifnot` asserts full coverage and no duplicates before the join. Both the reject-and-recompute and the reuse paths were then tested end to end |
| F3 | `10`/`11` carried dead guards written for the old 8,430-row matrix | `10_authentication.R:102`, `11_mucinous_authenticity.R:83` | See §4 (R1) | Both read through one shared loader |
| F4 | Stale row counts in comments | `06:21`, `06:34`, `19:23` | Documented the matrix as 8,430 proteins | Corrected to 8,427 |
| F5 | A pinned reference **input** lived in the regenerable output directory, so emptying `output/` silently triggered a **live network query** | `01_rna_load_qc.R:49-115` | `tx2gene_ensembl_rel105.csv` (14 MB, 266,615 transcripts) was only ever in `output/`. `01` loads it if present and otherwise re-queries `dec2021.archive.ensembl.org` via biomaRt. 8 scripts (`01, 03, 04, 12, 13, 17, 19, 21`) consume it, and it determines the 39,568-gene matrix, so a silent drift would move every downstream number with no warning | Relocated to `data/reference/`; `01` loads it from there, **asserts md5 `119dfe0ab4f856e6b81efc6ce78f4ba7` and 266,615 rows**, refreshes the `output/` working copy by byte copy so the other 7 consumers are unchanged, and keeps the biomaRt fallback but now warns loudly when it is used. Offline path verified (§4 R6) |
| F6 | My own first attempt at the F5 fix was wrong | `01_rna_load_qc.R:104-108` | I checksummed a `write_csv()` re-serialisation of the parsed tibble rather than the file. A `read_csv` -> `write_csv` round-trip is not byte-preserving here (empty `external_gene_name` values return as `NA` and re-serialise differently), so the assertion failed against the very file it was meant to verify | Assert on the authoritative file on disk and propagate to `output/` with `file.copy()`, never by re-serialising |

Interrupted-run artefact, for completeness: `18`'s 5 partially written CSVs (from my own
tooling error, §2) were deleted and regenerated, not patched.

---

## 4. Cross-half reconciliation

The two repair workstreams edited disjoint files, so four contracts had never been
exercised together. All four now hold.

### R1 — the protein matrix shape change never reached `10`/`11`

`10` and `11` consumed `prot_abundance_matrix.csv` when it still had 8,430 rows, three of
them named `NA` / `NA.1` / `NA.2`, with duplicate symbols as `SYMBOL.N`. Both carried a
private guard written for that shape:

```r
pa %>% filter(!is.na(protein), protein != "", !duplicated(protein))
```

After the RNA agent's change (8,427 rows, unique names, `SYMBOL|UNIPROT`) that guard
removes nothing. It never crashed — it just meant the two halves no longer agreed on what
the file was, and the dead code hid the change.

**Resolved:** both now call `read_prot_matrix()` (§R2). All five consumers report the same
matrix in their logs: `8427 proteins x 31 lines (8.8% NA) | 70 zero-plex rows retained`.

**Scope of the rename:** only 31 rows are `SYMBOL|UNIPROT`, from 13 duplicated symbols
(CDKN2A, CHTF8, CUX1, HLA-A, HLA-B, HLA-C, HLA-DRB1, MIEF1, POLR1D, RAB34, RABGAP1L,
TMPO, ZFP64). None is a SWI/SNF panel gene, a mucinous marker, or an ADC target, so no
marker-level number in `10`/`11`/`13` moves because of it. This is now **asserted in code**,
not assumed.

### R2 — one zero-plex approach across `10`, `11`, `12`, `13`, `19`

The RNA agent could not put a `zero_plex` flag inside `prot_abundance_matrix.csv` because
`10` and `11` read it as a bare numeric matrix.

**Decision: keep the deposited matrix numeric-only and make the contract enforced rather
than merely documented.** A flag column would coerce the whole matrix to character for
every consumer that does `as.matrix(x[, -1])` — which is exactly what all five scripts and
any external reuser do. The flag already has a natural home in `prot_qc.csv`, the
feature-level table keyed on `row`.

**Implementation — one function, `read_prot_matrix()` at `00_setup.R:136-181`**, which:

* asserts no NA, empty or duplicated row names (i.e. that `05` did its job);
* **derives** the zero-plex set from the matrix itself (`rowSums(!is.na(m)) == 0`);
* **cross-checks it against `prot_qc.csv$zero_plex`** with a `stopifnot`, so the two files
  cannot drift apart in either direction;
* returns it as `attr(m, "zero_plex")`, with an optional `drop_zero_plex`.

Per-consumer handling, chosen so that **no denominator moves**:

| Script | Treatment of the 70 zero-plex rows |
|---|---|
| `10` | retained; asserts no SWI/SNF panel protein is zero-plex (`10:134-135`) |
| `11` | excluded before `t(scale(t()))`, which would return `NaN` for an all-NA row; asserts no mucinous marker is zero-plex (`11:83`, `11:116-117`) |
| `12` | retained; they contribute 0 to `n_prot_per_gene` and are removed by the `n>=10` rule, preserving the 8,212 shared-symbol denominator |
| `13` | retained; asserts no ADC target is zero-plex — "absent protein" and "never quantified" are different claims for an ADC candidate (`13:95-96`) |
| `19` | removed **by name** before abundance-decile binning, and the two sets asserted equal (`19:409-421`) |

**Confirmed: no zero-plex row has ever contributed a value to any output.** The decile
denominator was already 8,357 = 8,427 − 70 — but it got there incidentally, because
`rowMeans(all-NA, na.rm = TRUE)` is `NaN` and a downstream `filter(is.finite(...))`
dropped it. Nothing would have noticed if a zero-plex row had acquired a value. It is now
deliberate and self-checking. `12`'s numbers are unchanged: 8,212 shared symbols, 7,894
per-gene correlations, median Spearman 0.397.

### R3 — `20`'s `rna_qc_metrics.csv` contract

Holds. The regenerated file carries `assigned_gene_counts` (median 56,497,657) and
`n_processed_fragments` (median 64.3 M, range 45.1–97.0 M); `lib_size` is gone.
`20:118`'s `else lib_size` fallback sits inside a lazily evaluated `if`, so it is never
evaluated — back-compat for an older file, not a stale read. `supplement_per_line.csv`
has **42 rows x 55 columns** and both RNA depth columns are populated; coverage reads
RNA 31 / protein 31 / WES-MAF 22 / tri-omic 13, and the `genomics_consistent` tally is
**19 consistent / 3 partial / 1 discordant / 19 not assessed** — i.e. the WES C2 fix
(OV90 moving to *partial*) is intact end to end.

`20`'s own family-map fallback was tested separately, because pass 2 ran `15` first and so
never exercised it: with `metadata/line_family_map.csv` deleted, `20` run standalone logs
`metadata/line_family_map.csv absent -> running 15_patient_family_map.R`, rebuilds the map
**byte-identically**, and writes a **byte-identical** `supplement_per_line.csv`
(`reports/integration_logs/20_guardtest.log`).

### R4 — no script reads a column that no longer exists

* `lib_size` — only in comments and the unevaluated fallback above.
* `fga_0.2` outside `08` — four sites (`10:494`, `10:497`, `10:481` comment, `20:142`).
  All are **valid**: `08` deliberately retains `fga_0.2` as the chrX-inclusive legacy
  column in `wes_cnv_fga.csv` (confirmed present in the regenerated file at column 5), and
  every consumer names its result `*_withX_legacy`.
* No `SYMBOL.N` row names survive: the regenerated matrix has 8,427 rows, 0 names ending
  in `.N`, 31 containing a pipe, 0 NA.

### R5 — family-map guard, verified with logged evidence

With `metadata/line_family_map.csv` removed, `08_wes_cnv.R`'s log
(`reports/integration_logs/pass1_08_wes_cnv_familymap_guard.log`, preserved for this
purpose — pass 2 ran `15` first, so its own `08` log cannot show the guard) reads:

```
metadata/line_family_map.csv absent -> running 15_patient_family_map.R
Wrote /Users/dpcook/Analysis/ovcan_human/metadata/line_family_map.csv
Generated lines: 42  ->  independent patients: 34  (collapsed 8 lines)
```

The regenerated map was **byte-identical** to the pre-clean copy (verified before
`output/_preclean/` was removed, and re-confirmed independently by the `20` test in R3).
Pass 2 additionally ran `15` standalone (exit 0, 1 s), so the direct path, the
`ensure_family_map()` path (`08`) and `20`'s own fallback are all tested.

### R6 — the tx2gene reference map: what actually happened, and the fix

**What happened, stated exactly:** in *both* passes `01_rna_load_qc.R` **re-queried the map
live** from `dec2021.archive.ensembl.org`. It was not restored by me. The pass-2 log
(`reports/integration_logs/01_rna_load_qc.log:6-7`) reads
`Querying Ensembl release 105 (pinned archive dec2021) for tx2gene ...` /
`Cached tx2gene -> ... (266615 transcripts)`.

**The result is a genuinely strong reproducibility datum.** Two independent live queries,
~30 minutes apart, produced a file byte-identical to the hand-cached copy and to the
baseline snapshot — md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`, 13,884,511 bytes, 266,615
transcripts, all four copies. The pinned Ensembl 105 archive reproduces exactly.

**But it made the wrong claim true.** "Reproducible" and "reproducible offline" are
different, and only the former was earned, for a file that no one intended to be
regenerable. Fixed as F5: the map now lives in `data/reference/`, is checksum-asserted on
load, and is copied (not re-serialised) into `output/` for the other seven consumers.

**Verified:** with `output/tx2gene_ensembl_rel105.csv` deleted, `01` logs
`Loaded pinned tx2gene: .../data/reference/...` and
`tx2gene verified: md5 119dfe0ab4f856e6b81efc6ce78f4ba7 | source: data/reference (pinned input)`
— **no network query** — then reports 39,568 genes x 31 samples and the 22,544-gene filter,
exit 0. All eight files `01` writes are **byte-identical** to the certified pass-2 copies, so
the relocation moved no number and did not perturb the certified run
(`reports/integration_logs/01_offline_test.log`).

---

## 5. Non-determinism

**None found.** No unseeded RNG was detected and no seed fix was required.

The single changed file, `wes_signature_refit_bootstrap.csv` (1,250 -> 1,804 rows), was
checked specifically for this and is not non-determinism:

* Pass 1 and pass 2 — two independent from-scratch 200-replicate bootstraps — produced
  **byte-identical** caches.
* A third recompute, forced by planting the short baseline cache, produced a
  byte-identical cache again, plus byte-identical `wes_signature_refit_exposures.csv` and
  `wes_signature_refit_summary.csv`.
* 1,804 = 22 x 60 (full) + 22 x 22 (restricted) — the complete grid. The baseline's
  916 + 334 held only the ever-selected subset, i.e. a superseded `boot_summary()`.
* Both derived result tables are identical to the baseline, so no reported statistic
  depends on the extra rows.

That investigation is what surfaced defect **F2**, which is the more important finding:
the cache was *reusable while incomplete*.

**A second, independent determinism result** comes from F5: the transcript-to-gene map was
fetched over the network twice, ~30 minutes apart, and both fetches were byte-identical to
each other and to the hand-cached copy (md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`). The pinned
Ensembl 105 archive is stable. That is now protected by an assertion rather than by luck.

Two benign sources of run-to-run churn, neither a number:

1. **Timestamps.** `output/session_info.txt` (rewritten on every `source()` of `00_setup.R`)
   and 11 `session_info_<script>.txt` files record when a script ran.
   `package_versions.csv` compares **identical** once its `recorded_utc` row is excluded —
   which is why it is the better single citation for the environment.
2. **Float printing across a CSV round-trip.** On the cache-reuse path,
   `wes_signature_refit_exposures.csv` is numerically identical (`all.equal` TRUE, same row
   order) but not byte-identical: a final significant digit is dropped on a few values,
   e.g. `96.72750068668061` -> `96.7275006866806`. The delivered `output/` holds the
   recompute-path files.

`sessionInfo()`'s "loaded via a namespace (and not attached)" list also varies with which
namespaces R lazy-loaded during a run (visible in `session_info_01_rna_load_qc.txt`). That
is session introspection, not analysis output.

---

## 6. Headline number verification

Re-derived from the regenerated `output/` only, by
`reports/integration_logs/verify_headlines.R` (machine-readable result in
`headline_check.tsv`). **25 checks, 25 PASS, 0 moved.**

| Check | Expected | Observed | Verdict | Source |
|---|---|---|---|---|
| raw Mutect2 records | 557,392 | 557,392 | PASS | `07` cascade |
| FILTER==PASS | 15,692 | 15,692 | PASS | `07` cascade |
| rare (pop AF > 0.001) | 15,609 | 15,609 | PASS | `07` cascade |
| coding non-synonymous | 6,036 | 6,036 | PASS | `07` cascade |
| `wes_mutations_filtered.csv` rows | 6,036 | 6,036 | PASS | file |
| WES lines | 22 | 22 | PASS | file |
| median variants / line | 206.5 | 206.5 | PASS | file |
| TOV21G coding candidates | 1,416 | 1,416 | PASS | file |
| variance-decomposition genes | 22,542 | 22,542 | PASS | `rna_variancepartition.csv` |
| protein matrix rows | 8,427 | 8,427 | PASS | `prot_feature_accounting.csv` |
| complete case | 6,855 | 6,855 | PASS | same |
| absent from >= 1 plex | 1,572 | 1,572 | PASS | same |
| **6,855 + 1,572 == 8,427** | 8,427 | 8,427 | PASS | same |
| `prot_abundance_matrix.csv` rows | 8,427 | 8,427 | PASS | file |
| zero-plex proteins | 70 | 70 | PASS | `prot_zero_plex_proteins.csv` |
| `supplement_per_line.csv` rows | 42 | 42 | PASS | file |
| `genomics_consistent == "consistent"` | 19 | 19 | PASS | file |
| RNA depth: assigned counts populated | TRUE | TRUE | PASS | file |
| RNA depth: sequenced fragments populated | TRUE | TRUE | PASS | file |
| arm 3q gain | 82 % | 82 % | PASS | `wes_cnv_arm_freq_patient.csv` |
| arm 20q gain | 91 % | 91 % | PASS | same |
| arm 17p loss | 82 % | 82 % | PASS | same |
| arm 8q gain | 73 % | 73 % | PASS | same |
| arm 13q loss | 64 % | 64 % | PASS | same |
| arm 19q gain | 55 % | 55 % | PASS | same |

Also reproduced exactly, spot-checked from the logs:

* `07`'s content guard: `sha256(identity cols) = 486051bc5277540b891bebd164246d97fec8314501e6e22d2f7e56ec69fe8d41`.
* chrX-exclusion scope: all 23 lines 0.606 -> 0.622; HGS only (n=18) 0.627 -> 0.635.
* BIN67 SMARCA4 mRNA rank 23/31 (log2 TPM 6.58, z +0.62) vs protein rank 2/31
  (log2 13.32, z −1.96), margin 0.105.
* TOV21G MMR-d relative exposure 0.733; SBS6 selected in 99.5 % of 200 bootstraps;
  reconstruction cosine 0.977 (others 0.892–0.960).
* Mucinous verdicts stable in 60/60 threshold combinations; VOA8771 splits
  48 colorectal-leaning / 12 GI-leaning.
* `12`: 8,212 shared symbols, 7,894 per-gene correlations, median Spearman 0.397.
* `05`: 8,430 search rows, 3 dropped (A6NIZ1, A6NNZ2, Q6ZSR9), 8,396 distinct symbols,
  31 non-representative rows, `>=50%` presence 7,733.

**No headline number moved.**

---

## 7. Runtime

Pass 2, sequential, one script at a time. Total **1,775 s ≈ 29.6 min** of compute.
Full table in `reports/integration_logs/timings_pass2.tsv`.

| Script | s | | Script | s |
|---|---|---|---|---|
| `00_setup.R` | 0 | | `22_wes_signature_refit.R` | **571** |
| `15_patient_family_map.R` | 1 | | `18_external_benchmarking.R` | 67 |
| `01_rna_load_qc.R` | 42 | | `12_rna_protein_concordance.R` | 7 |
| `05_proteomics_load_qc.R` | 4 | | `13_adc_atlas.R` | 2 |
| `02_rna_separation.R` | 4 | | `19_proteomics_dynamic_range.R` | 10 |
| `06_proteomics_separation.R` | 4 | | `14_hgs_heterogeneity.R` | 4 |
| `03_rna_de_signatures.R` | **520** | | `10_authentication.R` | 2 |
| `04_rna_markers_genesets.R` | 7 | | `11_mucinous_authenticity.R` | 2 |
| `07_wes_mutations.R` | 13 | | `17_variance_confounders.R` | **243** |
| `08_wes_cnv.R` | 4 | | `21_rna_sensitivity.R` | **253** |
| `09_wes_hrd.R` | 4 | | `20_supplement_table.R` | 0 |
| `16_wes_signatures_msi.R` | 11 | | | |

Three notes for anyone budgeting a rerun:

* The pipeline is **much cheaper than the pre-run estimate**. `22` takes 9.5 min from
  scratch with no cache (not ~40 min); `17` takes 4 min (not ~40 min). Only `03` (fgsea at
  `nPermSimple = 50000` x 6 contrasts), `22`, `21` and `17` exceed a minute.
* `22`'s runtime was 570 s in pass 1 and 571 s in pass 2 — both from an empty cache.
  With a valid cache present it is a few seconds.
* Whole-pipeline wall time is ~30 min on this machine, so the full chain is a practical
  pre-submission check rather than an overnight job.

---

## 8. Residual risks

1. **Concurrency.** `output/` is a single shared namespace with no lock. Two agents (or two
   terminals) running scripts against the *same* `OVCAN_PROJ` silently interleave writes —
   this happened during pass 1 and cost a full re-run. `reports/integration_logs/run_pass2.sh`
   implements an advisory guard, but nothing in `scripts/` enforces it. A real fix would be
   a lockfile acquired in `00_setup.R`. Note that the guard as written keys on the **script
   name**, which is over-broad: it blocks on any pipeline script, including ones pointed at a
   different `OVCAN_PROJ`. A correct check keys on the **target root** — treat a process as
   foreign only when its `OVCAN_PROJ` (`ps -Eww -p <pid>`) is unset or equals the real
   project root, or use `lsof +D output/`. `OVCAN_PROJ` is what makes a sandboxed run safe,
   so it must be set on every invocation that is not meant to touch the real `output/`.
2. **Figure scripts are outside this test.** The `3*` scripts were not run as part of the
   certified chain and are therefore **not covered by this certification**. That remains a
   real caveat: a future clean run should include them.

   *[Corrected 2026-07-24 after this section was first drafted. The original text stated
   that the figure work needed rebuilding and that required edits "remain unmade". Both
   claims were wrong — they were inferred from `fix_report_rna.md` §5, which was written
   before the figure revision finished. Verified in code and on disk:*

   * *The concurrent figure processes were running with `OVCAN_PROJ` pointed at a separate
     sandbox root holding a frozen copy of `output/`. They never read this clean room, and
     all 28 files in `docs/manuscript/figures/` were byte-identical to the pre-revision
     snapshot throughout the run — so nothing they produced was derived from half-written
     inputs.*
   * *The hard-coded labels are fixed. `35_fig3_biology.R:222-225` now reads `n_features`
     from `rna_variancepartition.csv` / `prot_variancepartition.csv` under a `stopifnot`
     (22,542 / 6,855); `22,544` and `6,856` survive only in an explanatory comment.
     `35:166-167` uses `n_pergene_reported` = 7,894, retaining 8,212 as a separately
     labelled quantity. `37_supp_rnaprot.R:449` reads
     `grp_lvl <- c("HGS","CC","EC","MC","MMMT","SCCOHT")`.*
   * *All eight figure scripts were re-rendered against the certified outputs after this
     run completed, in an order placing `34`/`35`/`37` after `17`/`19` (four
     `reports/assets/` filenames are written by both a pipeline script and a figure
     script, so the later writer wins). All 8 exited 0; 28 of 28 figure files updated.*

   *The re-render is a separate, uncertified step. The accurate residual risk is that the
   `3*` scripts have no clean-room coverage — not that their content is stale.]*
3. **`output/session_info.txt` records only the last script sourced**, by design, so its
   content depends on run order. `package_versions.csv` is the order-independent record and
   should be the one Methods cites. The two conventions (`00_setup.R`'s single
   `session_info.txt` plus 11 per-script `session_info_<script>.txt`) are still
   unconsolidated, as both fix reports noted.
4. **`.DS_Store`.** A Finder artefact appeared in `output/` during this work. It is not a
   pipeline product and no script writes it; it is simply not part of the 122.
5. **`diptest` is still not installed**, so the Hartigan dip test for FOLR1 bimodality
   remains substituted by `mclust` BIC + the bimodality coefficient (documented in
   `adc_folr1_bimodality.csv`). Unchanged by this work.
6. **`variancePartition` 1.38.1 is installed but non-functional** in this R (it calls
   `lme4::findbars`, which moved to `reformulas`). `17` uses its documented `lme4` REML
   equivalent. A future R/package upgrade could change which path runs — the
   `model` column in `rna_variancepartition.csv` records which one did.
7. **`output/external/` is reproducible only by `scripts/fetch_external_data.R` in
   `--download` mode**, which was not exercised here (it defaults to verify-only, and I
   deliberately never wrote to that directory). The 652 MB of DepMap/Cellosaurus data
   remains an input a clean checkout cannot regenerate without network access.
8. **`data/reference/tx2gene_ensembl_rel105.csv` is now a protected input** (F5). If it is
   deleted, `01` still runs — but it falls back to a live Ensembl query, which needs network
   access. The md5 assertion means a *drifted* map fails loudly, so the residual risk is
   availability, not silent corruption. It should be deposited alongside the data.
9. **Write coverage is fully accounted for.** A separate audit confirmed that the 23
   baseline outputs which look unreferenced to a literal grep are written by constructed
   filenames — `rna_de_<SUBTYPE>.csv` and `rna_signatures_<SUBTYPE>.csv` (12 files) via
   `sprintf()` in `03`'s subtype loop, and `session_info_<script>.txt` (11 files) via
   `write_session_info()`. Every other baseline file has an identifiable writer. With F5
   resolved, **no file in `output/` is written by no script**, which closes the question
   this test was commissioned to answer.
10. **No lockfile.** `renv.lock` still does not exist. `package_versions.csv` records
    versions but does not pin them.

---

## Execution order

Canonical order for Methods and Code Availability. Numeric order is **not** dependency
order; this is. Each script sources `scripts/00_setup.R`, which sets paths, the seed
(`SEED = 1234`), and writes `output/package_versions.csv`.

```
00_setup.R                 # sourced by all others; run once to confirm it writes package_versions.csv
15_patient_family_map.R    # writes metadata/line_family_map.csv (consumed by 07, 08, 16, 17, 20, 21)
01_rna_load_qc.R
05_proteomics_load_qc.R
02_rna_separation.R
06_proteomics_separation.R
03_rna_de_signatures.R
04_rna_markers_genesets.R
07_wes_mutations.R
08_wes_cnv.R
09_wes_hrd.R
16_wes_signatures_msi.R
22_wes_signature_refit.R
18_external_benchmarking.R
12_rna_protein_concordance.R
13_adc_atlas.R
19_proteomics_dynamic_range.R
14_hgs_heterogeneity.R
10_authentication.R
11_mucinous_authenticity.R
17_variance_confounders.R
21_rna_sensitivity.R
20_supplement_table.R      # pure join; safe to run last
```

Then the figure scripts (`30`–`37`), which read only from `output/` and are outside this test.

Notes on the order:

* `15` is listed first because `07`, `08`, `16`, `17`, `20` and `21` all consume
  `metadata/line_family_map.csv`. It is not strictly required: `ensure_family_map()` in
  `00_setup.R` regenerates the map from `metadata/samples.csv` if it is absent, and this
  was verified by deleting the map and running `08`.
* `18` must precede `10` and `20`, which consume `consensusov_calls.csv` and
  `cellosaurus_str_status.csv`.
* `05` must precede `12`, `13` and `19`, which read `prot_abundance_matrix.csv` through
  `read_prot_matrix()`.
* **`19` must follow both `12` and `13`**, not just `05`: it reads `integ_rnaprot_cor.csv`
  (from `12`) and `adc_expression.csv` (from `13`). Do not move it earlier to "optimise"
  the order — its position between `13` and `14` is load-bearing.
* `22` must follow `16` (it reads `wes_sbs_context.csv`); `18`, `10` and `11` must follow
  `07` (`wes_mutations_filtered.csv`).
* `20` must run last: it is a pure join over already-materialised results.
* Scripts must be run **one at a time** against a given `OVCAN_PROJ`. They share `output/`,
  and concurrent runs interleave writes (§8.1).

**Required inputs a clean checkout must have on disk** (none is regenerated by the chain):
`metadata/samples.csv`, `judy_archive/data/`, `output/external/` (652 MB, hand-downloaded;
`scripts/fetch_external_data.R` verifies it and can fetch it with `--download`), and
`data/reference/tx2gene_ensembl_rel105.csv` (md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`,
asserted at `01`). `metadata/line_family_map.csv` is *generated* by `15` and is regenerated
automatically if absent.

Reproduced in full on 2026-07-24 in **1,775 s** of sequential compute from an empty
`output/` (excluding the protected `output/external/`). That run reached the Ensembl 105
archive over the network for the transcript-to-gene map; with the map now pinned in
`data/reference/` (F5), the chain runs **offline**.
