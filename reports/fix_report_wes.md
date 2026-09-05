# Fix report — WES / authentication / external-benchmarking workstream

**Scope:** `scripts/00_setup.R`, `07_wes_mutations.R`, `08_wes_cnv.R`, `09_wes_hrd.R`, `10_authentication.R`, `11_mucinous_authenticity.R`, `15_patient_family_map.R`, `16_wes_signatures_msi.R`, `18_external_benchmarking.R`, `20_supplement_table.R`; new `scripts/22_wes_signature_refit.R` and `scripts/fetch_external_data.R`.
**Environment:** R 4.5.2, Bioconductor 3.21, aarch64-apple-darwin20. Every script below was executed to completion after editing.
**Not touched:** the manuscript, any `docs/` file, any `3*` figure script, `01`–`06`/`12`/`13`/`17`/`19` (RNA/proteomics agent), and `output/external/` (verified read-only, never rewritten).

---

## 1. Every fix

| # | Item | File:line | What changed | Before → after |
|---|---|---|---|---|
| C1 | Central mutation table written by no script | `07_wes_mutations.R:497-546` | Uncommented the `write_csv` for `output/wes_mutations_filtered.csv`; replaced the `nrow`-only drift guard with a **sha256 content hash** over the canonically-sorted identity columns, plus pinned row/line/median assertions | write commented out, guard = `nrow` only → file written every run; guard = `sha256 486051bc…` + `n_rows 6036` + `n_lines 22` + `median 206.5`. Cascade reproduces exactly: 557,392 raw → 15,692 PASS → 15,609 rare → **6,036** coding non-syn, median **206.5** (range 133–1,416) |
| C2 | Authentication used the chrX-contaminated FGA | `10_authentication.R:470-497` | `fga_0.2` → `fga_auto_0.2`; legacy column retained as `FGA_withX_legacy`, and a `stopifnot` asserts at most one line changes instability category | **OV90**: FGA 0.308 → **0.269**, `cnv_instability` high → **intermediate**, `genomics_consistent` consistent → **partial**. No other line changed category (verified) |
| C3 | Table S1 listed Tier-3 driver calls unflagged | `20_supplement_table.R:60-100, 140-192`; also `10_authentication.R:455-475` | Joined `wes_driver_tiers.csv`; emit `drivers_tier12`, `drivers_tier3`, `drivers_annotated` (`GENE(T3)`), `n_drivers_tier12/3`, and the hypermutator `driver_context` column that never previously reached Table S1. `stopifnot` asserts no BRCA1/2 call is Tier 1–2 | `key_drivers` only, no tier column → 4 tier-aware columns + context. **7 lines / 9 (line, gene) pairs** are Tier-3-only, matching the audit exactly. TOV81D: `key_drivers = "BRCA2"` → `drivers_tier12 = NA`, `drivers_tier3 = "BRCA2"`, `drivers_annotated = "BRCA2(T3)"` |
| 4 | ConsensusOV calls not computed | `18_external_benchmarking.R:405-505` | `consensusOV` **is** installed (1.30.0), so the calls are now genuinely computed with `get.consensus.subtypes()` on the 15 HGSC RNA lines. Emits all four RF class probabilities, the top-vs-second **margin**, the inherited label, agreement, and an input-set stability check (15 HGSC vs all 31 RNA lines) | string-parsed from `samples.csv` `notes` → computed; version + settings + seed recorded. See §2 for the number changes |
| 5 | TOV21G EPCAM call dropped from the narrative | `16_wes_signatures_msi.R:206-290` | Added `mmr_panel_hit`, `mmr_variants`, `mmr_any_truncating` and a per-line **`mmr_note`** carrying the tumour-only caveat; a dedicated console block prints every MMR-panel hit; figure panel-B subtitle states the caveat. Also fixed the codon extractor (`(?<=[A-Za-z])[0-9]+` → first integer), which had been silently returning `NA` for position-only HGVS forms | EPCAM buried in two boolean columns → unmissable, with the mechanism mismatch stated (see §2) |
| 6 | Signature evidence was cosine screening, not fitting | `16_wes_signatures_msi.R:66-88, 360-440`; new `22_wes_signature_refit.R` | 16: pinned `COSMIC_v3.2`/`GRCh38` + MutationalPatterns version into the outputs; report **which** signature carries each group maximum; report three **margins**; replaced the unjustified absolute 0.75 gate with a panel-relative margin rule (0.75 retained as a *reported* column `cos_mmr_d_ge_075`); wrote the full cosine profile to `wes_sbs_cosine.csv`. 22: `fit_to_signatures_strict` + `fit_to_signatures_bootstrapped` (200 replicates), bootstrap 95% intervals, per-signature selection frequency, reconstruction cosine, variants per spectrum, and a hypothesis-restricted reference set alongside the full 60-signature one | group max only → margins + full profile + a real refit. `sbs_call` unchanged for all 22 lines |
| 7 | chrX-exclusion statistic mis-scoped | `08_wes_cnv.R:245-265` | Compute and label **both** scopes explicitly in one table | "0.61→0.62" (actually all-23-line) → **all 23 lines 0.606 → 0.622**; **HGS only (n=18) 0.627 → 0.635**. Both audit values confirmed |
| 8 | Non-standard HGVS strings | `07_wes_mutations.R:110-194` | Rewrote the reconstruction: `fs` for frameshift, `*` for nonsense, `del`/`delins` with shared prefix/suffix trimming; added `hgvsp_reconstructed` (TRUE for every row — none came from the MAF) and `hgvsp_canonical` (safe to quote as an identifier) | `p.L639X` → **`p.L639fs`**; `p.R2845X` → `p.R2845fs`; `p.L123X` → `p.L123fs`; `p.CP1255-1256C` → `p.P1256del`; `p.N210X` → `p.N210fs`; `p.H214X` → `p.H214fs`; `p.-156-157P` → `p.156_157insP` (flagged non-canonical). `p.Q192*`, `p.R175H`, `p.I195T` unchanged. Non-canonical strings 159 → **103** of 6,036; 5,682 canonical, 251 not derivable |
| 9 | A literature value presented as a measurement | `11_mucinous_authenticity.R:56-71, 186-345` | Removed TOV2414's hard-coded verdict branch — all three lines now go through the same data-derived `expr_verdict()`. Literature claims explicitly attributed to Sauriol 2020 ("their IHC/sequencing, not measured here"). New `literature_vs_measured` column states the discordance. `ovarian_index` documented as an ad-hoc unweighted cross-assay composite (with `ovarian_index_definition` in the output). Thresholds named and swept (`auth_mucinous_sensitivity.csv`). Per-marker values against the **full** 31-line distribution (`auth_mucinous_marker_ranks.csv`). z-reference set stated in the output and printed | "SATB2−" read as in-house → measured **SATB2 z = +0.17, rank 19/31 (60th percentile)**, explicitly higher than VOA8762 (−0.59) and lower than VOA8771 (+0.80). All three verdicts survive **60/60** threshold combinations |
| 10a | `SMARCA2_loss` ignored somatic-confidence tier | `10_authentication.R:108-240` | Tier now travels with every SWI/SNF call (`*_wes` shows `[TierN]`); truncating status is admitted as evidence only at **Tier 1–2**; a new `swisnf_tier3_only_calls` column records what was excluded and why | OV2295's Tier-3 SMARCA2 truncation no longer makes it `swisnf_deficient`. TOV21G unaffected (its SMARCA2 call is in-frame; deficiency comes from Tier-1 ARID1A) |
| 10b | "26 expression-consistent" rests largely on absence | `10_authentication.R:560-580, 688-702` | New `expression_basis` column + a printed breakdown | 26 = **8** positive lineage program (5 CC + 3 MC) + **16** absence-of-competing-program (15 HGS + 1 MMMT) + **2** SCCOHT circular. Exactly the audit's split |
| 10c | `has_drv()` regex over a joined string | `10_authentication.R:574-582` | Replaced with a lookup over the structured `drv_struct` table | identical results, no string parsing |
| 10d | No SWI/SNF figure input exists | `10_authentication.R:290-345` | New long-format `auth_swisnf_long.csv` (248 rows): per line × gene × assay — raw value, z, rank, rank denominator, z-reference set, WES call, tier, whether admitted. Protein ranks and raw log2 values added to `auth_swisnf_panel.csv` | **BIN67 quantified**: SMARCA4 mRNA rank **23/31** (log2 TPM 6.58, z +0.62) vs protein rank **2/31** (log2 13.32, z −1.96), margin to the next-lowest protein value **0.105 log2** |
| 11 | Two hypermutator definitions | `00_setup.R:139-160`; `07:402-410`; `16:129-140` | One `hypermutator_stats()` helper in `00_setup.R` (robust_z > 5 **and** > 3× median); both scripts call it | `07` used ">3× median", `16` used the conjunction → one rule. TOV21G remains the only hypermutator |
| 12a | DepMap collapse/HVG asymmetry undocumented | `18:110-130`, `external_selfmatch_margin.csv` | Documented in the header comment and carried into the output as `collapse_note` / `hvg_note` | — |
| 12b | Self-match margin not emitted | `18:228-300` | New `external_selfmatch_margin.csv` (per line: self ρ, rank of 67, best non-self, margin, non-self median/SD, z vs the non-self distribution, reciprocal margin, all-gene sensitivity) and `external_depmap_spearman_all.csv` (all 31 × 67 = 2,077 pairs). `stopifnot` asserts rank 1, reciprocal-best, positive margin | **margins: OV90 0.434, TOV21G 0.231, TOV112D 0.204, COV434 0.139, BIN67 0.036**; all rank 1 of 67 and reciprocal-best. BIN67's margin is thin — the referee's question 6 now has a number |
| 12c | Cellosaurus fell back to the search top hit | `18:535-556, 600-612` | Exact name match only; no match now returns "not found". `stopifnot` asserts exact-match-only and that every RRID cross-check passes | 30/42 exact matches, 12 not found, 0 fuzzy — unchanged, but a future rename fails loudly instead of attaching a wrong accession |
| 12d | DepMap driver cross-check console-only | `18:340-405` | New `external_depmap_driver_crosscheck.csv` (55 line × gene rows with an explicit `agreement` verdict) and `external_depmap_burden.csv` | corroborated: TOV112D TP53 R175H + SMARCA4, OV90 SMAD4, all 5 TOV21G CC drivers. Not corroborated: TOV112D KRAS A59T (Tier 2), OV90 BRAF `p.N486_P490del` (Tier 2), TOV21G SMARCA2 (Tier 3). DepMap-only: TOV21G SMARCA4 damaging |
| 13 | `00_setup.R` recorded no versions | `00_setup.R:57-118, 163-170` | New `record_versions()`, called on every `source()`: writes `output/package_versions.csv` (R, Bioconductor, platform, OS, seed, timestamp + 40 analysis packages) and `output/session_info.txt` | claim was false → true |
| 14 | `output/external/` reproducible by no script | new `scripts/fetch_external_data.R` | Documented fetch script with a per-file manifest (bytes, md5, rows, cols). **Default mode is verify-only** and writes nothing; `--download` fetches only missing files; `--force` required to overwrite. DepMap URLs resolved through the figshare REST API by file name (article 27993248) rather than hard-coded file IDs; Cellosaurus via its REST API; the four derived subsets regenerated from the raw files | verified: **10/10 files match on size, md5 and row/column counts**; 42 Cellosaurus search JSON + 2 record JSON; 684 MB total. Nothing overwritten |
| 15 | Numeric order ≠ dependency order | `00_setup.R:113-137`; `07:249`; `08:100`; `16:77`; headers of `15`, `07`, `08` | `ensure_family_map()` regenerates `metadata/line_family_map.csv` from `samples.csv` if absent, by sourcing `15`. Chose the guard over renumbering because the manuscript, `reports/` and `docs/` all cite "script 15" — renaming would have broken those references (flagged in §5) | clean-checkout run of `07` aborted → runs |
| 15b | Patient-level FGA and rare-subtype labelling | `08:220-243` | Documented that a patient's value is the **mean** of its lines (clear cell = mean of 0.671 and 0.071); added a `label_note` naming the constituent lines for every n ≤ 2 subtype | "EC 0.226 (n=1)" now reads "labelled subtype, n≤2 lines (TOV112D) — see 10 for reclassification" |
| 16 | Arm gain/loss threshold unrecorded | `08:62-71, 248-340` | Named `ARM_LOG2_THRESH` (0.20) and `ARM_MAJORITY` (0.50), written into `wes_cnv_arm_freq_patient.csv` along with `centring` and `pct_per_patient`; new `wes_cnv_arm_freq_sensitivity.csv` sweeps both parameters | **All six manuscript frequencies reproduce exactly**: 3q gain 82%, 20q gain 91%, 17p loss 82%, 8q gain 73%, 13q loss 64%, 19q gain 55% (n = 11 patients, 9.1% each). Sensitivity: at \|log2\| 0.10–0.40 the six values move 82→64, 100→73, 82→64, 73→55, 73→45, 55→36 |
| 17 | Recoverable method parameters unwritten | `09_wes_hrd.R:70-133` | New `output/wes_pipeline_parameters.csv` (15 rows) parsed from the archived `commands.txt` | CNVkit **v0.9.10** (`cnvkit-0.9.10--pyhdfd78af_0.img`), GRCh38 `Homo_sapiens_assembly38.fasta`, `--drop-low-coverage --diagram --scatter`, `-t intervals_sorted.bed -a intervals.antitarget.bed`, pooled normals `normal_samples/all_bams/*.bam`, 23 samples. **12 of 15 recoverable; 3 require the PI** |

---

## 2. Manuscript-facing numbers that changed

| Claim | Before | After | Where |
|---|---|---|---|
| **"20 genomics-consistent"** | 20 consistent, 2 partial | **19 consistent, 3 partial** | C2. OV90 moves to partial because its autosome FGA (0.269) does not clear the stated FGA > 0.30 rule; with chrX it read 0.308 |
| OV90 `cnv_instability` | high | **intermediate** | C2 |
| OV90 `FGA` in `auth_perline_table.csv` | 0.308 (chrX-inclusive) | **0.269** (autosome) — now agrees with Table S1's `fga_autosome`, which was already correct | C2 |
| SMARCA4 variant identifier | `p.L639X` | **`p.L639fs`** | item 8. Also `p.R2845X`→`p.R2845fs` (TOV81D BRCA2), `p.L123X`→`p.L123fs` (CDK12), `p.N210X`→`p.N210fs`, `p.H214X`→`p.H214fs` |
| Table S1 drivers for TOV81D | `key_drivers = "BRCA2"` (its only driver) | `drivers_tier12 = NA`, `drivers_tier3 = "BRCA2"` | C3 |
| Table S1 drivers for TOV3133D | `"BRCA2, CDK12, TP53"` | `drivers_tier12 = "CDK12, TP53"`, `drivers_tier3 = "BRCA2"` | C3 |
| **ConsensusOV class counts** | DIF 6 / MES 4 / IMR 3 / PRO 2 (inherited) | **DIF 5 / MES 5 / IMR 3 / PRO 2** (computed, consensusOV 1.30.0) | item 4. 14 of 15 lines agree with the inherited labels; **OV1369-R2** flips DIF → MES |
| **"7 of 15 mesenchymal or immunoreactive"** | 7 | **8** | item 4, follows from the above |
| ConsensusOV reliability evidence | bare labels | **five lines have a top-vs-second margin below 0.10**: OV2085 **0.022** (IMR 0.318 vs PRO 0.296), TOV3041G 0.044, OV90 0.072, TOV3291G 0.078, OV1369-R2 0.088; median margin 0.160. **2 of 15 calls (OV2085, TOV3041G) change when the input set changes** from the 15 HGSC lines to all 31 RNA lines | item 4. This is much stronger support for the paper's own argument than the labels were |
| chrX-exclusion effect on FGA | "high-FGA HGSC genomes 0.61→0.62" | **all 23 lines 0.606→0.622; HGS only (n=18) 0.627→0.635** — the quoted pair was the all-line median | item 7 |
| "SBS6 cosine 0.88, plus SBS44/SBS15/**SBS20**" | text says SBS20; Fig. 5C says SBS15 | **the figure is right**: TOV21G's top three are SBS6 **0.877**, SBS44 **0.835**, SBS15 **0.809**; SBS20 is 4th at **0.565**. Text should read SBS6/SBS44/SBS15 | item 6 |
| Signature evidence for TOV21G | cosine 0.877 vs 0.51–0.70 in the other 21 lines (1.25× the best of the rest) | **refit MMR-d relative exposure 0.733 vs 0.000–0.292 (2.5× the next line)**; SBS6 alone 1,001 of 2,340 fitted variants (rel 0.428), bootstrap 95% CI **271–1,062**, selected in **99.5% of 200 bootstraps**; reconstruction cosine **0.977**, the best in the panel (others 0.892–0.960). MMR-d-minus-clock cosine margin **+0.225** for TOV21G, **−0.270 to −0.040** for all 21 others | item 6. The refit strengthens the claim rather than weakening it |
| "No MMR-enzyme coding mutation is present" (TOV21G) | asserted | **two MMR-panel hits exist and are now unmissable.** TOV21G **EPCAM `p.2_3fs` (Frame_Shift_Ins)** — the only truncating MMR-panel hit; TOV2881EP **MSH3 `p.Y334F`** (missense, ClinVar uncertain significance, VAF 0.649). *New nuance:* EPCAM is on the panel because **3′-end** deletions silence MSH2 in cis, but this call is a **5′ frameshift at codon 2–3**, which truncates EPCAM itself and does **not** invoke that mechanism. It is novel in dbSNP, has no VAF support (indel), and sits in the panel's one hypermutator | item 5 |
| DepMap self-match specificity | "rank 1 of 67, reciprocal-best" | rank 1 of 67 and reciprocal-best confirmed for all five, **with margins**: OV90 **0.434**, TOV21G **0.231**, TOV112D **0.204**, COV434 **0.139**, **BIN67 0.036** (reciprocal margin 0.056). BIN67's identity claim is much thinner than the other four and should be stated as such | item 12b |
| TOV2414 "KRT7+/PAX8+/MUC5AC+/**SATB2−**" | reads as four in-house measurements | three are (KRT7 z **+2.58**, PAX8 **+0.28**, MUC5AC **+4.25**); **SATB2 z = +0.17, rank 19/31, 60th percentile** — mid-panel, and *higher* than VOA8762 (−0.59). SATB2-negative is Sauriol's IHC and must be attributed | item 9 |
| Mucinous verdicts | thresholds unstated | robust: **60/60** threshold combinations give TOV2414 ovarian-compatible and VOA8762 GI-leaning; VOA8771 is GI in 60/60 but splits between two GI flavours (48 colorectal-leaning / 12 GI-leaning) | item 9 |
| "26 expression-consistent" | one number | **8 positive lineage program + 16 absence-of-competing + 2 circular** | item 10b |
| SWI/SNF-deficient lines | 7 (OV2085, **OV2295**, TOV21G, VOA4841, TOV112D, BIN67, COV434) | **6** — OV2295 drops out because its only SWI/SNF evidence was a **Tier-3** SMARCA2 truncation (germline not excludable in tumour-only WES). Its `flags` entry ("secondary SWI/SNF subunit variant") disappears, and `supplement_per_line.csv$swisnf_deficient` changes for that row. The two reclassification lines (COV434, BIN67) and TOV112D are unaffected | item 10a |

Numbers **verified unchanged**: the filtering cascade (557,392 / 15,692 / 15,609 / 6,036; median 206.5, range 133–1,416); 50 driver calls, Tier 1/2/3 = 27/11/12; all six arm frequencies; subtype FGA ordering (HGS 0.622 > CC 0.371 > MC 0.321 > EC 0.226 > LGS 0.021); TOV21G 1,416 coding candidates, 2,417 exome-wide SNVs; the 5-line DepMap overlap; 30/42 Cellosaurus STR records; `sbs_call` for all 22 lines; BIN67 SMARCA4 mRNA rank 23/31 and protein rank 2/31.

**Row-for-row check on the regenerated central table.** The audit could not confirm whether the archived `wes_mutations_filtered.csv` came from the current filter logic. I diffed the regenerated file against the pre-revision snapshot on the variant key across all 21 shared columns: **identical in every column except `HGVSp_Short`** (786 rows, all intentional per item 8), plus the two new provenance columns. The archived table was indeed the product of this logic.

---

## 3. New output files

| File | Rows | Contents |
|---|---|---|
| `output/wes_cnv_arm_freq_sensitivity.csv` | 60 | The six quoted arm frequencies across \|log2\| 0.10–0.40 × arm-majority 0.30–0.70, with `is_headline` |
| `output/wes_sbs_cosine.csv` | 1,320 | Full cosine profile: 22 lines × 60 COSMIC v3.2 signatures, with group, release and `n_snv_used` |
| `output/wes_signature_refit_exposures.csv` | 1,250 | Per line × signature × reference set: absolute and relative exposure, bootstrap mean/median/95% CI, selection frequency, method parameters, versions |
| `output/wes_signature_refit_summary.csv` | 22 | Per line: variants used, signatures selected, reconstruction cosine (both reference sets), group exposures, top MMR-d signature with its interval, side by side with the cosine screen |
| `output/wes_signature_refit_bootstrap.csv` | 1,250 | Cache of the 40-minute bootstrap, keyed on the fit parameters (delete or set `OVCAN_REFIT_FORCE=1` to recompute) |
| `output/wes_pipeline_parameters.csv` | 15 | Machine-readable WES pipeline provenance with `recoverable`, `evidence` and `action` per parameter |
| `output/auth_swisnf_long.csv` | 248 | **SWI/SNF figure input**: line × gene × assay value, z, rank, rank denominator, z-reference set, WES call, tier, admitted-as-evidence |
| `output/auth_mucinous_marker_ranks.csv` | 496 | Every mucinous marker for every line in both assays: value, z, rank, percentile, panel min/median/max, z-reference set |
| `output/auth_mucinous_sensitivity.csv` | 180 | Mucinous verdict across 60 threshold combinations × 3 lines |
| `output/external_selfmatch_margin.csv` | 5 | **DepMap identity figure input**: self ρ, rank of 67, best non-self, margin, z vs the non-self distribution, reciprocal margin, all-gene sensitivity |
| `output/external_depmap_spearman_all.csv` | 2,077 | Every our-line × DepMap-ovarian-line correlation (HVG and all-gene), with `is_selfpair` |
| `output/external_depmap_driver_crosscheck.csv` | 55 | Line × gene: our call and tier vs DepMap hotspot/damaging counts, with an explicit `agreement` verdict |
| `output/external_depmap_burden.csv` | 5 | DepMap total damaging/hotspot counts vs our PASS coding burden |
| `output/package_versions.csv` | 46 | R, Bioconductor, platform, OS, seed, timestamp + 40 analysis package versions |
| `output/session_info.txt` | — | `sessionInfo()` + seed + calling script |
| `figs/f_wes_signature_refit.pdf` | — | Refit exposures, TOV21G bootstrap intervals, reconstruction cosines (diagnostic; available for a figure agent) |

Columns added to existing outputs: `wes_mutations_filtered.csv` (+`hgvsp_canonical`, `hgvsp_reconstructed`); `wes_cnv_arm_freq_patient.csv` (+`log2_threshold`, `arm_majority_frac`, `centring`, `pct_per_patient`); `wes_msi_mmr.csv` (+`mmr_panel_hit`, `mmr_variants`, `mmr_any_truncating`, `mmr_note`, cosine margins, `cos_*_sig`, `cos_mmr_d_ge_075`, `cosmic_source`, `cosmic_genome`, `mutationalpatterns_version`); `auth_swisnf_panel.csv` (+protein ranks, raw log2 values, `*_trunc_tier12`, `swisnf_tier3_only_calls`); `auth_perline_table.csv` (+`expression_basis`, `key_drivers_tier12/tier3/annotated`, `FGA_metric`, `FGA_withX_legacy`, `swisnf_tier3_only_calls`); `auth_mucinous.csv` (+`SATB2_rna_rank`, `ovarian_index_definition`, `expression_verdict`, `literature_vs_measured`, z-reference and threshold columns); `consensusov_calls.csv` (+four class probabilities, margin, second call, inherited label, agreement, input-set stability, version, settings, provenance); `supplement_per_line.csv` (+the C3 tier columns, `driver_context`, `expression_basis`, `fga_withX_legacy`, ConsensusOV probability/margin/provenance, `rna_sequenced_fragments_M`).

---

## 4. Incomplete or deliberately not done

1. **`session_info.txt` duplication.** My `record_versions()` in `00_setup.R` writes one `output/session_info.txt`; the RNA agent independently added a local `write_session_info()` to nine of their scripts, producing `session_info_<script>.txt`. Both work; they should be consolidated to one convention before Methods describes it. `package_versions.csv` is script-independent and is the better single citation.
2. **Script `15` was not renumbered.** I implemented the `ensure_family_map()` guard instead, because the manuscript, `reports/05_scientific_data_descriptor_draft.md` and `reports/02_manuscript_outline.md` all cite it as "script 15". If the PI prefers numeric order to equal dependency order, renaming to `06b_patient_family_map.R` also requires updating those references — outside my scope.
3. **No exome-to-genome trinucleotide renormalisation** in the signature refit, because the capture kit is unrecoverable. Stated in the script header and in `wes_signature_refit_summary.csv$exome_renormalisation`. Relative exposures are therefore approximate; the ranking within a line and the contrast between identically-processed lines is what carries.
4. **`renv.lock` not created.** `package_versions.csv` is a substitute, not an equivalent. Creating a lockfile touches the whole project environment and should be done deliberately, last, by the PI.
5. **The mucinous panel was not extended** with external anchors (TCGA/CPTAC mucinous ovarian vs colorectal) as referee M13 suggests. That needs new external data and a new download workstream; I documented and stress-tested the existing panel instead.

---

## 5. What I need from files I do not own

- **Manuscript / `docs/`** (all of these follow mechanically from §2):
  - "20 genomics-consistent" → **19**; and OV90 must be described as *partial*, not consistent.
  - Every quotation of `p.L639X` → `p.L639fs` (and `p.R2845X`, `p.L123X`, `p.N210X`, `p.H214X` likewise). Table S1's `drivers_tier12`/`drivers_tier3` columns should be cited where "defensible somatic BRCA1/2 is zero" is claimed.
  - ConsensusOV: **DIF 5 / MES 5 / IMR 3 / PRO 2**, "8 of 15" mesenchymal-or-immunoreactive, described as **computed with consensusOV 1.30.0 in script 18** (not "generated … script 10"). Report the margins and the 2-of-15 input-set instability.
  - "SBS44/SBS15/**SBS20**" → **SBS6/SBS44/SBS15**; describe the cosine step as *screening* and cite `22_wes_signature_refit.R` for the fit.
  - The EPCAM sentence must be replaced: a truncating EPCAM call **is** present, with the 5′-vs-3′ mechanism caveat and the tumour-only germline caveat.
  - chrX statistic: state which scope (all 23 lines 0.606→0.622, or HGS-only 0.627→0.635).
  - TOV2414: attribute SATB2-negative to Sauriol 2020 and report the measured z (+0.17) alongside.
  - "26 expression-consistent" needs the 8 / 16 / 2 split.
  - Arm frequencies need the call rule stated: \|log2\| > 0.20 over > 50% of arm length, per-sample autosome median centring, n = 11 patients (9.1% each).
  - BIN67's DepMap self-match margin is 0.036 — the weakest of the five; "reciprocal-best" should not be stated without it.
- **Figure scripts (`3*`, not mine).** Three tables were written specifically for them: `auth_swisnf_long.csv` (the SWI/SNF panel the paper claims but no figure contains), `external_selfmatch_margin.csv` + `external_depmap_spearman_all.csv` (the rank-1-of-67 margin plot a 5×5 submatrix cannot show), and `consensusov_calls.csv` (per-class probabilities for Fig. S6). `wes_signature_refit_exposures.csv` and `figs/f_wes_signature_refit.pdf` are available if the signature panel is upgraded. Note `31_fig4_genomics.R:150-152` hard-codes 25,914/493 and `32_fig5_rare.R` hard-codes 1,416 and 6.9× — all still correct after this work, but they should read the CSVs.
- **RNA agent.** `20_supplement_table.R` now consumes `assigned_gene_counts` and `n_processed_fragments` from `rna_qc_metrics.csv` (with a fallback to the old `lib_size`). If those names change again, `20` needs a matching edit.
- **PI.** The three non-recoverable parameters in §7.

---

## 6. Missing packages

**None blocked any analysis.** Everything needed was present, including the two that mattered most:

- `consensusOV` **1.30.0** — installed, so item 4 took the preferred route (compute) rather than the fallback (rename the parse). A documented fallback path remains in the script for environments without it.
- `MutationalPatterns` **3.18.0**, `BSgenome.Hsapiens.UCSC.hg38` **1.4.5**, `digest` **0.6.39** — all present.
- `variancePartition` **1.38.1** is installed and fails at run time (used only by `17`, which I do not own). **Code Availability's claim that the package is "unavailable in this R build" is wrong** — it is installed; it errors when called, and `17` correctly substitutes an `lme4` REML equivalent. The accurate statement is "installed but non-functional in this environment; an equivalent per-gene lme4 REML decomposition is used".
- `scarHRD` and `copynumber` remain absent, consistent with `09`'s feasibility verdict — but that verdict is driven by the **data** (no allele-specific copy number), not by the tooling.
- Two masking hazards worth recording because both caused real failures during this work: `consensusOV` exports `margin`, which masks `ggplot2::margin` (so it must be called namespaced, never attached), and `S4Vectors` (pulled in by `MutationalPatterns`) exports `rename`/`first`/`second`, which mask the dplyr verbs. `matrixStats::count` masks `dplyr::count` in `11`.

---

## 7. Method parameters: recovered vs still missing

**Recovered and now machine-readable** in `output/wes_pipeline_parameters.csv` (12 of 15 rows):

| Parameter | Value |
|---|---|
| CNV caller | `cnvkit.py batch`, **v0.9.10** (container `cnvkit-0.9.10--pyhdfd78af_0.img`) |
| Reference genome | GRCh38 GATK bundle, `Homo_sapiens_assembly38.fasta` |
| CNVkit flags | `--drop-low-coverage --diagram --scatter` |
| Target / antitarget | `intervals_sorted.bed` / `intervals.antitarget.bed` |
| Panel of normals | `normal_samples/all_bams/*.bam` — 5 public healthy exomes, PRJNA339046 |
| Samples processed | 23 |
| Mutect2 matched normal | none (tumour-only; `n_depth`/`n_ref_count`/`n_alt_count` empty in every MAF) |
| Genome build actually used | GRCh38/hg38 (the MAF `NCBI_Build` field reads `GRCh37` — a spurious vcf2maf default; flag as a deposition item) |
| COSMIC release | **COSMIC_v3.2 / GRCh38**, the newest COSMIC bundled with MutationalPatterns 3.18.0 (v3.3+ unavailable here) |
| DepMap release | Public 24Q4, figshare article 27993248 |
| consensusOV | 1.30.0, `concordant.tumors.only = TRUE`, `remove.using.cutoff = FALSE`, seed 1234 |
| Arm-call rule | \|log2c_auto\| > 0.20 over > 50% of arm length, per-sample autosome probe-weighted median centring |

**Still missing — requires the PI (3 rows, flagged as `recoverable = FALSE` with an `action`):**

1. **Exome capture kit identity.** `intervals_sorted.bed` is referenced by every CNVkit call but is not itself archived, and nothing in the archive records the kit. This also blocks the pooled-normal validity question: the 5 public normals come from an unrelated skin-cancer study and **capture-kit concordance with the tumours is unconfirmed**, which is a validity concern for the CNV layer, not merely a hygiene one.
2. **Proteomics search parameters — genuinely unavailable.** `judy_archive/data/proteomics/Readme.md` is a **1-byte empty file**. Instrument, MS2-vs-MS3 acquisition, search engine and version, sequence database, enzyme, modifications, tolerances, FDR method and thresholds, interference filter, and peptide-to-protein rollup are all unrecoverable from the archive and must come from the Morin laboratory. Referee M2 is correct that the proteomics record is not reproducible as described, and no amount of code fixes it.
3. **Sarek / GATK versions** for the SNV arm are not in the archived `commands.txt` (only the CNVkit invocations are). Recoverable in principle from a Nextflow log or `pipeline_info/` directory if one survives on HPC scratch.
