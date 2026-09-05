# Revision delta — WES drivers & CNV (workstream: `wes-drivers`)

Addresses scientific review §3.1 (pseudoreplication), §3.3 (tumor-only somatic
confidence), §3.4 (CNV validity + chrX), and the minor reconciliations (§5).
All numbers below are recomputed and materialized; scripts rerun end-to-end.

**Scripts edited:** `scripts/07_wes_mutations.R`, `scripts/08_wes_cnv.R`
**Read-only preserved:** `output/wes_mutations_filtered.csv` is **not** rewritten
(07 recomputes it deterministically and `stopifnot`-checks the on-disk row count;
timestamp unchanged — safe for concurrent readers).

---

## 1. Patient-level re-tabulation (the headline fix)

Frequencies are now computed at the **independent-patient** level using
`metadata/line_family_map.csv` (a driver counts for a patient if ≥1 of that
patient's lines carries it). MAF set = **22 lines → 16 patients**; HGSC-with-MAF
= **17 lines → 11 patients**. `n_patients` is reported alongside `n_lines`
everywhere.

| Claim | Per **line** (old) | Per **patient** (new) | Verdict |
|---|---|---|---|
| **TP53 in HGSC** | 17/17 (100%) | **11/11 (100%)** | **Survives** — positive control intact |
| TP53 all-subtype | 18/22 (82%) | 12/16 (75%) | denominators now explicit |
| **CDK12 in HGSC** | 6/17 (**35.3%**) | **3/11 (27.3%)**; Tier1–2 only **2/11 (18.2%)** | **Collapses** — was one family's frameshift ×4 |
| BRCA2 (HGS) | 1/17 (5.9%) | 1/11 all-calls; **0/11 Tier1–2** | indefensible (see §2) |

- **CDK12 mechanism:** the 3133 family (OV3133-R/-R2, TOV3133D/G) all share the
  same `p.L123X` frameshift — **one patient event counted four times**. The only
  other CDK12 patients are OV3291 (splice, VAF 0.985) and TOV2929D (`p.S325I`,
  VAF 0.09 → Tier3). Defensible independent CDK12 signal ≈ **2 patients**, in
  line with the ~3% TCGA-HGSC rate — *not* a striking enrichment.
- **New file:** `output/wes_driver_freq_patient.csv` (gene × subtype: n_lines,
  n_lines_mut, pct_lines, n_patients, n_patients_mut, pct_patients,
  n_patients_mut_tier12, pct_patients_tier12).
- **Augmented:** `output/wes_driver_freq_by_subtype.csv` now carries the five
  patient-level columns (original columns preserved).

## 2. Somatic-confidence tiering (replaces `germline_like_vaf`)

Per-**call** tiers (`output/wes_driver_tiers.csv`, 50 calls;
cols: cell_line, patient_id, family, subtype, gene, protein_change,
variant_classification, vaf, germline_like_vaf, clin_sig, tier, rationale):

- **Tier1 canonical (27 calls):** known hotspot (TP53 DBD/truncating, KRAS
  12/13/61, PIK3CA 1047, CTNNB1 degron) or clear truncating LOF in a canonical
  TSG (PTEN, ARID1A, SMARCA4, NF1, RB1).
- **Tier2 plausible (11):** damaging, non-canonical position (KRAS A59, BRAF
  in-frame β3-αC, CDK12/ARID1B truncating, CDKN2A damaging+LOH).
- **Tier3 cannot-exclude-germline (12) — EXCLUDED from headline frequencies:**
  both BRCA2 calls, CDK12 in-frame/low-VAF, non-degron CTNNB1, uncharacterised
  ARID1A/NF1 missense, all three SMARCA2 calls.

**`germline_like_vaf` is not a somatic filter** and we no longer treat it as one:
21/50 driver calls are flagged TRUE, including bona-fide somatic TP53 hotspots
(R175H, Q192*, R249W) whose VAF→1 is **LOH**, not germline. In ~100%-pure lines a
VAF threshold cannot separate germline from somatic-with-LOH. Keep the column as
descriptive metadata only.

**Specific dropped/annotated calls (both → Tier3, excluded):**
- **TOV81D BRCA2 `p.R2845X`** — incoherent: TOV81D is the **quietest** genome in
  the panel (autosome FGA **0.021**, HR-proficient); a biallelic pathogenic BRCA2
  loss is incompatible with a near-diploid genome. Tumor-only WES cannot call a
  rare germline BRCA2 variant somatic. **Do not list BRCA2 as a "key driver" for
  TOV81D.**
- **TOV3133D BRCA2 `p.A22E`** — N-terminal missense, VAF **0.038** (at the noise
  floor), no established pathogenicity.
- **Add caveat (essential):** "tumor-only WES with population-AF filtering cannot
  classify a rare germline BRCA1/2 variant as somatic; BRCA1/2 statements carry
  this caveat." Defensible somatic BRCA1/2 rate in this panel = **0**.

**Headline driver claims restricted to robust events:** TP53 (universal, LOH-
supported), KRAS hotspots (CC: TOV21G G13C, TOV3392D G12S; MC: TOV2414 G12A —
all G12/G13, Tier1), and **SMARCA4 `p.L639X` truncation in TOV112D** (Tier1;
supports the SWI/SNF-null dedifferentiated call). CTNNB1 S37A (TOV112D, degron,
Tier1), PIK3CA H1047Y + PTEN L265X (TOV21G CC, Tier1) are also defensible.

**TOV21G hypermutator context** (cross-note from `wes-signatures`; flagged in the
new `context` column of `wes_driver_tiers.csv`). TOV21G is a candidate MSI-high /
MMR-deficient clear-cell line — **1,416 coding candidates = 6.9× panel median,
3.4× the next line** (SBS6/15/44). In a hypermutator, passengers land in driver
genes by chance, so its 8-gene list is inflated. Its **canonical CC drivers are
still real** (ARID1A truncating ×2, PIK3CA H1047Y, KRAS G13C, PTEN L265X — all
Tier1); its **secondary calls are already down-tiered** (ARID1B Tier2; CTNNB1 A5V,
SMARCA2 Q250- Tier3), so the tier framework itself guards against the inflation.
**Caution:** the two-line CC "driver frequency" (n=2 patients) is dominated by
TOV21G — read CC per-gene rates as descriptive, not recurrence, and interpret
TOV21G's SWI/SNF calls in the hypermutator context (full treatment in the
`wes-signatures` figure `f_wes_hypermutation`).

## 3. CNV: autosome-restricted FGA, chrX, and validity

- **chrX was INCLUDED in FGA (correction beyond the review).** The review inferred
  from `total_assessed_mb` ≈ 2.9 Gb that FGA was already autosome-restricted — it
  was **not**. Autosomes alone assess ≈**2746 Mb**; the ≈2902 Mb in the old file
  = autosomes **+ chrX** (≈2875 Mb autosome total < 2902, so chrX must be in).
  **FGA is now autosome-restricted** (`fga_auto_0.2` etc.).
- **chrX diagnostic** (`wes_cnv_fga.csv`: `chrX_median_log2c`, `chrX_frac_altered`):
  per-line chrX median log2c ranges **−1.23 → +0.83** (median 0.11), with a median
  **96%** of chrX called "altered" and the **sign flipping between lines** — a
  systematic sex-composition shift vs the sex-mixed pooled normal, not tumor CN.
- **Impact of removing chrX:** modest for high-FGA HGSC lines (median FGA 0.61→0.62)
  but large for quiet genomes — **TOV81D 0.073 → 0.021** (chrX inflated it ~3.5×).
  Autosome FGA by subtype (patient-level median, n patients): **HGS 0.62 (n=11) >
  CC 0.37 (n=2) > MC 0.32 (n=1) > EC 0.23 (n=1) > LGS 0.02 (n=1)**. Ordering holds
  but rests on **n=1** for MC/EC/LGS and the two CC lines differ wildly
  (TOV3392D 0.69 vs TOV21G 0.12) — state this.
- **chrX dropped from the figure** (chrY already was); it no longer contradicts the
  report's own caveat.
- **Median-centering caveat:** per-sample probe-weighted median-centering miscalls
  the neutral baseline when >50% of the genome is altered. High-FGA lines flagged
  `high_fga_flag=TRUE`: **TOV3121D, TOV2929D** (autosome FGA >0.70).
- **Unmatched public normals = VALIDITY concern, not hygiene.** The pooled reference
  is 5 unmatched public exomes (PRJNA339046, a skin-cancer study). If the capture
  kit differs, kit-specific coverage masquerades as CN genome-wide. Reframe the
  "strip before deposit" item as: **confirm capture-kit concordance (or rebuild the
  reference from a kit-matched panel-of-normals) before quantitative CNV use.**

**Arm-event frequencies are ROBUST to patient collapse** (new file
`output/wes_cnv_arm_freq_patient.csv`; HGS, 18 lines → 11 patients):

| Arm event | Per line | Per patient |
|---|---|---|
| 3q gain | 83% | **82%** |
| 8q gain | 67% | 73% |
| 19q gain | 44% | 55% |
| 20q gain | 78% | **91%** |
| 17p loss | 83% | **82%** |
| 13q loss | 50% | 64% |
| 10q loss | 28% | 36% |

Point-locus recovery (per line → per patient): 3q26/MECOM 89→91%, 3q26/SOX2
89→82%, 8q24/MYC 78→73%, 19q12/CCNE1 72→91%, 20q 72→82%, 17p/TP53 72→64%,
13q/RB1 56→73%, 10q/PTEN 33→45%.

**State the asymmetry explicitly:** *pseudoreplication distorts MUTATION
frequencies badly (CDK12 35%→18–27%) but CNV frequencies barely (3q gain 83%→82%,
17p loss 83%→82%)*, because HGSC arm-level events are near-universal trunk events
shared across all sublines, whereas a family-private mutation (the 3133 CDK12
frameshift) is one event inflated ×4 by subline counting.

## 4. Figures (regenerated; family tracks + tier encoding)

- **`figs/f_wes_oncoplot.pdf` + `reports/assets/f_wes_oncoplot.png`** — rebuilt with
  `ComplexHeatmap::oncoPrint`: (a) **patient-family annotation track** (the four
  3133 columns read as ONE patient); (b) **somatic-tier encoding** (fill = Tier1
  rust / Tier2 amber / Tier3 grey; box height = truncating-LOF vs missense/in-frame);
  (c) a **TMB track** (retained coding candidates/line) flagging the TOV21G
  hypermutator spike (rust bar, ~1,416 vs ~200 baseline); gene bars count
  **patients (of 16)**, per-line % suppressed.
- **`figs/f_wes_cnv.pdf` + `reports/assets/f_wes_cnv.png`** — autosome-only (chrX &
  chrY dropped), per-sample autosome-median-centred, **patient-family + subtype
  tracks** on the left, canonical HGSC event markers retained.
- Legacy `figs/07_oncoplot_drivers.pdf` and `figs/08_cnv_landscape.pdf` are also
  refreshed so existing references don't break; **report should switch F4 refs to
  `f_wes_oncoplot` / `f_wes_cnv`.**

## 5. Minor reconciliations

- **"BRCA2 68→9%":** 9% = 2/22 all-subtype lines (post-filter). By-subtype table
  shows 1/17 HGS = 5.9%. After tiering the **defensible somatic BRCA2 rate is 0**;
  report the 9%/5.9% as filtering-artefact residual, not HR-deficiency evidence.
- **TP53 denominators made explicit everywhere:** 18/22 all-subtype MAF lines (82%)
  = 12/16 patients (75%); **17/17 HGSC lines = 11/11 HGSC patients (100%)**.

## 6. Cross-workstream / deposition notes

- **MAF genome-build mislabel (deposition FYI, from `wes-signatures`):** the archived
  Mutect2 MAF `NCBI_Build` header reads **GRCh37**, but the data is actually **GRCh38**
  (verified via PoN/coordinates — a `vcf2maf` default mislabel). Correct the build field
  before deposition; all coordinates in the re-analysis are treated as GRCh38.
- **Hypermutator handled jointly with `wes-signatures`:** TOV21G is covered fully in
  `f_wes_hypermutation`; the oncoplot TMB track and the `context` column here are the
  driver-side tie-in (no duplicate figure).

---

## Report sections to update (`reports/01_multiomic_characterization_results.md`)

- **Abstract (l.13):** "TP53 in 17/17 HGSC" → "TP53 in 17/17 HGSC lines (**11/11
  independent patients, 100%**)".
- **WES methods (l.36):** add — frequencies re-tabulated per patient; somatic-
  confidence tiers replace the uninformative `germline_like_vaf`; FGA is
  autosome-restricted and chrX is excluded (sex-composition artifact, now
  quantified); median-centering caveat for FGA>0.7; unmatched public normals are a
  **validity** caveat.
- **Drivers & CN (l.55–56):** update figure refs to `f_wes_oncoplot`/`f_wes_cnv`;
  restate CDK12 as 3/11 patients (2/11 defensible); insert the patient-level CNV
  frequencies and the pseudoreplication-asymmetry sentence; FGA ordering →
  autosome values (HGS 0.62 > CC 0.37 > MC 0.32 > EC 0.23 > LGS 0.02) with the
  n=1 caveat; TOV81D FGA 0.073 → **0.021**.
- **Limitations (l.15, l.90–91):** add patient-level framing; BRCA1/2 germline
  caveat; defensible-BRCA2 = 0.
- **Deposition (l.100):** reframe the unmatched-normal item as a validity check.
- **F4 list (l.110):** `f_wes_oncoplot`, `f_wes_cnv`.

## Files (all under project root)

New: `output/wes_driver_freq_patient.csv`, `output/wes_driver_tiers.csv`,
`output/wes_cnv_arm_freq_patient.csv`.
Augmented (non-destructive): `output/wes_driver_freq_by_subtype.csv`,
`output/wes_cnv_fga.csv`, `output/wes_cnv_segments.csv` (adds `log2c_auto`).
Figures: `figs/f_wes_oncoplot.pdf`, `reports/assets/f_wes_oncoplot.png`,
`figs/f_wes_cnv.pdf`, `reports/assets/f_wes_cnv.png` (+ legacy names refreshed).
