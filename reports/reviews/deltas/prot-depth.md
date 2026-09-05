# Delta — Proteomics depth: TMT ratio compression, bridge design, block-missingness

**Workstream:** `prot-depth` (peer-review §3.7)
**Script:** `scripts/19_proteomics_dynamic_range.R` (runs top-to-bottom; verified)
**New outputs:** `output/prot_dynamic_range.csv` (7,895 genes), `output/prot_block_missingness.csv` (8,429 proteins)
**Figure:** `figs/f_prot_compression.pdf` + `reports/assets/f_prot_compression.png` (4 panels)
**Read-only sources:** `output/prot_abundance_matrix.csv`, `rna_tpm.csv` + `tx2gene_ensembl_rel105.csv`, `integ_rnaprot_cor.csv`, `prot_bridge_cor.csv`, `judy_archive/.../{protein_relative_abundance,tmt.layout}.xlsx`, `metadata/samples.csv`, `scripts/05,06`.

Bottom line: all three §3.7 concerns are real and quantified, but none is a data-quality defect — they are intrinsic TMT properties that the resource should **document and guide reuse around**. The reported ~0.40 RNA–protein concordance is a compression-limited *ceiling*, not weak data.

---

## 1. TMT ratio compression — cross-line dynamic range (protein vs RNA)

**Method.** Same paired matrices as the concordance analysis (script 12): RNA = log2(TPM+1) collapsed to gene symbol; protein = log2 relative abundance; **8,212 shared genes × 30 shared lines** (RNA-only VOA6861, protein-only VOA14993 excluded). Per gene, cross-line spread (range = max−min, IQR, SD) computed over the **same lines** in both assays (n-matched: for a gene, RNA and protein use the identical protein-present line set), requiring ≥10 paired lines → 7,895 genes; 6,686 fully matched at all 30 lines.

**Result — protein spread is ~⅓ of RNA (n-matched, complete-30 subset):**

| statistic | RNA (median) | protein (median) | protein/RNA ratio |
|---|---|---|---|
| range (max−min) | 3.38 | 1.16 | **0.34** |
| IQR | 1.11 | 0.34 | **0.30** |
| SD | 0.85 | 0.27 | **0.32** |

Across the full ≥10-line set the median IQR ratio is **0.30**, and **98.0%** of genes have smaller protein than RNA spread (Fig A–B).

**FOLR1 exemplar — confirmed.** Protein **10.47–13.40** (range 2.93; 10.47–**13.60** across all 31 lines) vs RNA **0–9.48** (range 9.48); ratio 0.31. This matches the reviewer's exemplar (protein ~11–13.6 vs RNA ~0–9.5); the exact matrix minimum is 10.5 rather than 11.2.

**Caveat (stated, not hidden).** log2(TPM+1) and log2 TMT relative abundance are not a common physical scale, and part of RNA's larger range is genuine on/off regulation (TPM→0) that a TMT **ratio to a pooled reference cannot structurally represent** (the pool contains the protein, so protein abundance is bounded away from zero). The ratio therefore blends true co-isolation ratio compression with this scale/structural difference. Read it as *"the cross-line spread **as reported and used** is ~3× smaller at the protein level"* — which is exactly what governs downstream reuse — rather than as a pure estimate of MS2 co-isolation compression. The bridge-noise ruler below is the cleaner mechanistic number.

### Consequence (a): compression CAPS achievable RNA–protein concordance

- **Empirical, not asserted:** per-gene concordance rises monotonically with protein dynamic range. Binning the 7,893 genes with a concordance value by protein-IQR tercile → median Spearman **0.30 (low) / 0.40 (mid) / 0.53 (high)**. Spearman(protein IQR, concordance) = 0.33.
- **Noise floor:** the bridge replicates (§2) give a per-measurement technical-noise SD of **0.171 log2**. The median gene's cross-line protein SD is **0.287 log2 — only ~1.7× the noise floor.** A typical protein's biological variation sits close to noise, so ranks are noise-corrupted and rank concordance is mechanically attenuated.
- **Framing:** the reported per-gene median **~0.40** is a *compression-limited ceiling*, not poor data — which is precisely why CPTAC (0.45), CCLE (0.48), ProCan (0.42) land in the same band (they share this TMT/scale limit) and why the reproducibility ceiling is ~0.72. Well-measured, high-spread genes already reach ~0.53.

### Consequence (b): cross-line protein shortlists (ADC atlas) are less discriminating than RNA

Even for the clinically salient ADC targets — which are *well*-measured (concordance 0.59–0.90) — the protein axis is compressed 3–5× vs RNA (Fig C):

| target | RNA range | protein range | ratio | concordance |
|---|---|---|---|---|
| MSLN | 9.30 | 1.76 | 0.19 | 0.85 |
| TACSTD2 | 9.86 | 2.41 | 0.24 | 0.90 |
| CDH6 | 9.12 | 2.34 | 0.26 | 0.72 |
| CD276 | 3.61 | 1.00 | 0.28 | 0.69 |
| FOLR1 | 9.48 | 2.93 | 0.31 | 0.89 |
| ERBB2 | 3.16 | 0.98 | 0.31 | 0.59 |
| VTCN1 | 6.64 | 2.48 | 0.37 | 0.81 |
| SLC34A2 | 8.67 | 3.60 | 0.42 | 0.78 |

The protein axis has far less headroom to separate the top-expressing line from the pack, so protein-based "highest expresser" calls are more sensitive to noise/ties than the RNA equivalents. **Recommendation:** lead ADC shortlists with RNA and use protein as confirmatory (rank-consistency), rather than ranking on protein magnitude.

---

## 2. Bridge design + inter-plex normalization ACTUALLY implemented

**Inter-plex normalization is a common-reference (PIS) HUB — the bridge is NOT the normalization path.** From the TMT layout, **ch1 = "PIS" (Pooled Internal Standard) in every one of the 5 plexes**; the Morin lab expressed every channel as a log2 value normalized to its plex's PIS (→ "relative abundance", scale ~2–20 log2). Because the *same pooled standard* is present in all plexes, this anchors all five to one common reference — a hub for normalization. **Scripts 05/06 apply no further inter-plex normalization** — no median-centering, no bridge-based realignment; they consume the Morin matrix as-is (complete-case proteins for multivariate work; presence ≥50% for the reusable matrix). *This should be stated accurately in the methods: normalization is PIS-anchored, not bridge-anchored.*

**The daisy-chain bridge (ch10) is a separate technical-replicate QC device, and it IS a chain** (confirmed programmatically — every link connects consecutive plexes):

```
plex1 ─VOA10816─ plex2 ─TOV1369─ plex3 ─VOA3993*─ plex4 ─TOV3133G─ plex5
```

4 consecutive links; per-link Pearson **0.9914–0.9943** (Spearman 0.9932–0.9950). `*`VOA3993 is the **external Carey LGS** line (the plex3↔4 link). COV434 occupies plex1's ch10 as a *unique biological sample* (no prior plex to bridge from), not a replicate. Max separation plex1↔plex5 = 4 links.

**Calibrated correction of the review framing.** Because normalization is PIS-anchored, cross-plex *quantities* do **not** propagate through bridge links — the bridge is QC, not the alignment path, so distant-plex comparisons are not chained through 3–4 links as feared. **However**, the QC *evidence* for distant-plex comparability is itself chained: we hold direct technical replicates only for *consecutive* plexes (no single sample was run in both plex1 and plex5). Confidence that distant plexes are comparable therefore rests on (i) the shared PIS anchor and (ii) transitively chaining the consecutive bridge agreements. Per-link ~0.99 is excellent, but a common-reference bridge (same replicate in every plex) would have delivered *direct* distant-plex QC; the daisy-chain does not, and small per-link biases could in principle accumulate. Worth one honest sentence in the manuscript.

---

## 3. Structural block-missingness (per-plex)

TMT missingness is **all-or-nothing within a plex** (a protein seen in any line of a plex is seen in *all* its lines), so per-protein plex-coverage fully determines line-coverage. Of 8,429 proteins:

- **6,856 (81.3%) present in all 5 plexes** — the complete-case / `na.omit` set (matches the report's ~6,857).
- **1,573 (18.7%) absent from ≥1 whole plex.** Plex-coverage distribution: `5→6,856 · 4→512 · 3→366 · 2→360 · 1→264 · 0→71`. (The 71 with 0-coverage are all-NA across the 31 analysis lines — quantified only in external/bridge/spike channels; drop for reuse.)
- Median lines lost per non-complete protein: **13 of 31.**
- Per-plex quantified counts (lines): plex1 7,739 (7) · plex2 7,759 (6) · plex3 7,563 (5) · plex4 7,695 (6) · plex5 7,654 (7) — matches script 05.

**Bias on cross-line comparison.** A block-missing protein loses **whole plexes, i.e. whole line-blocks**, and because subtypes are distributed unevenly across plexes (e.g. plex3 carries the SCCOHT line BIN67 + 2 CC + 2 HGS; plex4 carries the EC line TOV112D, MC VOA8762, MMMT VOA5217 + 3 HGS; plex5 carries MMMT VOA5436 + 2 MC + 1 CC + 3 HGS), a missing plex can drop an entire rare-subtype representative. For these 1,573 proteins, "absent" is confounded with plex identity — present/absent and cross-subtype claims must be read as **plex-conditional, not biological**. The recommended reuse filter (presence ≥50% lines, 7,734 proteins) still admits proteins missing from 1–2 whole plexes, so `present_n_plex` (now in `prot_block_missingness.csv`) should be checked before cross-subtype comparisons.

---

## 4. Report sections to update

- **§ Proteomics QC (report ¶ "8,430 proteins / 146,830 peptides…").** Add: (i) block-missingness is *structural/per-plex* — 18.7% of proteins absent from ≥1 whole plex, coverage tiers as above; (ii) accurate normalization statement — inter-plex normalization is the **common PIS reference**, the ch10 bridge is a *daisy-chain technical-replicate QC* (per-link r≈0.99), with the chained-QC caveat. Current text implies the bridge is the batch-linking mechanism; correct to PIS-anchored.
- **§ RNA–protein concordance (report ¶ "30 lines / 8,212 shared genes… median 0.41").** Reframe the ~0.40 median as a **compression-limited ceiling**: protein cross-line spread ≈0.30× RNA and only ~1.7× the technical-noise floor; concordance rises with protein dynamic range (tercile 0.30/0.40/0.53). This turns "moderate concordance" from a worry into an expected assay property, consistent with the benchmarks sitting at the same value.
- **§ ADC-target atlas (report ¶ "Subtype-resolved RNA + protein…").** Add the shortlist caveat: protein ranges for ADC targets are compressed 3–5× vs RNA → protein shortlists are less discriminating; **lead with RNA, confirm with protein**. FOLR1 exemplar can cite the exact spans (protein 10.5–13.6 vs RNA 0–9.5).
- **Figure.** `f_prot_compression` slots into the proteomics-QC/concordance figure group (F2/F3) or as an ADC-atlas companion (F5). Panels: A paired per-gene IQR (RNA vs protein), B protein/RNA ratio distribution, C ADC-target ranges overlaid, D per-plex coverage (block-missingness).
