# Delta — Hypermutation / MMR / MSI / mutational context (WES)

**Workstream:** `wes-signatures` · **Peer-review item:** §3.8 (TOV21G carries a mutation load >3× the next line, previously unflagged)
**Script:** `scripts/16_wes_signatures_msi.R` · **Date:** 2026-07-23
**Inputs (read-only):** `output/wes_mutations_filtered.csv` (PASS coding-nonsyn candidates), raw Mutect2 MAFs in `judy_archive/data/wes - old/mutect2/` (SBS substrate), `metadata/line_family_map.csv`
**Outputs:** `output/wes_mutation_load.csv`, `output/wes_msi_mmr.csv`, `output/wes_sbs_context.csv`; `figs/f_wes_hypermutation.pdf` + `reports/assets/f_wes_hypermutation.png`

---

## Headline

**TOV21G is the panel's sole hypermutator and is best described as a candidate MMR-deficient / MSI-high ovarian clear-cell carcinoma (OCCC) model — not a POLE-ultramutator and not a tumor-only artifact.** This is a *reuse feature* (rare MSI-high OCCC model, relevant to immunotherapy / MMR biology), not a QC defect. The call is **qualitative** (tumor-only WES; no matched normal) and rests on four converging, independent proxies plus a subtype-internal negative control.

## TOV21G exact load

| Metric | TOV21G | Panel median | 2nd-highest line |
|---|---|---|---|
| Coding non-syn variants (filtered MAF) | **1,416** | 206 | TOV2929D = 413 |
| Fold over median / over 2nd | **6.9× / 3.4×** | — | — |
| Robust outlier z (MAD-based) | **21.2** | 0 | 3.6 (next) |
| SNV / indel / MNV | 1,015 / 392 / 9 | — | — |
| **Indel fraction** (coding; exome-wide) | **0.28 / 0.30** | ~0.09 | 0.13 |
| **Ts/Tv** | **2.72** | ~0.8 | 0.57–1.21 range |
| Exome-wide PASS SNVs (SBS substrate) | 2,417 | — | — |

Hypermutator flag (robust-z > 5 **and** > 3× median) fires on **TOV21G only**; the definition is insensitive to the exact cutoff because TOV21G's robust-z (21.2) is ~6× the next line (3.6). Per-patient collapse (family map → **16 WES patients**): TOV21G is 1/16; no multi-line family approaches it (family load ranges 166–321). The **other clear-cell line, TOV3392D, is not elevated** (196 variants, indel fraction 0.08) — hypermutation is line-specific, not a CC-annotation or pipeline artifact.

## MMR / proofreading gene findings

- **No MMR-enzyme coding/splice mutation in TOV21G** (MLH1, MSH2, MSH6, PMS2 all wild-type in the filtered MAF).
- **No POLE mutation anywhere in the 22-line panel** (0 exonuclease-domain [codon 268–471] hits; POLD1 also clean). This alone excludes the POLE-ultramutator hypothesis.
- TOV21G's only MMR-locus hit is an **EPCAM frameshift** (HIGH impact, pop-AF 0; HGVSp reconstructed as `p.-2-3X`). EPCAM acts on MMR only *indirectly* (3′ deletions silence MSH2); a coding frameshift is not that canonical mechanism, and in tumor-only data its significance is **uncertain** (possible germline/artifact; VAF unresolved). Do not over-interpret.
- Absence of an MMR-enzyme lesion is **expected and consistent with MSI**: MSI in OCCC is usually driven by **MLH1 promoter hypermethylation** (epigenetic), which WES cannot see. The causal lesion is therefore *inferred*, not observed.
- Unrelated: TOV2881EP (HGS, not hypermutated) carries an MSH3 missense VUS — noise, not MSI.

## Mutational context (SBS-96) — computed

Signatures **were computed** (feasible). MutationalPatterns 3.18 + BSgenome.Hsapiens.UCSC.hg38 were installed; the 96-context matrix was built on **GRCh38** from exome-wide PASS SNVs and compared to COSMIC v3.2 by cosine similarity.

- **Genome build was verified GRCh38/hg38**, *not* the "GRCh37" written in the archived MAF `NCBI_Build` column (a spurious vcf2maf default). Evidence: Mutect2 PoN `1000g_pon.hg38.vcf.gz`; VCF `##contig chr1 length 248,956,422` (= GRCh38); driver hotspots at GRCh38 coordinates (KRAS G12 chr12:25,245,350; PIK3CA chr3:179.2 Mb; TP53 chr17:7.67 Mb). **→ deposition item: the MAF build header is mislabeled.**
- **TOV21G top COSMIC v3.2 matches: SBS6 (0.88), SBS44 (0.84), SBS15 (0.81), SBS20 (0.57) — all defective-MMR / MSI signatures.** POLE signatures (SBS10a/b/c/d, SBS28) all **≤0.32**. Spectrum is C>T-dominant (55%) with elevated T>C (20%).
- **Panel-wide specificity (negative control): TOV21G is the only MMR-d/MSI-like line.** Every other line's best match is a flat/clock signature (SBS5/40/3) with MMR-d cosine ≤0.70 and no POLE signal. Full per-line table in `output/wes_msi_mmr.csv`; 96×22 matrix in `output/wes_sbs_context.csv`.

## Verdict + evidence each way

**Verdict: candidate MMR-deficient / MSI-high OCCC** (qualitative; tumor-only).

**For genuine MMR-d/MSI:** (1) load 6.9× median / 3.4× next, line-specific (robust-z 21 vs 3.6); (2) indel fraction 0.28–0.30, ~3× panel — microsatellite indels are the MSI hallmark; (3) Ts/Tv 2.72 (transition-heavy) vs 0.6–1.2 elsewhere; (4) SBS6/44/15/20 top, uniquely in the panel; (5) TOV3392D (other CC) unaffected → not a subtype/pipeline artifact; (6) coherent with ARID1A-biallelic + PIK3CA H1047Y + KRAS G13C OCCC genotype (a subset of OCCC is MSI-high).

**Against POLE:** no POLE exonuclease mutation panel-wide; SBS10a/b cosine ≤0.32; POLE ultramutators are SNV-dominated, indel-poor, C>A/C>T at Tp*C*pT/Tp*C*pG — the opposite of TOV21G's indel-rich, SBS6-dominant profile.

**Against tumor-only artifact:** tumor-only germline inflation is ~constant across all lines (all ~206 median); the excess is TOV21G-specific. Germline contamination would not selectively raise indel fraction, transition fraction, *and* produce an MSI SBS signature in one line. Four independent metrics align — an artifact would not.

**Not proven (tumor-only limits):** the causal MMR lesion is inferred (MLH1 methylation is WES-invisible); MSI cannot be *definitively* called without a matched normal. **Recommended orthogonal confirmation:** MMR IHC (MLH1/PMS2/MSH2/MSH6), MSI-PCR or MSIsensor/MSIsensor-pro on the WES BAMs, and an MLH1 promoter-methylation assay.

## Proposed NEW report subsection (drop into F4 — Genomics and authentication)

> **### Hypermutation and MMR/MSI status — TOV21G is a candidate MSI-high clear-cell model (`figs/f_wes_hypermutation`)**
> One line, **TOV21G (clear cell), is a clear hypermutation outlier: 1,416 coding candidate variants — 6.9× the panel median and 3.4× the next line (TOV2929D, 413)** (robust-z 21 vs 3.6 next; the only line flagged at any reasonable cutoff, and the only outlier among 16 WES patients). The excess is enriched for indels (fraction 0.28 vs ~0.09 panel; exome-wide 0.30) and transitions (Ts/Tv 2.72 vs 0.6–1.2), and its GRCh38 SBS-96 spectrum matches the defective-mismatch-repair / MSI signatures **SBS6 (cosine 0.88), SBS44, SBS15, SBS20** — uniquely in the panel, with no POLE signal (SBS10a/b ≤0.32). No MMR-enzyme (MLH1/MSH2/MSH6/PMS2) or POLE coding mutation is present (consistent with the epigenetic — MLH1-methylation — mechanism typical of MSI-high clear-cell carcinoma, which WES cannot detect; TOV21G also carries biallelic ARID1A truncation, PIK3CA H1047Y, KRAS G13C). Because the WES is tumor-only, this is a **qualitative, converging-evidence call**: TOV21G is best described as a **candidate MMR-deficient / MSI-high OCCC model** — a valuable, rare reuse feature (immunotherapy / MMR biology) rather than a data defect — pending MMR IHC / MSI-PCR confirmation. The other clear-cell line (TOV3392D) is not hypermutated, confirming the signal is line-specific.

Also add to *Outstanding items / deposition*: **the archived MAF `NCBI_Build` header reads "GRCh37" but the calls are GRCh38** (verified from the PoN, contig lengths, and hotspot coordinates) — correct on deposition.

## Figure

`figs/f_wes_hypermutation.pdf` + `reports/assets/f_wes_hypermutation.png` — 4 panels: **A** ranked coding-load bar (fill = indel fraction; TOV21G flagged; median line); **B** MMR/proofreading gene-status tile (top-10 lines; only EPCAM in TOV21G); **C** TOV21G SBS-96 spectrum; **D** TOV21G cosine to diagnostic COSMIC v3.2 signatures (MMR-d vs POLE vs clock).

## Cross-workstream note (for `wes-drivers`)
TOV21G's hypermutation inflates its passenger burden; its non-canonical "driver" calls carry more noise than other lines. Its canonical OCCC drivers (ARID1A/PIK3CA/KRAS) remain real. Consider annotating TOV21G's hypermutator status alongside its driver tier in the oncoplot/driver tables.
