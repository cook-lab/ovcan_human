# Methods Reverse-Engineering: Judy's Analysis Notebooks

**Assessment target:** Three "All inclusive Final Code" R notebooks (RNA, Protein, WES) in `judy_archive/notebooks/`
**Purpose:** Document the complete methods as-run, cross-referenced to `judy_archive/output/` and `judy_archive/data/`, to inform a clean, reproducible re-implementation for a *Scientific Data* resource article.
**Prepared for:** Cook Lab clean re-analysis
**Date:** 2026-07-23

> Scope note: These three `.Rmd` files are the **only** analysis code in `judy_archive/` (confirmed by `find`). No proteomics-processing script, no LGSOC-specific script, and no HRD/CNV-scoring script is present, yet outputs from all of those exist in `output/`. Several consumed intermediate files (`tpm_matrix_human.csv`, `txi_human.rds`, `abundance_matrix_best.csv`, `protein_abundance_matrix_best.csv`, `*_LGS.csv`) are **produced by code not archived here**. Any faithful re-implementation must reconstruct those steps from scratch.

---

## 1. Per-notebook methods digest

### 1.1 RNA notebook — `All inclusive Final Code RNA.Rmd` (1,942 lines)

**Language/tools:** R. Packages: `tximport`, `tidyverse`/`dplyr`/`tidyr`/`readr`, `biomaRt`, `DESeq2`, `pcaExplorer` (topGO wrapper), `org.Hs.eg.db`, `GO.db`, `topGO`, `Rtsne`, `factoextra`, `cluster`, `matrixStats`, `pheatmap`, `ComplexHeatmap`, `ggrepel`, `RColorBrewer`, `viridisLite`, `progeny`, `singscore`, `fgsea` (via `gmtPathways`), `consensusOV`, `nichenetr` (**used but commented out in `library()` — see bug below**), `Hmisc`. No versions pinned; `sessionInfo()` not captured.

**Inputs consumed:**
- `../data/rna_seq/<sample>/abundance.h5` — kallisto per-sample output (37–39 sample dirs; passage suffix in dir name, e.g. `OV90p70`).
- `../data/human_sample_annotation_cleaned.csv` — subtype source (contains `Adenocarcinoma`, `Low-grade serous`, etc.).
- `../data/kallisto_qc.tsv` — per-sample QC (`n_processed`, `n_unique`, `p_pseudoaligned`).
- `../data/h.all.v7.4.symbols.gmt.txt` — MSigDB Hallmark v7.4.
- `../data/Han_NatComm_2020_Signalling_Metagenes_Cleaned.csv` — Han et al. 2020 signalling metagenes (mouse symbols).
- `../data/ovarian_cancer_HRD_peng_NatComm` — Peng et al. HRD table (**column names used as a gene list**, not a genomic score — see §4).
- Ensembl via `biomaRt` live query (`host='https://www.ensembl.org'`) for the transcript→symbol map.
- Mixed reads from `/Volumes/Eevee/projects/ovcan_human/output/differential_expression_*.csv` (hardcoded external drive) **and** `../output/*` — inconsistent within the same notebook.

**Outputs produced:** `../output/txi.rds`, `tpm_matrix.csv`, `annotation.csv`/`.tsv`, per-subtype `differential_expression_<subtype>_vs_Other.csv` (written to notebook CWD, no dir prefix, line 504), `*_signature_detailed.csv`, `Go_terms/*_goterms_all.csv`, `subtype_assignments_HGSC.csv` (to `/Volumes/Eevee/...`; the archived copy is `output/TCGA_subtype_assignments_HGSC.csv`), and many `../figs2/*.pdf`.

**Pipeline as-run:**
1. **t2g map:** `biomaRt::getBM` (transcript_id, gene_id, external_gene_name) → reduced to `target_id`→`ext_gene` (maps transcripts directly to symbol, skipping stable gene ID).
2. **Import:** `tximport(type="kallisto", txOut=FALSE, tx2gene=t2g, ignoreTxVersion=TRUE)`. `tpm_mat` drops row 1 ("sum of transcripts missing from tx2gene", line 110) and all-zero genes (line 111).
3. **Annotation:** `cell_line` = sample name split before "p" (passage). `subtype_dat$Subtype[cell_line %in% c("OV3331","OV90")] <- "High-grade serous"` — **manual reclassification from "Adenocarcinoma"** (lines 88–90). LGSOC (8 lines) hard-excluded (lines 134–137). Subtype releveled to a fixed 6-class order.
4. **Sample correlation:** Spearman on `log1p(TPM)`, `pheatmap` (mako, ward.D2, breaks 0.825–1).
5. **QC boxplots:** processed reads / unique / pseudoalignment %, jittered by subtype.
6. **PCA:** `log1p(TPM)` → zero-variance filter → top 2,000 most-variable genes → `prcomp(scale=FALSE, center=TRUE)`. Axis % variance **hardcoded** in labels ("PC1 (20.1% Variance)", line 328). tSNE: `Rtsne(pca$x[,1:5], perplexity=5)`.
7. **Cluster separation:** silhouette on subtype labels over PC1–PC10 (`cluster::silhouette`); Euclidean-distance `pheatmap` over PCs (**PC list has duplicated PC3–PC7, line 368 — copy-paste bug**); per-subtype PCA + distance heatmaps for CC/MC/HGS.
8. **Differential expression:** `DESeqDataSetFromTximport(txi_filtered, colData=annotation, design=~Group)` in a **one-vs-rest loop** (`Group = subtype vs "Other"`), `DESeq()`, `results(contrast=c("Group",subtype,"Other"))`, 6 subtypes. (Lines 483–489 contain leftover/duplicate `DESeqDataSetFromMatrix`/`~Subtype` scaffolding that never runs.)
9. **Signatures:** filter `padj<0.05 & log2FC>1` (up-only) per subtype → `*_signature_detailed.csv`. **Bug:** CC (line 703) and MC (line 609) read the DE CSV with `row.names=1` commented out, then set `Gene <- rownames(df)` → Gene becomes `"1","2",…` (integer strings), silently corrupting those two signatures (they later fail the `Gene %in% rownames(tpm)` filter and contribute ~0 genes downstream). HGS/EC/SCCOHT/MMMT use `row.names=1` and are correct.
10. **GO:** `pcaExplorer::topGOtable(fg = signature genes, bg = rownames(dds), ontology="BP", mapping="org.Hs.eg.db")` (topGO **elim** algorithm; `p.value_elim`). **Plotted terms are hand-picked by name** (5 curated terms per subtype for CC/EC/HGS/MC; top-5-by-p for SCCOHT/MMMT) — cherry-picking. `gene_list <- na.omit(gene_list)` (e.g. line 522) references an undefined object — a no-op/error line.
11. **Volcano plots:** per subtype, threshold `|log2FC|>0.6 & pvalue<0.05` (note: **unadjusted p** for coloring), label top-20 up by padj.
12. **Signature heatmaps:** top-10 and top-100 by padj per subtype, z-scored `log1p(TPM)` clipped ±2, `ComplexHeatmap` (no clustering, ordered by subtype).
13. **PROGENy:** `progeny(log1p(TPM), organism="Human", top=500)`, z-scored across samples; all-subtype and HGS-only heatmaps.
14. **singscore gene-set scoring:** `rankGenes(tpm_mat)` (raw TPM) → `simpleScore(upSet=…)` for: 10 Hallmark sets, PD-L1 (`CD274`), MHC-Ia (`HLA-A/B/C`), YAP targets (Wang 2018, hardcoded), chemoresistance/chemosensitivity (Sun 2019, hardcoded), **HRD = column names of the Peng table** (expression signature, not genomic), and Han et al. **"Targets" metagenes** (`convert_mouse_to_human_symbols`, nichenetr). Z-scored, `ComplexHeatmap`.
15. **ADC-antigen panels (hardcoded, used twice):** Panel A `FOLR1, MSLN, TACSTD2, ERBB2, TF, DPEP3, SLC34A2`; Panel B `FAP, EGFR, MSLN, AXL, ERBB2, CA9, MUC1, CEACAM6, CDH17`. Z-scored heatmap + size/colour dot plot faceted by subtype.
16. **Molecular subtyping (HGSC only):** `consensusOV::get.subtypes(method="consensusOV")` on symbol→ENTREZ-mapped `log1p(TPM)` for 15 HGSC lines → `subtype_assignments_HGSC.csv` (`DIF/IMR/PRO/MES_consensus`). A follow-up marker heatmap **retypes the assignments into a hardcoded `subtype_map` named vector** (lines 1892–1908) rather than joining the CSV — fragile and de-synchronizable.

### 1.2 Protein notebook — `All inclusive Final Code Protein.Rmd` (1,038 lines)

**Language/tools:** R. Packages as RNA minus DESeq2/progeny/singscore/consensusOV, plus `readxl`. Statistics are base-R `t.test` + `p.adjust`.

**Inputs consumed:**
- `../data/proteomics/protein_relative_abundance.xlsx` — **8,431 proteins × 56 columns**; per-sample values are **already TMT-normalized, log2-scale relative abundances** (~7–13). Columns also carry quant metadata (`Number of peptides`, `Number of unique peptides`, `qvalue`, `Sum PEP`, `Average SN`, `CV replicates`, `Npeptides_quant`).
- `../output/abundance_matrix_best.csv`, `../output/tpm_matrix_human.csv` — **provenance not in any notebook.**
- Protein t-test outputs and RNA `*_signature_detailed.csv` for integration.
- (`tmt.layout.xlsx`, `peptide_ratio.xlsx` present in `data/proteomics/` but **never read by the notebook**.)

**Outputs produced:** `filtered_abundance_table.csv`, `abundance_matrix.csv`, `protein_annotation_df.csv`, `ttest/ttest<subtype>_vs_OTHERS.csv`, `Protein signatures/ttest_sig/*_signature_protein.csv`, `*_merged_sig.csv` + `protein_rna_signatures.csv`, `protein_rna_genes.csv`, `figs2/*.pdf`.

**Pipeline as-run:**
1. **Column subset:** hardcoded 34-column keep-list (includes duplicate channels `VOA10816.1`, `TOV1369.1`, `TOV3133G.1`) → `Symbol` + samples. `make.unique` on symbols; **`na.omit` = listwise deletion** (drops any protein with a missing value in *any* sample — severe with TMT missingness). PCA uses `prcomp(scale=FALSE, center=TRUE)`; % variance hardcoded.
2. **Annotation:** a **second hardcoded `data.frame`** of CellLine→Subtype (lines 77–95), independent of the RNA `annotation.csv` (divergence risk).
3. **RNA–protein correlation:** intersect genes/lines; per-line Pearson `cor` on `log2(TPM+1)` vs protein; **only `TOV21G` is actually plotted** (single example, line 129).
4. **PCA / tSNE / silhouette:** same recipe as RNA.
5. **Protein signatures (t-test):** for a **single hardcoded `target_subtype` (="SCCOHT")**, per-gene **Student's t-test (`var.equal=TRUE`)** subtype-vs-rest, `log2FC = mean(group1) − mean(group2)` (data already log), BH `p.adjust`. **`annotation_table` is undefined in the notebook** (should be `protein_annotation_df_new`). Six per-subtype `ttest*` CSVs exist in `output/`, so this block was manually edited and rerun once per subtype. Signature = `padj<0.05 & log2FC>1`.
6. **RNA×Protein "merged signatures":** per subtype, `intersect()` of the protein t-test signature gene list and the RNA DESeq signature gene list → `*_merged_sig.csv`, combined `protein_rna_signatures.csv`; and a full-gene merge of protein t-test × RNA DE tables → `protein_rna_genes.csv` (used for log2FC-vs-log2FC scatterplots per subtype). Final merged-signature heatmaps rendered in both RNA and protein space.

### 1.3 WES notebook — `All Inclusive Final Code WES.Rmd` (265 lines)

**Language/tools:** R (`maftools`, `GenVisR`, `vcf2mafR`, `VariantAnnotation`, `data.table`, `dplyr`, `ggplot2`). Upstream calling done on an HPC with **Nextflow (nf-core/sarek-style; `work/<hash>/…recal.bam`)** + **Mutect2 + VEP** and **CNVkit 0.9.10 via Apptainer** (see `data/cnvkit wes - new/commands.txt`, `command_diagram.txt`).

**Inputs consumed:**
- `/Volumes/Eevee/…/data/wes/mutect2/<sample>/<sample>.maf` — hardcoded external drive; maps to `judy_archive/data/wes - old/mutect2/` (23 samples; each dir has `.maf` + `.mutect2.filtered_VEP.ann.vcf`).
- `/Volumes/Eevee/…/data/cnvkit/All samples/**/*.cns` — CNVkit segments; maps to `judy_archive/data/cnvkit wes - new/`.

**Outputs produced:** `figs2/HGSOC_oncoplot_paper_dec2.pdf`, `NONHGSOC_oncoplot_paper.pdf`, `cnvkit_heatmap_ALL_noY.pdf`, `cnSpec_heatmap_HGS_chr17_tp53.pdf`.

**Upstream calling (from `commands.txt`):** `cnvkit.py batch <tumor.recal.bam> --drop-low-coverage -n normal_samples/all_bams/*.bam -f GRCh38(GATK Homo_sapiens_assembly38.fasta) -t intervals_sorted.bed -a intervals.antitarget.bed --diagram --scatter`. `diagram` re-run at `-t 0.8`.

**Old vs new CNVkit reference (important — per PI):** there were two CNVkit runs.
- **"New" run (the trustworthy set, and the only one archived here):** the reference is built from **healthy diploid normal BAMs** via `cnvkit.py batch -n normal_samples/all_bams/*.bam` — a **pooled panel-of-normals reference**, which is CNVkit's recommended practice. The five recurring SRA BAMs (`SRR4039087/88/89/96/97`) are the normal panel (their target/antitarget coverage `.cnn` are copied into each sample dir; only `SRR4039087.bam` is physically retained in this archive copy). **Verified as a real, non-flat reference:** `All samples/*/reference.cnn` carries genuine varying `log2`/`depth`/`spread` values (a flat reference would be all-zero), and the file is **byte-identical (md5 `f6de5ce50095823c8fe9818acb0630ab`) across samples** — i.e. one shared pooled-normal reference built once and reused cohort-wide. This is a sound reference construction.
- **"Old" run:** used a **flat/pooled reference (or none)** and is the less trustworthy set. It is **not present in this archive** (no `.cns/.cnn` exist outside `cnvkit wes - new/`); flag it so no figure is regenerated from old, flat-reference calls.
- **Caveat that remains (minor):** the diploid normals are **not donor-matched** to the cell lines and are external/public samples — fine for cohort CN *profiling*, but per-line germline-vs-somatic CNV cannot be separated, and one should confirm the normals used the **same capture kit** as the tumor libraries. This is a much smaller issue than the matched-normal problem in mutation calling (§ below).

**Pipeline as-run (R):**
1. **Mutations:** `read.maf` per sample; `merge_mafs` into HGS cohort (**9 samples listed, but title says "12 lines"** — mismatch, lines 74–96) and non-HGS cohort (4 samples). `plotmafSummary`; `oncoplot(top=25, removeNonMutated=TRUE, draw_titv=TRUE)`; and a curated-gene oncoplot (`TP53, BRCA1/2, NF1, RB1, PTEN, PIK3CA, KRAS, CDK12, ATM, ATR, ARID1A, SMARCA4, NRAS, BRAF`).
2. **CNV:** reads from `/Volumes/Eevee/…/data/cnvkit/All samples` → maps to the archived **"new"** tree (`cnvkit wes - new/All samples/`), i.e. the **diploid-normal-reference (trustworthy) calls**. Reads all `.cns` (exclude `call.cns`/`bintest.cns`), keeps `log2` as `segmean`, **centers per sample by subtracting the median** (line 202), drops `chrY`. `GenVisR::cnSpec` (relative scale, blue-white-red, limits ±2) with a **hardcoded 13-sample order**; `GenVisR::cnFreq` for chr17 with a TP53 line at 7,668,421 bp. `.call.cns` (integer copy-number calls) are explicitly **not** used (viz is on median-centered log2 segment means).
3. **Provenance mismatch to resolve:** mutation MAFs are read from `…/data/wes/mutect2` → the archived **`wes - old/mutect2/`** tree, while CNV uses the **"new"** CNVkit tree. So the notebook combines **old-run Mutect2 calls with new-run CNVkit calls**. Confirm whether a "new" Mutect2 run exists that should replace the old MAFs before finalizing the oncoplots.
4. **No HRD, no mutational-signature, no driver-enrichment, no CN-burden statistics.**

---

## 2. Consolidated "Analyses Performed" table

| # | Analysis | Assay | Tool / function | Key params | Inputs | Outputs | Concern |
|---|----------|-------|-----------------|-----------|--------|---------|---------|
| 1 | Transcript→gene map | RNA | `biomaRt::getBM` | live Ensembl, unversioned | Ensembl | (in-memory `t2g`) | High (network/version) |
| 2 | Quant import | RNA | `tximport` (kallisto) | `txOut=F, ignoreTxVersion=T`, drop row 1 | `abundance.h5` | `txi.rds`, `tpm_matrix.csv` | Med |
| 3 | Subtype relabel + LGS exclusion | RNA | manual | OV3331/OV90 Adeno→HGS; drop 8 LGSOC | annotation csv | `annotation.csv/.tsv` | Med (document) |
| 4 | Sample correlation | RNA | `cor` Spearman + `pheatmap` | ward.D2, log1p | TPM | `spearman_cor_heatmap.pdf` | Low |
| 5 | Library QC | RNA | `ggplot2` | boxplots by subtype | `kallisto_qc.tsv` | QC pdfs | Low |
| 6 | PCA / tSNE | RNA | `prcomp` / `Rtsne` | top2000 HVG, scale=F; perplexity=5 | log1p TPM | `pca/tsne_plot_rna.pdf` | Med (hardcoded %, perplexity) |
| 7 | Cluster separation | RNA | `cluster::silhouette`, dist | PC1–10; labels=subtype | PCA | silhouette/euclid pdfs | Med (PC dup bug L368) |
| 8 | Differential expression | RNA | **DESeq2** one-vs-rest | `~Group`, Wald, default | `txi` counts | `differential_expression_*_vs_Other.csv` | **High (no replication, n=2 groups)** |
| 9 | Subtype signatures | RNA | filter | `padj<0.05 & log2FC>1`, up-only | DE csv | `*_signature_detailed.csv` | **High (CC/MC Gene bug)** |
| 10 | GO enrichment | RNA | `pcaExplorer::topGOtable` | BP, elim, org.Hs.eg.db | signatures | `Go_terms/*_all.csv` | Med (hand-picked terms) |
| 11 | Volcano plots | RNA | `ggplot2` | `|log2FC|>0.6 & p<0.05` (unadj) | DE csv | volcano pdfs | Med |
| 12 | Signature heatmaps | RNA | `ComplexHeatmap` | top10/100 padj, z ±2 | TPM+sig | `sigs_heatmap_*` | Low |
| 13 | Pathway activity | RNA | **PROGENy** | `top=500`, Human, z-score | log1p TPM | heatmaps + `pathway_summary_LGS.csv`* | Med |
| 14 | Gene-set scoring | RNA | **singscore** `simpleScore` | rank TPM, upSet only | Hallmark/Han/HRD/ADC/YAP… | heatmaps + `geneset_score_matrix_LGS.csv`* | Med |
| 15 | HRD (expression) | RNA | singscore of Peng **column names** | up-set | Peng table | (score row) | **High (not genomic HRD; mislabel risk)** |
| 16 | Han signalling metagenes | RNA | singscore + `convert_mouse_to_human_symbols` | Targets sets only | Han csv | heatmap | Med (nichenetr commented out) |
| 17 | ADC-antigen panels | RNA+Prot | z-score heatmap/dotplot | 2 hardcoded panels | TPM / protein | antigen pdfs | Low (descriptive) |
| 18 | Molecular subtyping | RNA | **consensusOV** `get.subtypes` | `method="consensusOV"`, ENTREZ | log1p TPM (15 HGSC) | `TCGA_subtype_assignments_HGSC.csv` | **High (tumor-tool on cell lines; hardcoded remap)** |
| 19 | Protein matrix build | Protein | `readxl` + `na.omit` | listwise deletion; 34 cols | `protein_relative_abundance.xlsx` | `abundance_matrix.csv` | **High (missingness handling)** |
| 20 | Protein PCA/tSNE/silhouette | Protein | `prcomp`/`Rtsne`/`silhouette` | scale=F, perplexity=5 | abundance | `*_protein.pdf` | Med |
| 21 | RNA–protein correlation | Protein | `cor` (Pearson) | per-line; **only TOV21G plotted** | TPM+protein | `rna_vs_protein_TOV21G.pdf` | Med (single example) |
| 22 | Protein DE | Protein | **Student's t-test** (`var.equal=T`) | one-vs-rest; BH | abundance | `ttest/ttest*_vs_OTHERS.csv` | **High (n=2, undefined var, manual loop)** |
| 23 | Protein signatures | Protein | filter | `padj<0.05 & log2FC>1` | ttest csv | `*_signature_protein.csv` | Med |
| 24 | RNA×protein merged sigs | Protein | `intersect`/`merge` | by Gene | RNA+prot sigs | `protein_rna_signatures.csv`, `*_merged_sig.csv` | Med |
| 25 | Somatic mutation calling | WES | **Mutect2** (sarek) + VEP → MAF; **"old" run** | **tumor-only, no matched normal / no PoN documented** | `.recal.bam` | `wes - old/mutect2/<sample>.maf` | **High (germline vs somatic in cell lines)** |
| 26 | Oncoplots | WES | `maftools` | HGS(9→"12") vs non-HGS(4); curated genes | MAFs | oncoplot pdfs | Med (n mismatch; old-run MAFs) |
| 27 | CNV calling | WES | **CNVkit 0.9.10** batch; **"new" run** | **pooled diploid-normal reference** (5 SRA normals, verified non-flat, shared); GRCh38; capture BED | BAMs | `cnvkit wes - new/*.cns/.cnr/.call.cns` | Low–Med (sound ref; normals unmatched/public) |
| 28 | CNV visualization | WES | `GenVisR` cnSpec/cnFreq | per-sample median-centered log2, drop chrY | "new" `.cns` | CNV heatmap pdfs | Med (uses log2 not integer `.call.cns`) |

\* `*_LGS.csv` and `_human`/`_best` files are consumed/exist but the generating code is **not in these notebooks**.

---

## 3. Sample structure (critical for a resource paper)

- **RNA (post-LGS exclusion, n=31 lines, 1 sample each — no replicates):** Clear cell **7**, Endometrioid **2**, High-grade serous **15**, MMMT **2**, Mucinous **3**, SCCOHT **2**. Three subtypes at **n=2**.
- **Protein (n≈34 columns incl. 3 duplicate channels):** subtype counts similar; LGS present in raw data (`tmt.layout.xlsx`) but excluded in-notebook.
- **WES:** ~23 lines with MAFs; oncoplot cohorts of 9 (HGS) and 4 (non-HGS).
- **Passage mismatch across assays** is real and unaddressed, e.g. OV90 = RNA `p70`, WES/CNV `P63`; TOV3133G = RNA `p61`, protein channel, WES `P65`. Passage is a plausible confounder and must be tracked in the resource metadata.
- **Proteomics batch/site structure (from `tmt.layout.xlsx`, ignored by the notebook):** multiple **TMT plexes**, a **`PIS` pooled-internal-standard bridge channel per plex**, and **two source sites (`Mes-Masson` [CHUM] vs `Carey`)**. Plex and site are unmodeled batch variables.

---

## 4. Clarifications the task specifically asked about

- **TCGA subtype assignment = consensusOV.** The four labels in `TCGA_subtype_assignments_HGSC.csv` (`DIF/IMR/PRO/MES_consensus`) are consensusOV's TCGA/Verhaak-style calls, run on 15 HGSC lines only. Assignments are then **manually retyped** into a hardcoded vector for the marker heatmap.
- **Differential expression tool = DESeq2** (RNA, one-vs-rest, `~Group`) and **Student's t-test** (protein). Not limma/edgeR. **No replicates and no batch/passage/site covariates** in either model.
- **HRD "scoring" is expression-based**, via singscore of the *column names* of `ovarian_cancer_HRD_peng_NatComm`. **There is no genomic HRD (LOH/TAI/LST) computed anywhere**, and the WES notebook does not touch HRD despite the data folder's presence. Do not describe this as genomic HRD.
- **Mutect2:** tumor-only (cell lines have no matched normal); no panel-of-normals or explicit germline filtering is documented. The notebook uses the **"old"-run MAFs** (`wes - old/mutect2/`). This is where the unmatched-normal problem genuinely bites — germline and cell-line-private variants can masquerade as somatic.
- **CNVkit:** the **"new" run** built a proper **pooled diploid-normal reference** (`cnvkit.py batch -n <normal BAMs>`; verified non-flat and shared cohort-wide) and is the trustworthy set the notebook actually visualizes; the **"old" run** used a flat/pooled reference and is not archived. `.call.cns` integer calls exist but the figure uses per-sample median-centered `log2` segment means. Do not regenerate CNV figures from "old"/flat-reference calls.

---

## 5. Prioritized re-implementation & concerns

### KEEP (sound in principle; re-run cleanly with rigor)
- kallisto→`tximport`→DESeq2 backbone (proper counts+offset path); **but** re-scope the DE question (see below).
- PROGENy, singscore gene-set scoring, consensusOV, maftools/oncoplots, CNVkit + a CNV visualization — all appropriate tools for a resource. Retain as descriptive layers with explicit caveats.
- MSigDB Hallmark v7.4, Han metagenes, ADC panels, PROGENy — reusable, but move all gene sets to versioned, documented resource files.

### FIX (methodological corrections required)
1. **Reframe/limit differential expression.** One-vs-rest DESeq2/t-tests across **unreplicated cell lines with n=2 groups** cannot support generalizable "subtype markers." Either (a) present as **descriptive line-level contrasts** with heavy caveats, (b) drop formal p-values for n≤2 subtypes, or (c) restrict statistical DE to HGS (n=15) only. Never call these validated subtype signatures.
2. **Proteomics batch/site/plex.** Model or at least diagnose **plex and site (Mes-Masson vs Carey)**; use the `PIS` bridge channels; assess the duplicate channels (VOA10816/TOV1369/TOV3133G) as reproducibility controls instead of independent samples. Replace `na.omit` listwise deletion with a presence threshold + principled missing-value handling.
3. **CC and MC signature `Gene` bug** (RNA lines 609, 703): fix `row.names`/`Gene` assignment so symbols aren't overwritten by row indices. Re-generate all downstream merged signatures and heatmaps.
4. **Stop hardcoding derived values:** PCA % variance in axis labels; the consensusOV `subtype_map`; hand-picked GO terms; the duplicated PC list (L368). Derive everything programmatically; if curating GO terms, disclose it.
5. **WES cell-line caveats (mutation vs CNV are different):**
   - *Mutations (the real problem):* tumor-only Mutect2 with **no matched normal** conflates germline and cell-line-private variants with somatic events. Apply gnomAD/germline + PoN-style filtering, restrict oncoplots to well-supported drivers, and state the limitation prominently. Confirm whether a **"new" Mutect2 run** should replace the archived "old" MAFs (the notebook currently mixes old mutations with new CNV). Resolve the HGS "9 vs 12 lines" discrepancy.
   - *CNV (largely fine — keep):* use only the **"new"** CNVkit calls (pooled diploid-normal reference, verified sound). Do **not** use "old"/flat-reference calls. Minor follow-ups: confirm the normal panel used the same capture kit as the tumors; consider using integer `.call.cns` (or `-t 0.8` thresholds) rather than raw log2 for the CN heatmaps.
6. **Unadjusted p-values** used for volcano thresholds — switch to padj or clearly label.
7. **consensusOV on bulk cell lines** is off-label (designed for HGSOC *tumors*). Keep as exploratory; do not over-interpret IMR (immunoreactive = stroma/immune signal absent in pure lines).

### DROP or REBUILD (not reproducible as archived)
- **Hardcoded `/Volumes/Eevee/...` paths** and mixed `../output` vs external-drive reads throughout RNA + WES — rebuild with project-relative paths / a config.
- **Live `biomaRt` query** — replace with a pinned annotation package (e.g. versioned `EnsDb`/`org.Hs.eg.db`) so the transcript→gene map is frozen.
- **Undocumented intermediates** (`tpm_matrix_human.csv`, `txi_human.rds`, `abundance_matrix_best.csv`, `protein_abundance_matrix_best.csv`, all `*_LGS.*`) — their generating code is **missing**; reconstruct and document, or drop.
- **Proteomics upstream normalization is entirely undocumented** (empty `Readme.md`; `.xlsx` already "relative abundance"). Obtain and record the facility's TMT pipeline (search/normalization/rollup/reference) or re-process from `peptide_ratio.xlsx` transparently — mandatory for *Scientific Data*.
- **Protein t-test loop** that requires manual `target_subtype` edits and references an undefined `annotation_table` — rebuild as a single reproducible function over all subtypes.
- **nichenetr `convert_mouse_to_human_symbols`** used while the `library(nichenetr)` line is commented out (won't run as written) — fix dependency loading or replace with a documented ortholog map.

### ADD (for a defensible resource)
- A single source-of-truth sample sheet keyed by cell line with **subtype, site, passage-per-assay, plex/channel**, and cross-assay availability.
- `sessionInfo()`/`renv` lockfile and pinned reference/gene-set versions.
- Genomic **HRD/CN burden and (optionally) mutational signatures** from WES if HRD is to be claimed — the current "HRD" is expression-only.
- End-to-end runnable scripts (the RNA notebook is **not runnable top-to-bottom** as written — e.g. the commented-out `colnames(tpm_mat) <- annotation$cell_line` at line 112 breaks all later `cell_line`-indexed subsetting; it clearly relied on manual, out-of-order execution).

---

*Cross-referenced against `judy_archive/output/` (Differential expression, Go_terms, Protein signatures, RNA detailed signatures, ttest, and matrix/annotation CSVs) and `judy_archive/data/` (rna_seq kallisto dirs, proteomics xlsx + tmt.layout, cnvkit wes - new commands, wes - old mutect2 MAFs). All line numbers refer to the `.Rmd` files as archived.*
