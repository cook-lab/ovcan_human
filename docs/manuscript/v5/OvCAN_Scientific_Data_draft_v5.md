# A multi-omic resource of human ovarian cancer cell models spanning common and rare histotypes

*Working draft v5 for internal review. Text in square brackets identifies information required before submission. This draft is a new document and does not replace the existing manuscripts.*

[TO ADD: complete author list, affiliations, ORCIDs and corresponding-author details]

## Abstract

Ovarian cancer comprises molecularly distinct histotypes, but the cell models used to study these diseases have been characterised unevenly across laboratories and assay platforms. Here we describe transcriptomic, proteomic and whole-exome data for 42 human ovarian cancer cell models derived from 34 patients and representing seven histotypes. RNA sequencing and tandem mass tag proteomics profile 31 models each; whole-exome sequencing provides mutation calls for 22 models and copy-number profiles for 23; 13 models have all three layers. All data were reprocessed through version-controlled workflows, linked by a common model identifier and accompanied by per-model quality metrics, donor-family annotations, feature-level missingness and analysis-ready matrices. Technical validation assesses RNA-sequencing quality, cross-plex protein agreement, recovery of expected histotype-associated expression and genomic features, transcript-protein concordance and consistency with external cell-line references. The resource supports the selection and comparison of ovarian cancer models by histotype, genotype, transcript abundance and protein abundance.

## Background & Summary

Ovarian cancers are classified into histotypes with distinct cells of origin, genomic alterations and clinical behaviour. High-grade serous carcinoma is characterised by near-universal disruption of *TP53* and extensive copy-number alteration, whereas clear cell, endometrioid, mucinous, low-grade serous and small cell carcinomas follow different molecular trajectories^1^. These differences make model selection a biological decision rather than a matter of convenience.

Cell lines remain central to mechanistic studies and preclinical screening, but ovarian cancer models have often been used without adequate consideration of histotype or provenance. Genomic comparisons have shown that several frequently used lines are poor models of high-grade serous carcinoma^2^, and systematic characterisation has revealed discordance between historical labels and molecular features^3^. The problem is greatest for rare histotypes, for which only a small number of models are available and measurements are distributed across studies. As a result, researchers may need to compare a transcript measurement from one study with a protein or genomic measurement generated elsewhere under different conditions.

Several studies have improved this foundation. Uniform molecular profiling has been reported for ovarian cancer cell-line panels^3,4^, including proteomic analysis of predominantly common histotypes^5^ and matched genomic, transcriptomic and proteomic profiling of low-grade serous models^6^. These resources establish the value of measuring models consistently, but no single dataset covers the present panel across common epithelial histotypes, carcinosarcoma and small cell carcinoma of the ovary, hypercalcaemic type.

Here we describe a collection of 42 ovarian cancer cell models derived from 34 patients at three Canadian centres. The resource combines RNA sequencing, tandem mass tag proteomics and whole-exome sequencing, with 13 models represented in all three layers (Fig. 1). Thirteen models are sublines from five patients; donor-family membership and a patient-representative flag are therefore included with the data. We provide processed matrices, per-feature and per-model quality metrics, candidate driver calls, copy-number profiles and a master model table that connects every record through a common identifier. Together, these data allow researchers to compare models using measurements generated and processed within a common framework.

## Methods

### Cell models and study design

The resource contains 42 models marked `analysis_include = Y` and `provenance = generated` in `metadata/samples.csv` (Table 1). Twenty-nine models originated from the Centre de recherche du Centre hospitalier de l'Université de Montréal, 12 from BC Cancer-OVCARE and one from the Ottawa Hospital Research Institute. Thirteen additional low-grade serous samples in the source metadata were generated in a separate study^6^ and were excluded from the present resource. One excluded proteomic sample, `VOA3993`, was used only as an inter-plex bridge.

Histotype annotations were harmonised as high-grade serous, low-grade serous, clear cell, endometrioid, mucinous, carcinosarcoma or small cell carcinoma of the ovary, hypercalcaemic type. Conflicting historical annotations were retained in the metadata. In particular, `TOV112D` was historically described as endometrioid but was subsequently reassigned to dedifferentiated ovarian carcinoma with SWI/SNF loss^7^.

To account for related sublines, we derived patient-family membership and selected one representative model per patient using a deterministic score based on assay coverage. The resulting resource contains 34 patient representatives, including 28 with RNA data. We used models as the unit for descriptive displays and model-selection tables, but used patient representatives for marker inference and sensitivity analyses.

[TO ADD: model-specific derivation references and complete culture conditions for each contributing centre, including medium, serum, supplements, incubator conditions and passage at harvest.]

### RNA sequencing and quantification

[TO ADD: RNA extraction method, library-preparation kit, sequencing instrument and paired-end read length.] The libraries were sequenced in one run. We quantified transcripts with kallisto 0.46.0 against an index containing 185,299 targets^8^ and summarised transcript estimates to genes with tximport 1.36.1^9^. The transcript-to-gene map was obtained from the Ensembl release 105 archive and is deposited as `tx2gene_ensembl_rel105.csv`.

The kallisto index predates the deposited transcript map. Consequently, 3,529 index targets were absent from the map and were excluded during gene-level summarisation. These transcripts accounted for a median of 2.22% of sample TPM. Their contribution varied with histotype and contributing centre, although these factors are correlated in the study design. We therefore deposit the exact map used to create the matrices and recommend it for reproduction. Re-quantification against an index built from the same annotation remains the appropriate sensitivity analysis.

Gene-level TPM values were retained for 39,568 genes. For count-based analyses, we retained genes with at least 10 estimated counts in at least two models, leaving 22,544 genes, and calculated variance-stabilised values with DESeq2 1.48.2^10^. We recorded both the number of fragments processed by kallisto and the sum of assigned gene-level counts because these quantities describe different stages of the workflow.

### Tandem mass tag proteomics

Protein abundance was measured in five 11-plex tandem mass tag sets. Each set contained a pooled internal standard, and four samples were re-run in adjacent sets to connect the five plexes. The bridge samples differed between links and therefore assess cross-plex agreement across four representative samples rather than repeated measurement of a single common reference.

[TO ADD: sample preparation, digestion, tandem mass tag labelling, fractionation, liquid chromatography, mass spectrometer and acquisition mode, database-search software and version, sequence database, enzyme specificity, modifications, mass tolerances, peptide-spectrum-match and protein false-discovery thresholds, reporter-ion interference filtering, peptide-to-protein roll-up and pre-existing normalisation. Confirm the complete channel map and the definition of the vendor-reported coefficient of variation.]

The source table contained 8,430 protein rows. We removed three rows without a gene symbol, leaving an 8,427-feature matrix. When more than one row shared a symbol, we selected a representative using measurement completeness, peptide support and q-value; non-representative rows were retained with a `SYMBOL|UNIPROT` identifier. Protein values are reported as log2 abundance relative to the pooled internal standard. Missing values were retained because protein detection was structured by plex.

To estimate bridge repeatability, we calculated the difference between measurements of each bridge sample in adjacent plexes. Assuming equal and independent measurement error, the per-measurement standard deviation is the standard deviation of the paired differences divided by √2. We converted this value to a lognormal coefficient of variation and evaluated it by mean-abundance decile.

### Whole-exome sequencing and variant processing

Whole-exome data were processed against GRCh38 with nf-core/Sarek^11^, and small variants were called with GATK Mutect2 in tumour-only mode^12^. No matched normal was available. [TO ADD: nf-core/Sarek and GATK versions, library-preparation and sequencing details, exome capture kit and target BED.]

We filtered variant records in three stages: Mutect2 `FILTER == PASS`; maximum population allele frequency no greater than 0.001 across available gnomAD and 1000 Genomes fields; and membership in the maftools non-synonymous variant classes. Protein-change strings were reconstructed when the source field was empty and were flagged to distinguish canonical from non-canonical representations. Retained variants remain candidate coding alterations because tumour-only calling cannot reliably distinguish rare germline from somatic events.

We assigned candidate alterations in 14 ovarian cancer genes to three evidence tiers based on variant class and prior gene-level knowledge. The tiers prioritise calls for reuse but are not a validated somatic classifier. Tier 3 variants were deposited but were not used to support derived molecular calls or frequency statements.

Copy number was estimated with CNVkit 0.9.10^13^ using five public healthy exomes as a pooled reference. We centred each profile on its autosomal probe-weighted median and defined arm-level gain or loss as an absolute log2 copy ratio greater than 0.20 across more than half of an arm. Fraction of genome altered was calculated across autosomes. These profiles provide relative total copy number; they do not provide allele-specific copy number, loss of heterozygosity or a genomic homologous-recombination-deficiency score. Compatibility between the tumour and reference capture designs cannot be confirmed until the target BED is recovered.

### Analyses used for technical validation

We used principal component analysis to evaluate the structure of the RNA matrix. We transposed the 2,000 most variable genes, centred them and did not scale them. We compared cross-model spread between layers using the interquartile range and standard deviation of each gene measured in both datasets.

Histotype and contributing centre were partially confounded. To separate the components that could be estimated, we fitted joint models for each principal component and partitioned the explained variance into unique histotype, unique centre and shared contributions. We compared unique components with 1,000 label permutations. We also repeated the RNA analysis using one representative per patient.

To evaluate canonical expression patterns, we summarised established histotype markers across models and calculated effect sizes on the 28 RNA patient representatives. A marker was scored as recovered when the intended histotype was among the two highest group means for an expected high marker, or the two lowest for an expected loss marker; expected high markers also had to exceed 1 log2(TPM + 1). We evaluated the total recovered count using 20,000 joint permutations of histotype labels. This procedure preserves marker correlation and group size.

Transcript-protein concordance was calculated for the 30 models with both layers. For each gene quantified in at least 10 models, we calculated Spearman correlation across models and adjusted asymptotic two-sided P values with the Benjamini-Hochberg method. We repeated the analysis on the 27 dual-layer patient representatives.

For external comparison, we matched five models to DepMap Public 24Q4 expression profiles^14^. We compared each model with 67 DepMap ovarian models over 2,000 variable genes selected from the shared expression matrix. Cellosaurus records were retrieved by exact name^15^; these records provide external reference profiles and do not authenticate the stocks analysed here.

All analyses were performed in R 4.5.2 with Bioconductor 3.21. The project seed was 1234. The package environment is recorded in `renv.lock`, and `scripts/run_all.sh` defines the analysis order. Analysis scripts write intermediate checks and the tables used in each figure.

### Generative artificial intelligence assistance

Large-language-model coding assistants were used during code review, scripted analysis development and manuscript drafting. Numerical results were generated by the deposited R workflows rather than by a language model. The authors reviewed the analyses, verified the reported outputs and take responsibility for the final manuscript. [CONFIRM wording and tools used against the journal policy before submission.]

### Ethics

The models analysed here are pre-existing cell lines, and no new human material or identifiable participant data were collected for this study. [TO ADD: research ethics board approvals and consent framework under which the source tumours were originally collected, or the relevant institutional determination. Complete the Scientific Data Human Data Checklist.]

## Data Records

The resource is organised around `cell_line`, which is the shared identifier in all matrices and annotation tables. Table 2 lists the principal records, their units and intended repositories. The complete per-model inventory is provided in `supplement_per_line.csv`, with assay coverage, patient-family membership, quality metrics, candidate driver calls, fraction of genome altered, external-reference status and annotation-support fields.

The RNA record contains unfiltered gene-level TPM, filtered estimated counts, variance-stabilised values and per-model quality metrics. The protein record contains the 8,427-feature relative-abundance matrix, feature annotations, peptide support, plex presence and the zero-plex flag. The exome record contains filtered candidate variants for 22 models and copy-number segments and summary values for 23 models. Data dictionaries define each column and distinguish model, patient, feature and assay-specific denominators.

Protein missingness is structural by tandem mass tag plex. The files therefore record the number of models and plexes in which each feature was measured; missing values should not be interpreted as absent protein or imputed under a missing-at-random assumption. The family map similarly records which lines derive from the same patient so that population summaries can be calculated at the patient level.

[TO ADD: repository accessions and anonymous reviewer links. Raw RNA and exome reads should be deposited in GEO/SRA, raw proteomic spectra and search output in PRIDE through ProteomeXchange, and processed matrices, metadata, data dictionaries, checksums and the standalone browser in a persistent general-purpose repository.]

## Technical Validation

### RNA and protein measurements passed assay-level quality checks

To assess RNA-sequencing quality, we compared pseudoalignment, processed fragments and gene detection across the 31 models. Median pseudoalignment was 91.1%, median depth was 64.3 million paired-end fragments and a median of 20,119 genes was detected. No model was an outlier across these measures (Fig. 2a,b). Models contributed by BC Cancer detected approximately 5% more genes than those from the Centre de recherche du CHUM despite lower pseudoalignment. This difference is associated with contributing centre and should be modelled in depth-sensitive reuse analyses.

We next evaluated proteomic identification, missingness and cross-plex agreement. The final matrix contains 8,427 features supported by 146,830 quantified peptides. Every model within a plex has the same feature count, confirming that most missingness was introduced at the plex level rather than by individual models (Fig. 2d). Across the four bridge links, mean differences ranged from -0.028 to 0.009 log2, with per-measurement repeatability of 10.4-13.2% coefficient of variation.

Protein precision depended strongly on abundance. The bridge-derived coefficient of variation decreased from 21.3% in the lowest abundance decile to 5.2% in the highest, in parallel with the vendor-reported values (Fig. 2e). We therefore provide abundance-stratified precision estimates rather than a single threshold for all proteins. Together, these checks support comparison across plexes while identifying protein abundance and plex membership as necessary inputs to interpretation.

### Expression profiles recover histotype-associated structure

To determine whether the expression matrices retained expected ovarian cancer structure, we examined the leading RNA components. Models separated by recorded histotype, and histotype retained a unique contribution after joint modelling with contributing centre (Fig. 3a,b). On RNA principal component 1, the adjusted unique histotype component was 0.393, compared with -0.025 for centre. Repeating the analysis with one model per patient increased the unique histotype component to 0.523, indicating that the separation did not depend on duplicated sublines.

This design does not fully distinguish histotype from centre. All high-grade serous RNA models originated from one centre, and the shared variance component on principal component 1 was 0.311. The analysis therefore shows that centre adds little information beyond histotype, but cannot exclude culture or processing effects aligned with histotype. We consequently interpret the observed structure as histotype-associated rather than histotype-specific.

We next evaluated established histotype markers using patient representatives. Fifteen of 25 markers were recovered in the expected group, compared with a permutation mean of 6.8 and a 95th percentile of 12 (*P* = 0.0046; Fig. 3e,f). The result supports the biological coherence of the matrix without implying that every small histotype group is independently validated.

Finally, we assessed agreement between the two expression layers. Across 7,894 genes, the median transcript-protein Spearman correlation was 0.397, and the median within-model correlation was 0.408 (Fig. 3c). These values were maintained after selecting one model per patient and are comparable with correlations reported for other cancer cell-line proteomic resources^16,17^. Protein abundance varied less across models than transcript abundance, consistent with isobaric ratio compression and post-transcriptional regulation. The concordance estimates therefore support joint use of the matrices while preserving the distinction between transcript and protein measurements.

### Genomic profiles recover expected ovarian cancer features

To evaluate the tumour-only variant record, we first examined the effect of sequential filtering. Across 22 models, 557,392 Mutect2 records were reduced to 15,692 PASS variants and 6,036 rare coding candidates (Fig. 4a). Filtering removed widespread candidate calls in genes such as *ATM* and *ATR* while retaining *TP53* alterations in all 17 high-grade serous models. At the patient level, *TP53* was altered in all 11 high-grade serous patients with exome data, providing the primary positive control for the mutation record^1^.

Copy-number profiles showed the expected contrast between copy-number-high high-grade serous models and the low-grade serous model. Median autosomal fraction of genome altered was 0.635 among 18 high-grade serous models and 0.021 in the low-grade serous model. Common high-grade serous events, including 20q and 3q gain and 17p and 13q loss, were present across the patient-level profiles (Fig. 4d).

These genomic records are intended for candidate prioritisation and coarse copy-number comparison. Without matched normals, retained sequence variants cannot be assumed to be somatic. Likewise, the copy-number profiles are relative, total-copy estimates generated with pooled external normals of unconfirmed capture compatibility. The deposited tiers and threshold-sensitivity tables make these constraints explicit while preserving the expected ovarian cancer positive controls.

### External comparisons support a subset of model identities

To compare the expression data with an independent reference, we evaluated the five models with namesakes in DepMap. Each model matched its namesake at rank 1 among 67 ovarian models, and all five matches were reciprocal best matches (Fig. 4b). The smallest separation from a non-self model occurred between `BIN67` and `COV434`, two models of small cell carcinoma of the ovary, hypercalcaemic type, rather than between models of unrelated histotype.

The multi-omic profiles were also consistent with recognised SWI/SNF alterations. In `BIN67`, *SMARCA4* transcript abundance was near the middle of the panel whereas protein abundance was second-lowest, illustrating why the two expression layers should not be treated as interchangeable. `TOV112D` carried a truncating *SMARCA4* candidate together with low SMARCA4 protein, low *SMARCA2* RNA and *TP53* R175H, consistent with its published reassignment^7^.

No short-tandem-repeat or mycoplasma testing was performed on the stocks used for this study. Cellosaurus contains an external reference STR profile for 30 of the 42 model names, but these records do not establish a match to the analysed stocks. The DepMap comparisons and molecular profiles therefore support the identities of a subset of models; they do not authenticate the full panel. We provide the external-reference status per model so that users can prioritise confirmatory testing.

## Usage Notes

The master per-model table, `supplement_per_line.csv`, is the recommended entry point for model selection. It combines histotype, contributing centre, assay coverage, patient-family membership, quality metrics, candidate driver calls, copy-number burden and external-reference status. A standalone browser, `ovcan_viewer_standalone.html`, returns RNA and protein abundance for a supplied gene symbol across the panel. [TO ADD: persistent repository and browser URL.]

Several properties determine how the data should be reused. Population frequencies should use one representative per patient rather than counting related sublines independently. Protein analyses should retain plex membership and feature-detection information, particularly for low-abundance proteins. Variant calls are tumour-only candidates, and copy number is relative and not allele-specific. Histotype comparisons in the expression data remain partly confounded with contributing centre. Passage records are not matched across assays, with differences of up to 20 passages in the 13 models for which both records were available.

Users should confirm STR identity and mycoplasma status before initiating experiments with a model, especially for the 12 models without a Cellosaurus reference profile. Model-selection decisions involving a particular alteration or protein should also be checked against the corresponding per-feature quality fields rather than inferred from the histotype label alone.

## Data Availability

RNA-sequencing reads and processed expression data are available at [GEO/SRA ACCESSION]. Raw proteomic data, search output and processed protein matrices are available through ProteomeXchange at [PRIDE ACCESSION]. Whole-exome reads and processed variant and copy-number records are available at [SRA ACCESSION]. Processed matrices, per-model metadata, data dictionaries, checksums and the standalone browser are available at [REPOSITORY DOI]. Anonymous reviewer-access links are [TO ADD].

## Code Availability

Analysis code is available at [GITHUB URL, RELEASE OR COMMIT] and archived at [ZENODO DOI]. The numbered R scripts are orchestrated by `scripts/run_all.sh`; package versions are recorded in `renv.lock` and `output/package_versions.csv`.

## Author Contributions

[TO ADD: CRediT contribution statement for every author, distinguishing model derivation and culture, generation of each assay layer, computational analysis, data curation, supervision, funding acquisition and manuscript preparation.]

## Competing Interests

[TO ADD: competing-interest declaration for every author, or state that the authors declare no competing interests.]

## Acknowledgements

[TO ADD: acknowledgements that do not meet authorship criteria.]

## Funding

[TO ADD: funding agencies and grant numbers supporting model derivation, each assay layer, data analysis and publication.]

## References

1. Cancer Genome Atlas Research Network. Integrated genomic analyses of ovarian carcinoma. *Nature* **474**, 609-615 (2011). https://doi.org/10.1038/nature10166
2. Domcke, S., Sinha, R., Levine, D. A., Sander, C. & Schultz, N. Evaluating cell lines as tumour models by comparison of genomic profiles. *Nature Communications* **4**, 2126 (2013). https://doi.org/10.1038/ncomms3126
3. Beaufort, C. M. *et al.* Ovarian cancer cell line panel clinical importance of in vitro morphological subtypes. *PLoS ONE* **9**, e103988 (2014). https://doi.org/10.1371/journal.pone.0103988
4. Sauriol, S. A. *et al.* Modeling the diversity of epithelial ovarian cancer through ten novel well characterized cell lines covering multiple subtypes of the disease. *Cancers* **12**, 2222 (2020). https://doi.org/10.3390/cancers12082222
5. Coscia, F. *et al.* Integrative proteomic profiling of ovarian cancer cell lines reveals precursor cell associated proteins and functional status. *Nature Communications* **7**, 12645 (2016). https://doi.org/10.1038/ncomms12645
6. Shrestha, R. *et al.* Multiomics characterization of low-grade serous ovarian carcinoma identifies potential biomarkers of MEK inhibitor sensitivity and therapeutic vulnerability. *Cancer Research* **81**, 1681-1694 (2021). https://doi.org/10.1158/0008-5472.CAN-20-2222
7. Karnezis, A. N. *et al.* Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines. *Gynecologic Oncology* **160**, 568-578 (2021). https://doi.org/10.1016/j.ygyno.2020.12.004
8. Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology* **34**, 525-527 (2016). https://doi.org/10.1038/nbt.3519
9. Soneson, C., Love, M. I. & Robinson, M. D. Differential analyses for RNA-seq transcript-level estimates improve gene-level inferences. *F1000Research* **4**, 1521 (2015). https://doi.org/10.12688/f1000research.7563.1
10. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
11. Garcia, M. *et al.* Sarek a portable workflow for whole-genome sequencing analysis of germline and somatic variants. *F1000Research* **9**, 63 (2020). https://doi.org/10.12688/f1000research.16665.2
12. Benjamin, D. *et al.* Calling somatic SNVs and indels with Mutect2. *bioRxiv* 861054 (2019). https://doi.org/10.1101/861054
13. Talevich, E., Shain, A. H., Botton, T. & Bastian, B. C. CNVkit genome-wide copy number detection and visualization from targeted DNA sequencing. *PLoS Computational Biology* **12**, e1004873 (2016). https://doi.org/10.1371/journal.pcbi.1004873
14. Ghandi, M. *et al.* Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* **569**, 503-508 (2019). https://doi.org/10.1038/s41586-019-1186-3
15. Bairoch, A. The Cellosaurus, a cell-line knowledge resource. *Journal of Biomolecular Techniques* **29**, 25-38 (2018). https://doi.org/10.7171/jbt.18-2902-002
16. Nusinow, D. P. *et al.* Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* **180**, 387-402.e16 (2020). https://doi.org/10.1016/j.cell.2019.12.023
17. Upadhya, S. R. & Ryan, C. J. Experimental reproducibility limits the correlation between mRNA and protein abundances in tumour proteomic profiles. *Cell Reports Methods* **2**, 100288 (2022). https://doi.org/10.1016/j.crmeth.2022.100288

## Tables

### Table 1 | Composition of the resource

| Histotype | Models | Patients | RNA | Protein | Copy number | Mutations | All three |
|---|---:|---:|---:|---:|---:|---:|---:|
| High-grade serous | 24 | 16 | 15 | 15 | 18 | 17 | 9 |
| Clear cell | 8 | 8 | 7 | 7 | 2 | 2 | 2 |
| Mucinous | 3 | 3 | 3 | 3 | 1 | 1 | 1 |
| Endometrioid | 2 | 2 | 2 | 2 | 1 | 1 | 1 |
| Carcinosarcoma | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| Small cell carcinoma of the ovary, hypercalcaemic type | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| Low-grade serous | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| Total | 42 | 34 | 31 | 31 | 23 | 22 | 13 |

### Table 2 | Principal data records

| Record | Contents | Unit | Repository |
|---|---|---|---|
| RNA sequencing | Raw reads; gene-level TPM and estimated counts; variance-stabilised matrix; per-model quality metrics | 31 models, 28 patients | GEO and SRA [ACCESSION] |
| Tandem mass tag proteomics | Raw spectra and search output; relative-abundance matrix; peptide and feature annotations; plex-presence metrics | 31 models, 28 patients | PRIDE [ACCESSION] |
| Whole-exome sequencing | Raw reads; filtered candidate variants; copy-number segments and summary profiles | 23 models for copy number, 22 for mutations, 16 patients | SRA [ACCESSION] |
| Metadata and processed analyses | Master per-model table; sample and family metadata; transcript map; data dictionaries; tables underlying figures | 42 models, 34 patients | [REPOSITORY DOI] |
| Gene-level browser | Self-contained browser for transcript and protein abundance by gene symbol | 31 expression-profiled models | [REPOSITORY OR URL] |

## Figure legends

### Fig. 1 | Resource composition and data generation

(a) Overview of data generation and processing. (b) Assay coverage for 42 models, ordered by histotype and patient family. Filled cells indicate available RNA, protein, mutation or copy-number data; 13 models have all three principal layers. (c) Numbers of models and independent patients by histotype.

### Fig. 2 | Quality assessment of the RNA and protein records

(a) Genes detected and pseudoalignment rate across 31 RNA-sequenced models. Symbols denote contributing centre and open symbols show centre medians. (b) Genes detected in relation to processed fragments. (c) Cross-model transcript and protein spread for 7,896 paired genes. (d) Protein features detected across zero to five tandem mass tag plexes. (e) Vendor-reported and bridge-derived coefficient of variation by protein-abundance decile. (f) Bland-Altman agreement for four adjacent-plex bridge samples.

### Fig. 3 | Expression profiles recover histotype-associated structure

(a) RNA principal components for 31 models, coloured by histotype. (b) The same models coloured by contributing centre, with the principal component 1 commonality decomposition. (c) Distribution of per-gene transcript-protein Spearman correlations across 30 dual-layer models. (d) Per-feature variance decomposition for RNA and protein. (e) Canonical marker abundance summarised by histotype. (f) Marker abundance across individual models. The marker-recovery permutation uses the 28 patient representatives; principal component displays use all 31 models.

### Fig. 4 | Genomic filtering, copy-number profiles and external expression matching

(a) Variant counts before and after Mutect2 and population-frequency filtering. (b) Expression correlation of five models with 67 DepMap ovarian models; each namesake self-match ranks first. (c) Candidate alterations in selected ovarian cancer genes, with evidence tier shown for each call. (d) Arm-level copy-number gain and loss frequencies calculated over 11 high-grade serous patients.
