# Multiomic data from 42 human ovarian cancer cell models

[TODO A01: Add the agreed author list, affiliations, and corresponding author details.]

## Abstract

Ovarian cancer cell models are used to study diseases with distinct histological and molecular features, but molecular data for individual models are often distributed across studies. Here we describe transcriptomic, proteomic, and whole-exome data for 42 human ovarian cancer cell models derived from 34 patients. We generated RNA-sequencing and tandem mass tag proteomic profiles for 31 models each and reprocessed existing whole-exome data to provide candidate coding variants and copy-number profiles for 23 models. Thirteen models have all three data types. The resource includes abundance matrices, variant annotations, copy-number segments, and metadata identifying models derived from the same patient. Technical validation evaluated sequencing quality, protein detection and repeatability, established ovarian cancer markers, agreement between RNA and protein profiles, and correspondence with external expression references. These data support selection of cell models with documented histotype, genotype, transcript abundance, and protein abundance for ovarian cancer research.

## Background & Summary

Ovarian cancer includes histotypes with different molecular characteristics. High-grade serous carcinoma commonly carries *TP53* alterations and extensive copy-number changes (Cancer Genome Atlas Research Network, 2011; DOI: 10.1038/nature10166). The choice of a cell model therefore influences which features of the disease can be studied. Comparisons with tumour profiles have shown that several widely used ovarian cancer cell lines poorly represent high-grade serous carcinoma (Domcke et al., 2013; DOI: 10.1038/ncomms3126). Systematic characterisation of cell-line panels has also documented variation in molecular features, morphology, and historical annotations (Beaufort et al., 2014; DOI: 10.1371/journal.pone.0103988).

Previous resources include characterisation of patient-derived ovarian cancer cell lines (Sauriol et al., 2020; DOI: 10.3390/cancers12082222), proteomic profiles across ovarian cancer histotypes (Coscia et al., 2016; DOI: 10.1038/ncomms12645), and matched genomic, transcriptomic, and proteomic data for low-grade serous carcinoma (Shrestha et al., 2021; DOI: 10.1158/0008-5472.CAN-20-2222). For many models, particularly those representing rare histotypes, RNA, protein, and genomic data are available from separate studies with different experimental conditions. Linked profiles and sample metadata facilitate comparisons while making these differences identifiable.

The present resource combines newly generated RNA-sequencing and proteomic data with reprocessed whole-exome data from a collection of 42 ovarian cancer cell models contributed by three Canadian centres (Fig. 1). The models represent 34 patients and seven recorded histotype groups. Thirteen models are sublines from five patients; the metadata identifies their common patient of origin and indicates which model was selected when an analysis required one model per patient. Thirty models have both RNA and protein profiles, and 13 have RNA, protein, and exome profiles. The data include quantitative quality metrics and annotations for candidate coding alterations, relative copy number, and protein detection. These records allow researchers to select models on the basis of measured molecular features and to account for differences in assay coverage and patient origin.

[TODO A02: Confirm and cite any previous publication that used these exact RNA, protein, or exome datasets, and specify which data are first released with this descriptor.]

## Methods

### Cell models and study design

The collection contains 42 models: 29 contributed by the Centre de recherche du Centre hospitalier de l'Université de Montréal, 12 by BC Cancer-OVCARE, and one by the Ottawa Hospital Research Institute (Table 1). The accompanying model metadata contains one record for each included model. Histotype groups are high-grade serous, low-grade serous, clear cell, endometrioid, mucinous, carcinosarcoma, and small cell carcinoma of the ovary, hypercalcaemic type. Historical annotations are retained alongside documented revisions. In particular, TOV112D remains identifiable within the historical endometrioid group but was subsequently reassigned to dedifferentiated ovarian carcinoma with SWI/SNF loss (Karnezis et al., 2021; DOI: 10.1016/j.ygyno.2020.12.004).

Related sublines were assigned a common patient identifier. For expression analyses requiring independent patient observations, we selected one model per patient using one point each for RNA, protein, and candidate-variant availability and 0.1 points for copy-number availability, with alphabetical model-name order breaking ties. This gave 34 representatives overall, 28 with RNA profiles, and 27 with both RNA and protein profiles. Descriptive plots retain all models. Mutation and copy-number summaries identify their patient-level aggregation rule separately below.

[TODO M01: Confirm culture and harvest conditions for each centre and model. Published CRCHUM protocols describe complete OSE medium with 10% fetal bovine serum, 0.5 µg/mL amphotericin B, and 50 µg/mL gentamicin, at 37 °C and 5% CO2. Reported oxygen concentrations differ: 5% for the 1369, 2295, and 3133 series, and 7% in the later derivation studies. Confirm the conditions used for the present samples, including passage, before replacing this TODO with the final methods. Source-specific draft text and derivation references are supplied in the author confirmation document.]

### RNA sequencing and quantification

[TODO M02: Add the RNA extraction protocol and RNA integrity assessment.]

Libraries were prepared from 250 ng of total RNA. Polyadenylated RNA was enriched with the NEBNext Poly(A) Magnetic Isolation Module. First-strand cDNA synthesis used the NEBNext RNA First Strand Synthesis Module, followed by the NEBNext Ultra Directional RNA Second Strand Synthesis Module. Library preparation was completed with the NEBNext Ultra II DNA Library Prep Kit for Illumina, using adaptors and PCR primers from New England BioLabs. Libraries were normalised, pooled, denatured in 0.02 N NaOH, and neutralised with HT1 buffer. The pool was loaded at 200 pM onto an Illumina NovaSeq S4 lane using the Xp protocol, with 1% PhiX control. Sequencing comprised 2 × 100 cycles in paired-end mode in one run. Base calling used RTA v3, and bcl2fastq2 v2.20 generated demultiplexed FASTQ files.

Transcript abundances were estimated with kallisto 0.46.0 (Bray et al., 2016; DOI: 10.1038/nbt.3519). The reference contained 185,299 transcript targets matching the Ensembl release 93 human cDNA reference in both versioned identifier and sequence length. We derived the transcript-to-gene map from this reference and summarised the estimates with tximport 1.36.1 (Soneson et al., 2015; DOI: 10.12688/f1000research.7563.1). All transcript targets were represented in the map. The reference identity, mapping table, and file checksums accompany the processed data.

Gene-level transcripts per million (TPM) were retained for 39,733 genes. For count-based analyses, genes with at least 10 estimated counts in at least two models were retained, yielding 22,678 genes. Variance-stabilised expression values were calculated with DESeq2 1.48.2 (Love et al., 2014; DOI: 10.1186/s13059-014-0550-8). Gene detection was defined as at least one assigned count among the genes retained by the count filter. Quality metrics distinguish fragments processed by kallisto from the assigned count total after gene summarisation and filtering.

### Tandem mass tag proteomics

Protein abundance was measured in five 11-plex tandem mass tag sets. Each set contained a pooled internal standard. Four samples were measured in adjacent sets, connecting the five plexes. Each connection used a different sample, so repeatability was evaluated across four paired samples. VOA3993 was an additional technical bridge sample and is not part of the 42-model biological collection.

[TODO M03: Confirm the proteomics acquisition protocol with the generating laboratory. A draft based on the 11-plex cell-pellet protocol reported by Orlando et al. with Negri and Morin (2020; DOI: 10.7554/eLife.59073) is supplied in the author confirmation document. Confirm sample preparation, digestion, fractionation, instrument and acquisition mode, search engine and database version, peptide and protein false-discovery thresholds, reporter-ion interference filtering, protein summarisation, channel assignment, pooled-standard composition, prior normalisation, and the definition of the supplied coefficient of variation.]

The source protein table contained 8,430 rows. Three rows without a gene symbol were removed, leaving 8,427 features. Where rows shared a gene symbol, a representative was selected using measurement completeness, peptide support, and q-value. The remaining rows were retained with identifiers combining the gene symbol and UniProt accession. We retained the supplied log2 protein-abundance values, normalised using the pooled internal standard. [TODO M03: Confirm the upstream normalisation formula and scaling.] Missing values were retained, and detection was recorded by model and plex.

For each paired bridge sample, we calculated the difference in log2 abundance between adjacent plexes. Under equal, independent measurement error, the standard deviation of the differences divided by the square root of two estimates the standard deviation of an individual measurement. We converted this estimate to a lognormal coefficient of variation and summarised it by mean-abundance decile. This calculation estimates technical repeatability of the bridge measurements.

### Whole exome sequencing and variant processing

Existing whole-exome data were processed against GRCh38 with nf-core/Sarek (Garcia et al., 2020; DOI: 10.12688/f1000research.16665.2). Small variants were called with GATK Mutect2 and filtered with FilterMutectCalls, both version 4.5.0.0, in tumour-only mode without matched normal samples (Benjamin et al., 2019; DOI: 10.1101/861054).

Variants were annotated with Ensembl Variant Effect Predictor 113.0 against GRCh38.p14 (McLaren et al., 2016; DOI: 10.1186/s13059-016-0974-4). Archived headers identify 1000 Genomes phase 3, gnomAD exome and genome 4.1, ClinVar 202404, and dbSNP 156 resources.

[TODO M04: Retrieve the cluster run records for library preparation, sequencing instrument and read length, nf-core/Sarek version and configuration, alignment-reference checksum, original capture design, and run-level sequencing QC.]

Variant filtering first retained records that passed Mutect2 quality filters. It then removed alleles with a maximum annotated population frequency greater than 0.001 across the available gnomAD and 1000 Genomes fields and retained coding and splice alterations in the maftools non-synonymous classes. Calls without an available population frequency were retained and remain identifiable in the annotations. The resulting variants are candidates for follow-up; rare inherited variants can remain after tumour-only filtering.

Where the source protein-change annotation was missing, we reconstructed a readable consequence label from the available amino-acid and position fields. Reconstructed labels are flagged, and incomplete labels should be interpreted with the accompanying genomic coordinates and transcript annotations. For example, H214fs denotes a frameshift beginning at histidine 214. Complex or insertion labels lacking sufficient residue context were flagged as non-canonical.

We prioritised candidates in 19 selected ovarian cancer genes using three tiers. Tier 1 prioritised established activating hotspots, *TP53* missense alterations supported by DNA-binding-region location or ClinVar annotation, and truncating candidates in selected tumour suppressors. Tier 2 contained other candidates supported by gene and variant class. Tier 3 contained uncertain calls, including the *BRCA1* and *BRCA2* candidates. Tier 1 and Tier 2 calls were used for candidate-alteration frequency summaries; Tier 3 calls were retained for inspection. These tiers organise evidence for functional relevance, and do not establish whether a variant is somatic. For patient-level alteration frequencies, a patient was counted once when a qualifying call was present in any of their sequenced sublines.

Copy-number profiles were estimated with CNVkit 0.9.10 using five public healthy exomes (SRR4039087, SRR4039088, SRR4039089, SRR4039096, and SRR4039097; BioProject PRJNA339046) as a pooled reference (Talevich et al., 2016; DOI: 10.1371/journal.pcbi.1004873). Each profile was centred on its autosomal probe-weighted median. Segments were intersected with chromosome-arm boundaries defined by GRCh38 cytobands, excluding centromeric intervals. An arm was called gained or lost when more than half its annotated non-centromeric length had a centred log2 copy ratio above 0.20 or below −0.20, respectively. The fraction of genome altered was the fraction of represented autosomal segment length above either threshold. Acrocentric short arms were excluded. Patient-level arm frequencies counted each patient once if the event occurred in any of their profiled models. For patient-level copy-number burden, per-model fractions were averaged within each patient before summarising across patients. These measurements describe relative total copy number; allele-specific copy number, loss of heterozygosity, and homologous-recombination-deficiency scores were not estimated.

The archived CNVkit target-bin definition contains 290,475 bins and is included with the processing records. [TODO M05: Add the shared capture-kit name and design version, and the original vendor target BED. The authors report that reference and cell-model exomes used the same kit; the recovered CNVkit bins describe the analysis footprint rather than the original vendor design.]

### Technical validation analyses

Principal component analysis summarised transcriptional patterns across models. The variance-stabilised values of the 2,000 most variable genes were centred, without scaling, and models were treated as observations. We fitted principal-component scores to histotype, centre, or both and partitioned the unadjusted R-squared into unique histotype, unique centre, and shared components. Each unique component was the increase in explained variance when that factor was added to a model containing the other factor; the shared component was the sum of the separate R-squared values minus the joint R-squared. The analysis was repeated after selecting one model per patient. Because histotype and contributing centre were partly confounded, these analyses assess their joint association with expression patterns.

Established histotype markers were evaluated in the 28 RNA-profiled patient representatives. For an expected high marker, recovery required the intended group to rank among the two highest group means and to exceed 1 log2(TPM + 1). Expected loss markers required the intended group to rank among the two lowest means. We compared the total recovered count with 20,000 joint permutations of histotype labels within contributing centre, preserving group sizes, the centre composition, and correlations among markers. An unrestricted label permutation was retained as a sensitivity analysis. This is a panel-level check of expected expression patterns.

RNA–protein agreement was assessed in the 30 models with both data types. Gene-level TPM values sharing a gene symbol were summed, then transformed as log2(TPM + 1), and matched to the representative protein feature for that symbol. We calculated the Spearman correlation across models for each gene quantified in at least 10 paired models. The analysis was repeated in the 27 patient representatives with both profiles. Asymptotic two-sided P values were adjusted separately in each analysis by the Benjamini–Hochberg method. Patient-representative tests were used for inference; tests retaining related sublines were descriptive. We also described within-model correlations across genes and the spread of each gene across models using the interquartile range and standard deviation. The supplied protein values retain feature-dependent abundance baselines; their upstream normalisation and scaling affect interpretation of correlations across genes within a model.

Five models were compared with 67 DepMap 24Q4 Public version 1 ovarian models over 2,000 variable genes selected from the shared expression matrix (Broad DepMap, 2024; DOI: 10.25452/figshare.plus.27993248.v1). Cellosaurus name queries identified external short tandem repeat profiles in the documented 23 July 2026 snapshot (Bairoch, 2018; DOI: 10.7171/jbt.18-2902-002). Input checksums, cached responses, and entry versions accompany the records.

[TODO M06: Confirm the Cellosaurus database release or original retrieval documentation, if available.]

Analyses used R 4.5.2 and Bioconductor 3.21, with a project seed of 1234. The versioned workflow records the script order and package environment and produces the tables underlying each figure.

### Ethics

This study profiled pre-existing cell models. [TODO M07: Confirm the ethical approvals and consent framework for the originating patient samples, the institutional determination for this study, and whether the shared genomic data require controlled access. Complete the journal Human Data Checklist using the applicable institutional records.]

### Use of generative artificial intelligence

Coding assistants were used for code review, analysis development, and manuscript editing. Numerical summaries were produced by the analysis workflows. [TODO M08: Confirm the tools and versions used and the final author verification statement before submission.]

## Data Records

The resource links each assay record to a common cell-model identifier. The model metadata lists the 42 included models, their recorded histotype, contributing centre, assay availability, and patient of origin. Separate fields identify related sublines and representatives used for patient-level summaries. Table 2 describes the principal records and their measurement units. The release manifest identifies the corresponding files, formats, and checksums, and the data dictionary defines their fields and missing-value conventions.

RNA records comprise gene-level TPM, filtered estimated counts, variance-stabilised expression, transcript annotations, and per-model quality metrics. Protein records contain 8,427 relative-abundance features with gene and UniProt annotations, peptide support, and detection by plex. Exome records contain filtered candidate variants, copy-number segments, and summary profiles for 23 models. Candidate-variant records include quality filters, population annotations, consequence-label status, and the evidence tier where applicable.

Protein missingness is strongly structured by tandem mass tag plex. A missing entry indicates that a feature was not quantified in the relevant plex, rather than demonstrating absence of the protein in that model. The detection annotations allow users to select features with adequate coverage for a particular comparison. Similarly, patient identifiers allow related sublines to be accounted for when constructing independent observations.

[TODO D01: Add repository accession numbers, persistent record links, licences, and anonymous reviewer access. Deposit raw reads in an appropriate sequence repository, raw spectra and search results through ProteomeXchange/PRIDE, and processed records with versioned metadata and dictionaries in a persistent repository. Resolve any human genomic access requirements before deposition.]

## Technical Validation

### Sequencing quality and protein repeatability

Across the 31 RNA libraries, median pseudoalignment was 91.1%, median depth was 64.3 million paired-end fragments, and a median of 20,270 genes was detected (Fig. 2a,b). Gene detection differed by contributing centre, which should be considered in analyses sensitive to sequencing depth or culture conditions.

The protein matrix contains 8,427 features. The supplied peptide table contains 146,830 quantified peptide records. Seventy source features were not quantified in any included model and are identified by their detection annotations. Models within a plex shared the same feature count, indicating that detection differences arose mainly between plexes (Fig. 2d). Mean log2 differences across the four bridge pairs ranged from −0.028 to 0.009, and bridge-derived per-measurement coefficients of variation ranged from 10.4% to 13.2% (Fig. 2f).

Repeatability depended on protein abundance. The bridge-derived coefficient of variation declined from 21.3% in the lowest abundance decile to 5.2% in the highest (Fig. 2e). The bridge results provide an empirical estimate of agreement between adjacent plexes and support the use of abundance and detection information when comparing individual proteins. These precision estimates apply to the technical bridge comparisons.

### Expected expression patterns and agreement between assays

RNA principal components showed separation among recorded histotype groups (Fig. 3a). For the first component, which explained 20.6% of expression variance, histotype accounted for a unique 42.4% of the variance in component scores, centre for 0.2%, and their shared component for 31.0%. After selecting one model per patient and refitting the variable-gene selection and principal components, the corresponding fractions were 52.0%, 0.3%, and 26.9% (Fig. 3b). All high-grade serous RNA models originated from one centre, so effects of culture or processing aligned with histotype cannot be separated completely from histotype-associated expression.

Among 25 prespecified marker checks, 15 were recovered, compared with a within-centre permutation mean of 8.46 (permutation P = 0.0101) (Fig. 3d,e). This aggregate assessment supports the presence of established expression patterns in the panel. Histotype groups with few independent patients remain limited for inference about their broader populations, and the TOV112D annotation requires particular care in comparisons involving the historical endometrioid group.

Across 8,033 genes, the median correlation between RNA and protein abundance across models was 0.400 (Fig. 3c). After selecting one model per patient, the median was 0.418 across 7,969 genes in 27 paired models. Gene-specific differences may reflect regulation of translation or protein turnover, assay error, and variation in culture or passage between assays. The narrower spread of the supplied protein abundances compared with RNA abundance (Fig. 2c) may also reflect isobaric ratio compression. These descriptive comparisons do not quantify the contribution of any individual mechanism.

### Genomic quality checks and established alterations

The unfiltered tumour-only calls included variants that did not pass caller quality filters and common population variants. Sequential filtering reduced 582,474 variant records to 16,081 records passing Mutect2 filters and 6,194 population-filtered coding or splice candidates across 23 models (Fig. 4a). The filtering steps remove different classes of unsupported calls and are retained in the data record so that users can inspect the evidence for an individual candidate.

A qualifying *TP53* alteration was present in all 18 high-grade serous models with exome data, representing 11 patients (Fig. 4c). The high-grade serous copy-number profiles showed recurrent arm-level changes reported in ovarian tumours (Cancer Genome Atlas Research Network, 2011; DOI: 10.1038/nature10166). Among 11 patients, gains of 20q and 3q were present in 10 and nine patients, respectively, and 17p loss was present in nine (Fig. 4d). The median fraction of represented autosomal genome altered was 0.622 after averaging subline values within each high-grade serous patient. Other models showed heterogeneous copy-number burdens, including marked alteration in the clear cell model TOV3392D. The small non-serous groups do not support a general ranking of copy-number burden across histotypes.

### Comparison with external model references

Each of the five models with a DepMap namesake had its highest expression correlation with that namesake among 67 ovarian models, and all five were reciprocal best matches (Fig. 4b). The smallest margin over another model involved BIN67 and COV434, both annotated as small cell carcinoma of the ovary, hypercalcaemic type. These expression comparisons support correspondence with the named reference models while allowing for similarity among models of the same histotype.

TOV112D carried a truncating *SMARCA4* candidate alongside low SMARCA4 protein, low *SMARCA2* RNA, and *TP53* R175H, consistent with its published reassignment to dedifferentiated ovarian carcinoma (Karnezis et al., 2021; DOI: 10.1016/j.ygyno.2020.12.004). Such established alterations provide context for model annotations; additional functional conclusions require targeted validation.

Cellosaurus contained external short tandem repeat reference profiles for 30 of the 42 model names. The availability of a reference profile does not establish that the corresponding study stock was matched to it. [TODO M09: Obtain the contributing laboratories' short tandem repeat and mycoplasma test records, including dates, methods, results, and their relationship to the harvested stocks. Replace this TODO with the documented authentication and contamination-control procedures.]

## Usage Notes

Model metadata provides context on histotype, centre, patient origin, assay coverage, quality metrics, candidate alterations, and relative copy-number burden. A gene browser at https://www.cooklab.ca/ovcan_viewer supports inspection of transcript and protein abundance by gene symbol. A versioned copy accompanies the processed data. [TODO D02: Verify that the hosted browser serves the final release and add its archived version identifier.]

Analyses across independent patients should account for related sublines. The supplied representative annotations provide one reproducible selection; questions about disease progression or variation among a patient's sublines may require retaining the related models with an appropriate paired or hierarchical analysis. Protein comparisons should retain plex and detection annotations, especially for low-abundance features. Missing protein entries should not be imputed under an unexamined assumption that values are missing at random.

Variant calls should be interpreted with their tumour-only origin and annotation evidence. Relative copy-number segments describe gains and losses on the scale of the centred profile, and should not be treated as allele-specific or absolute copy number. RNA and protein comparisons also require attention to the partial association between histotype and centre. Passage records differed by up to 20 passages in the 13 models with both RNA and exome passage information; assay profiles therefore need not represent the same culture harvest.

Before follow-up experiments, users should consult current authentication, mycoplasma, and culture records for the stock they obtain. Decisions involving a particular gene or protein should be supported by its assay coverage and quality annotations as well as the model's histotype label.

## Data Availability

RNA reads and processed RNA data: [TODO D01: GEO/SRA or appropriate sequence-repository accession]. Proteomic spectra, search results, and processed protein data: [TODO D01: ProteomeXchange/PRIDE accession]. Whole-exome reads, candidate variants, and copy-number records: [TODO D01: sequence-repository accession and any access conditions]. Processed matrices, model metadata, data dictionaries, checksums, and the archived browser: [TODO D01: repository DOI]. Anonymous reviewer access: [TODO D01: links and access instructions].

## Code Availability

The analysis workflow, package environment, and scripts generating the figures and release records are available at [TODO D03: public code URL and release or commit identifier], with a persistent archive at [TODO D03: archive DOI].

## Author Contributions

[TODO A03: Add contributions for every author, including model derivation and culture, data generation, analysis, curation, supervision, funding, and manuscript preparation.]

## Competing Interests

[TODO A04: Add the agreed competing-interest declaration.]

## Acknowledgements

[TODO A05: Add acknowledgements.]

## Funding

[TODO A06: Add funders and grant numbers.]

## References

1. Cancer Genome Atlas Research Network. Integrated genomic analyses of ovarian carcinoma. *Nature* **474**, 609-615 (2011). https://doi.org/10.1038/nature10166
2. Domcke, S., Sinha, R., Levine, D. A., Sander, C. & Schultz, N. Evaluating cell lines as tumour models by comparison of genomic profiles. *Nature Communications* **4**, 2126 (2013). https://doi.org/10.1038/ncomms3126
3. Beaufort, C. M. *et al.* Ovarian cancer cell line panel (OCCP): clinical importance of in vitro morphological subtypes. *PLoS ONE* **9**, e103988 (2014). https://doi.org/10.1371/journal.pone.0103988
4. Sauriol, S. A. *et al.* Modeling the diversity of epithelial ovarian cancer through ten novel well characterized cell lines covering multiple subtypes of the disease. *Cancers* **12**, 2222 (2020). https://doi.org/10.3390/cancers12082222
5. Coscia, F. *et al.* Integrative proteomic profiling of ovarian cancer cell lines reveals precursor cell associated proteins and functional status. *Nature Communications* **7**, 12645 (2016). https://doi.org/10.1038/ncomms12645
6. Shrestha, R. *et al.* Multiomics characterization of low-grade serous ovarian carcinoma identifies potential biomarkers of MEK inhibitor sensitivity and therapeutic vulnerability. *Cancer Research* **81**, 1681-1694 (2021). https://doi.org/10.1158/0008-5472.CAN-20-2222
7. Karnezis, A. N. *et al.* Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines. *Gynecologic Oncology* **160**, 568-578 (2021). https://doi.org/10.1016/j.ygyno.2020.12.004
8. Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology* **34**, 525-527 (2016). https://doi.org/10.1038/nbt.3519
9. Soneson, C., Love, M. I. & Robinson, M. D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* **4**, 1521 (2015). https://doi.org/10.12688/f1000research.7563.1
10. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
11. Garcia, M. *et al.* Sarek: a portable workflow for whole-genome sequencing analysis of germline and somatic variants. *F1000Research* **9**, 63 (2020). https://doi.org/10.12688/f1000research.16665.2
12. Benjamin, D. *et al.* Calling somatic SNVs and indels with Mutect2. *bioRxiv* 861054 (2019). https://doi.org/10.1101/861054
13. Talevich, E., Shain, A. H., Botton, T. & Bastian, B. C. CNVkit: genome-wide copy number detection and visualization from targeted DNA sequencing. *PLoS Computational Biology* **12**, e1004873 (2016). https://doi.org/10.1371/journal.pcbi.1004873
14. Broad DepMap. DepMap 24Q4 Public (version 1). *figshare* (2024). https://doi.org/10.25452/figshare.plus.27993248.v1
15. Bairoch, A. The Cellosaurus, a cell-line knowledge resource. *Journal of Biomolecular Techniques* **29**, 25-38 (2018). https://doi.org/10.7171/jbt.18-2902-002

16. McLaren, W. *et al.* The Ensembl Variant Effect Predictor. *Genome Biology* **17**, 122 (2016). https://doi.org/10.1186/s13059-016-0974-4

## Tables

### Table 1 Composition of the resource

| Histotype group | Models | Patients | RNA | Protein | Copy number | Variants | All three |
|---|---:|---:|---:|---:|---:|---:|---:|
| High-grade serous | 24 | 16 | 15 | 15 | 18 | 18 | 9 |
| Clear cell | 8 | 8 | 7 | 7 | 2 | 2 | 2 |
| Mucinous | 3 | 3 | 3 | 3 | 1 | 1 | 1 |
| Endometrioid | 2 | 2 | 2 | 2 | 1 | 1 | 1 |
| Carcinosarcoma | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| Small cell carcinoma of the ovary hypercalcaemic type | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| Low-grade serous | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| Total | 42 | 34 | 31 | 31 | 23 | 23 | 13 |

Assay columns count models. All three denotes availability of RNA, protein, and exome data. Thirty models have both RNA and protein profiles. Groups retain the source annotations for traceability; TOV112D is counted in the historical endometrioid group and is annotated separately as reassigned to dedifferentiated ovarian carcinoma.

### Table 2 Principal data records

| Record | Contents | Coverage | Repository |
|---|---|---|---|
| RNA sequencing | Reads; TPM and estimated counts; variance-stabilised expression; transcript map; quality metrics | 31 models from 28 patients | Sequence repository and processed-data archive [TODO D01] |
| Proteomics | Spectra and search results; relative protein abundance; feature and peptide annotations; plex detection | 31 models from 28 patients | PRIDE [TODO D01] |
| Whole exome | Reads; candidate coding variants; copy-number segments and summaries | 23 models from 16 patients | Sequence repository and processed-data archive [TODO D01] |
| Model metadata | Histotype; centre; assay availability; patient identifiers; representatives; quality and external-reference annotations | 42 models from 34 patients | Processed-data archive [TODO D01] |
| Browser and release documentation | Gene browser; manifest; data dictionaries; checksums; figure source tables | Release-wide | Processed-data archive [TODO D01] |

TPM, transcripts per million. Protein values are supplied log2 abundances normalised using the pooled internal standard. Copy-number values are centred log2 ratios. Data dictionaries specify the units and missing-value conventions for each field.

## Figure legends

### Figure 1 Resource composition and data generation

(a) Newly generated RNA and protein data and reprocessed whole-exome data from 42 models representing 34 patients. Counts indicate models with each assay. (b) Model and independent-patient counts by recorded histotype; bars indicate model counts. (c) Assay coverage, with high-grade serous models on the left and other histotypes on the right. Filled slate cells indicate available profiles; outlined cells indicate assays not included for that model. The exome column represents both variant and copy-number records, which cover the same 23 models. Thirteen models have all three assays. Coloured strips identify histotype using the key in b, and brackets label patients represented by multiple models. Historical group labels are retained as described in Table 1. HGS, high-grade serous; CC, clear cell; MC, mucinous; EC, endometrioid; MMMT, carcinosarcoma; SCCOHT, small cell carcinoma of the ovary, hypercalcaemic type; LGS, low-grade serous.

### Figure 2 Quality assessment of the RNA and protein records

(a,b) Detected genes versus pseudoalignment and processed paired-end fragments in 31 RNA libraries. Symbols indicate centre; large hollow symbols indicate centre medians. Detection is an estimated count of at least one among the 22,678 genes retained by the global count filter. (c) Cross-model interquartile ranges for 8,035 genes with at least ten paired measurements. RNA units are log2(TPM + 1); protein units are supplied log2 normalised abundance. Boxes span the interquartile range, central lines indicate medians, and whiskers extend to the most extreme observations within 1.5 interquartile ranges. Numbers give medians. The display ends just above the pooled 99th percentile. (d) Protein-feature counts by number of quantified tandem mass tag plexes, on a logarithmic count axis. The rust point identifies the 70 features without measurements in model channels. (e) Supplied coefficient-of-variation medians and interquartile ranges (slate line and band) and bridge-derived per-measurement estimates (rust) by abundance decile. The bridge estimate assumes equal, independent measurement errors. (f) Primary-minus-bridge log2 differences versus mean abundance for four adjacent-plex links. Each point represents a paired protein feature; facet labels give counts. The four links use different samples; the asterisk identifies the external proteomics-only bridge VOA3993. Solid lines indicate mean differences; dashed lines indicate mean differences ± 1.96 standard deviations. The vertical limits exclude 0.191% of observations from view; all observations contribute to the agreement statistics.

### Figure 3 Expression patterns and agreement between assays

(a) Principal components of the 2,000 most variable variance-stabilised RNA genes across 31 models, centred without scaling. Fill indicates recorded histotype and shape indicates centre. The two historically endometrioid models are labelled. (b) Unadjusted partition of first-component score variance into unique histotype, unique centre, shared and unexplained components, comparing all models with 28 patient representatives. Variable-gene selection and principal components were refitted in the representative subset. (c) Per-gene Spearman correlations between RNA and protein abundance across 30 paired models. The 8,033 genes have at least ten paired measurements and nonzero variance in both assays; the rust line marks the median. (d) Means of marker log2(TPM + 1) across 28 patient representatives, standardised within each marker across the six histotype means. Column labels give patient counts. Solid boxes mark the expected high-expression histotype; dashed boxes and downward arrows mark expected low expression. (e) Individual-model log2(TPM + 1) values in the same marker order. Coloured points mark the intended histotype and pale points mark other models; vertical jitter separates points. The dashed line marks the high-expression recovery threshold, applied to histotype means. Fifteen of 25 marker checks were recovered using patient representatives; the aggregate assessment used 20,000 joint label permutations within centre. MKI67 is a proliferation control excluded from that count. Historical labels and small histotype groups require the qualifications described in the text.

### Figure 4 Genomic filtering and external expression matching

(a) Variant counts before filtering, after caller-quality filtering and after population-frequency and coding-consequence filters, for all 23 models and the example OV2295. (b) Expression correlations with 67 DepMap ovarian models for five models with namesake references. Filled rust points identify the namesake, open points identify the highest-correlating non-namesake, and pale ticks show other comparisons. Every namesake ranks first. (c) Candidate alterations in selected ovarian cancer genes. Shapes identify variant classes and slate shades identify priority tiers. Multiple classes within a gene-model cell use a diamond. Counts above the matrix show coding candidates, histotype strips identify recorded groups, and brackets label shared patients. Right-hand bars count patients with Tier 1 or Tier 2 calls, once per patient if present in any sequenced subline. Empty cells do not establish a wild-type genotype; zero patient counts can occur for genes with Tier 3 calls only. (d) Arm-level gain and loss frequencies among 11 high-grade serous patients, counting an event once per patient if present in any of their models. These tumour-only candidates and relative copy-number profiles retain the qualifications described in Methods.
