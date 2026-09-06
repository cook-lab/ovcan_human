# Proposed Data Descriptor additions for author review

These paragraphs are prepared from the expanded analyses and have not been inserted into manuscript v9. Final table/figure numbering should follow the author's choice of supplementary material. They deliberately describe characterization and technical limitations rather than clinical assignment or treatment prediction.

## Methods: supplementary molecular annotation

We constructed supplementary model-level annotations for selected cancer-associated loci and protein targets. Relative copy-number summaries were calculated for 52 prespecified genes across 23 exome-profiled models from 16 patients using the corrected target-only segmentation. For each locus, we retained overlap-weighted segment estimates, individual gene-overlapping target ratios, zero-depth counts and the fraction of the gene span covered by segments. Separate review flags identified discordance between segments and target bins, incomplete locus coverage and extensive zero-depth measurements. Selected NF1, CDKN2A and AKT2 loci were additionally evaluated against pinned canonical-exon annotations and the five recovered reference-normal coverage profiles. These summaries describe relative dosage and do not determine absolute or allele-specific copy number.

For 19 expression targets, RNA and protein measurements were paired by exact model before patient-level comparisons. We assessed Spearman correlations using one preselected model per patient and sensitivity analyses using within-patient means and histotype-specific subsets. We also compared four combinations of within-cohort RNA and protein ranks to identify model-selection leads while preserving raw measurements, missing-assay status and supplied isoDoping annotations. Rank-based measures describe this collection and do not define clinical expression thresholds.

Selected variant candidates were traced to exact archived VCF alleles to recover read counts, caller allele fractions, strand/orientation summaries and local phase identifiers. Clinical assertions, when available from reviewed public records, were retained separately from research prioritization tiers, tumour-only origin uncertainty and evidence of biallelic inactivation. Nearby locally phased variants were flagged for combined-haplotype interpretation.

## Technical validation: interpretation of exploratory annotations

The extended locus analysis identified concordant AKT2 DNA/RNA/protein evidence in OV3331 and highlighted NF1 and CDKN2A regions with profound or partial target depletion. These examples illustrate why segment summaries should be interpreted alongside underlying targets: a broad segment can obscure a localized peak or partial deletion, while a median computed only from positive-depth bins can conceal dropout. Related sublines can also differ substantially, as illustrated by FOLR1 expression within the 2295 family. Model-specific results should therefore retain sample identity and passage context rather than being generalized to all models from a patient.

SBS3 estimates were evaluated across genome- versus capture-opportunity-adjusted reference profiles, full and restricted signature dictionaries, and baseline versus stricter population/read-support screens. Across all 23 models, the target-adjusted full-dictionary 95% bootstrap ranges included zero under both baseline and stricter screens. The findings therefore do not support robust SBS3-positive model assignments from these tumour-only exomes. Genomic scar scoring would require suitable allele-specific copy-number inputs and fit validation; neither relative copy-number burden nor a signature point estimate establishes homologous-recombination deficiency.

## Usage notes

The supplementary annotations support the selection of models for targeted validation. For example, concordant expression of a putative surface target can guide subsequent membrane-protein assessment, while low bulk abundance or an inherited disease-associated variant cannot establish protein loss, present pathway deficiency or therapeutic response. MSI, allele-specific genomic scars and absolute focal amplification remain distinct measurements requiring additional analysis or orthogonal evidence. Updated model annotations should preserve these evidence types separately.

## Source / insertion notes

- Methods, exact denominators and external figure legends: module READMEs in this report directory.
- The 52-gene screen includes analyst-defined review thresholds, not validated amplification/deletion cutoffs; retain that distinction if a binary screen table is included.
- The 19-target matrix includes PGR and DPEP3 with no quantified protein; they must remain unscored for joint analyses.
- Exact ClinVar assertions are supported only for the specifically curated records, with germline-disease, somatic oncogenicity and functional evidence distinguished. Cite the variant report's primary sources if examples are retained.
- The gene-CN and expression tables are suitable supplementary characterization candidates. The SBS3 diagnostic belongs in technical validation/method sensitivity, not a headline HRD classification figure.
