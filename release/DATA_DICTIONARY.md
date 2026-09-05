# Data dictionary

CSV files use a header row, UTF-8 text, and a decimal point. Empty cells or `NA` mean missing or not applicable, not a biological zero. Boolean values are `TRUE`/`FALSE` unless a source field is explicitly retained. Percentages span 0–100; proportions span 0–1. Every matrix has one identifier column followed by model columns.

## Cohort metadata

`metadata/models.csv` contains the complete 42-model cohort. `cell_line` is the join key; `histotype_code` preserves the analysis grouping, and `histotype_label` and `histotype_note` explain its meaning and known reassignments. Codes are HGS (high-grade serous), LGS (low-grade serous), CC (clear cell), EC (historical endometrioid grouping), MC (mucinous), MMMT (carcinosarcoma), and SCCOHT (small cell carcinoma of the ovary, hypercalcaemic type). These codes describe the supplied models, not a new histological diagnosis.

`contributing_centre` is the supplying centre. `patient_id` groups models from the same patient; `models_from_same_patient` counts that patient's models. `selected_patient_model` chooses one model per patient by assay coverage with a lexical tie-break. Assay-specific analyses can require another selection restricted to the available assay intersection.

`rna_available`, `protein_available`, `copy_number_available`, and `variants_available` indicate that the corresponding processed record is present. `rna_sample_id` retains the library identifier; RNA and WES passage values are assay-specific source annotations. `tmt_plex` and `tmt_channel` identify the protein measurement channel. The two `stock_*_documentation` fields distinguish pending laboratory confirmation from negative test results.

The same table includes per-model RNA and protein detection, sequencing depth, retained coding-candidate count, and autosomal FGA. Fields ending `_M` are in millions; `rna_assigned_gene_counts_M` is calculated after the global expression filter. Cellosaurus accession, reference STR availability, and external problematic-line annotations describe published reference records and do not authenticate the analysed stocks.

## RNA records

`rna_tpm.csv` contains unfiltered gene-level transcripts per million, summed across the transcripts assigned to each Ensembl `gene_id`. `rna_counts.csv` contains tximport estimated counts for genes meeting the declared expression filter (at least 10 estimated counts in at least two models). `rna_vst.rds` is the corresponding variance-stabilised matrix in R binary format, readable with `readRDS`.

The matched transcript map links each quantified transcript to its gene and reference sequence region. Version suffixes are retained or removed only as declared by the reference script. Primary-assembly and alternative-locus transcripts are annotated; the gene matrix preserves the quantification target set. The reconciliation table documents target coverage and abundance conservation. Gene IDs and symbols are different identifiers: use the supplied map and the recorded symbol-aggregation rule for cross-layer comparisons.

`rna_qc_metrics.csv` reports one row per RNA model: `pseudoalign` is the percentage of kallisto-processed fragments pseudoaligned, `n_processed_fragments` is the number of fragments processed, `assigned_gene_counts` sums gene-level estimated counts after the global expression filter, and `detected` counts globally retained genes with estimated counts of at least 1. `subtype` and `site` are the original analysis labels, mapped to the cohort metadata above. A paired read pair counts as one fragment, not two reads.

## Protein records

`prot_abundance_matrix.csv` contains supplied log2 protein abundance normalised using the pooled internal standard. Rows use a representative gene symbol or `SYMBOL|UNIPROT` for secondary rows sharing a symbol. The annotation table's `row` column matches the matrix's `protein` column exactly. No missing-value imputation is applied to the released matrix.

In `prot_qc.csv`, `symbol` and `uniprot` retain source identities; `symbol_representative` identifies the row selected for symbol-based analyses using completeness, peptide support, and search q-value. `n_peptides`, `n_unique_peptides`, and `npeptides_quant` are identified, unique, and quantified peptide counts. `qvalue`, `cv_replicates`, and `isoDoping` retain source-search/precision/isotope-doping annotations; their detailed upstream definitions must be confirmed with the proteomics laboratory. The reported CV is expressed as a percentage.

`present_n_lines` and `present_n_plex` count measured resource models and plexes. `pass_presence50` requires at least 16 of 31 protein models. `complete_case` requires all 31. `zero_plex` flags identified features without any quantified value. `prot_block_missingness.csv` adds lists of present/absent plexes, `absent_from_any_plex`, and `pct_lines_missing`. `prot_sample_qc.csv` reports measured feature count and percentage missing for each model. `prot_feature_accounting.csv` defines each feature set (`level`), its size (`n`), and the precise `definition`.

## Variant records

`wes_mutations_filtered.csv` contains retained tumour-only variant candidates, not validated somatic variants. `Hugo_Symbol` is the annotated gene; `Chromosome`, `Start_Position`, and `End_Position` are GRCh38 MAF coordinates (one-based); `Reference_Allele` and `Tumor_Seq_Allele2` encode the reference and alternate alleles. Indel coordinates and dash alleles follow MAF convention as implemented by the archived conversion and the documented recovery script.

`Variant_Classification` describes coding consequence; `Variant_Type` distinguishes SNP/insertions/deletions/multinucleotide changes. `Consequence` and `IMPACT` preserve VEP annotation. `HGVSp_Short` is the annotated or reconstructed protein-change label; the `hgvsp_*` flags distinguish canonical notation from a reconstruction. Such a label is not an independently validated protein consequence.

`dbSNP_RS`, `Existing_variation`, and `CLIN_SIG` retain source database annotations. Clinical-significance annotations can be absent or outdated and are not clinical classifications made by this project. `pop_af_max` is the maximum of available source population-frequency fields; a missing source frequency must not be interpreted as proven absence from the population. `vaf` is the tumour alternate-allele fraction. Any VAF-based flag is a heuristic annotation, not proof of germline origin or loss of heterozygosity. `is_driver` denotes membership of the selected gene panel; `FILTER` records upstream filtering.

`wes_driver_tiers.csv` reports the prioritised subset using `gene`, `protein_change`, `variant_classification`, `vaf`, source clinical annotations, `tier`, `rationale`, and `context`. Patient and model identifiers specify the observation unit. Tier 1 represents the strongest gene/consequence prior used here, Tier 2 a weaker candidate prior, and Tier 3 unresolved candidates retained for inspection; consult the explicit per-row rationale and analysis script. These are heuristic tiers, not calibrated probabilities. Any summary that combines models from one patient must declare whether it uses a single selected model or the union across related models.

## Copy-number records

`wes_cnv_segments.csv` uses GRCh38 zero-based, half-open intervals (`chromosome`, `start`, `end`). `probes` and `weight` describe supporting CNVkit bins. `log2_raw` is the archived segment ratio; `log2c` is the legacy all-chromosome-centred value; `log2c_auto` uses the autosomal probe-weighted median and is preferred. These ratios do not supply absolute or allele-specific copy number.

`wes_cnv_fga.csv` reports segment count and assessed length in megabases. Fields beginning `fga_auto_` use autosomal segment length and the threshold in the field suffix; `frac_gain_auto` and `frac_loss_auto` use the primary absolute log2 threshold of 0.20. The corresponding fields without `_auto` retain legacy inclusion of chromosome X. `chrX_median_log2c` and `chrX_frac_altered` describe chromosome X separately. The operational `high_fga_flag` is a project threshold, not a calibrated genomic-instability classification.

`wes_cnv_arm_freq_patient.csv` reports HGSC gain/loss counts and percentages for each `arm`. `log2_threshold`, `arm_majority_frac`, `centring`, `arm_reference`, `fraction_denominator`, and `patient_rule` define the calculation. A segment contributes only its overlap with a reference arm; centromeric sequence is excluded. The arm-call denominator is the full annotated non-centromeric arm length. This differs from the assessed-segment-length denominator used for FGA. `n_lines` and `n_patients` specify denominator sizes; associated gain/loss fields count or give percentages of those units. Acrocentric p arms are omitted from the manuscript display. The threshold-sensitivity outputs remain available in the analysis project.

## Provenance and validation

`file_catalog.csv` records each packaged file's source, unit, and purpose. `field_inventory.csv` inventories literal columns, row counts, and missing values. Reference files retain target mapping, copy-number boundaries, and variant input/filter accounting when available. `validation.json` and `SHA256SUMS` document automated cross-file checks and exact packaged content. These checks establish internal consistency; they do not substitute for experimental assay validation or stock authentication.
