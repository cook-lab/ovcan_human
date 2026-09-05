# Author confirmations for manuscript version 6

This document collects information that needs laboratory or author confirmation before submission. The manuscript uses matching TODO identifiers at the relevant locations. Published protocols below are draft source material, not evidence that the same conditions were used for the samples in this study.

## Culture methods and model derivation

**M01 — Confirm model-specific culture and harvest conditions.** The relevant publications were located and their methods inspected. The main uncertainty is which published protocol was followed at each assay harvest. Oxygen and serum conditions differ across papers, including papers from the same group.

Draft for CRCHUM models, pending confirmation:

> Cells were maintained in OSE medium supplemented with 10% fetal bovine serum, 0.5 µg/mL amphotericin B, and 50 µg/mL gentamicin at 37 °C in 5% CO2 and [confirm oxygen concentration for each model]. Cells were passaged before confluence, and assay-specific passage numbers were recorded in the model metadata.

The nine models from patients 1369, 2295, and 3133 are described by [Létourneau et al. 2012](https://doi.org/10.1186/1471-2407-12-379): OV1369-R2, TOV1369, OV2295, OV2295-R2, TOV2295-R, OV3133-R, OV3133-R2, TOV3133D, and TOV3133G. The paper reports 5% O2 and 5% CO2 with the supplemented OSE medium above. Solid-tumour cultures were established by scraping, and ascites cultures from the centrifuged cellular fraction.

[Fleury et al. 2015](https://doi.org/10.18632/genesandcancer.76) describes OV866-2, OV4453, OV4485, TOV3041G, and TOV3291G among its new models. It reports 7% O2 and 5% CO2 with supplemented OSE medium. [Sauriol et al. 2020](https://doi.org/10.3390/cancers12082222) describes OV2085, TOV2414, TOV2835EP, TOV2881EP, TOV2929D, TOV3121EP, OV3291, OV3331, and TOV3392D among the models relevant here, with 37 °C, 7% O2, and 5% CO2. The latter paper describes OV3291 as not immortalised at the reported passages; confirm its status for the present samples. TOV3121D shares the patient number with TOV3121EP, but its specific derivation should be confirmed rather than inferred from the name.

[Ouellet et al. 2008](https://doi.org/10.1186/1471-2407-8-152) is the derivation reference for OV1946. [Provencher et al. 2000](https://doi.org/10.1290/1071-2690%282000%29036%3C0357%3ACOFNEO%3E2.0.CO%3B2) describes OV90, TOV81D, TOV21G, and TOV112D. The later molecular reassignment of TOV112D and COV434 is documented by [Karnezis et al. 2021](https://doi.org/10.1016/j.ygyno.2020.12.004). Historical histotype labels should not be presented as current pathological consensus.

Draft for BIN67, pending confirmation:

> BIN67 cells were cultured in [confirm medium formulation and serum concentration] at [confirm temperature, oxygen, and CO2 concentrations].

The Ottawa study by [Gamwell et al. 2013](https://doi.org/10.1186/1750-1172-8-33) reports DMEM supplemented with 20% fetal calf serum and enriched with 20% Ham's F12 medium. A later study by [Orlando et al. 2020](https://doi.org/10.7554/eLife.59073) instead used RPMI-1640 with 10% fetal bovine serum for BIN67 and COV434 at 37 °C and 5% CO2. These are documented alternatives, not interchangeable assumptions about the current dataset. [Zhang et al. 2000](https://doi.org/10.1093/molehr/6.2.146) is an earlier COV434 characterisation paper under its historical granulosa-cell annotation.

The 11 VOA models (VOA10816, VOA12539, VOA14993, VOA295, VOA4841, VOA6861, VOA4395, VOA8762, VOA8771, VOA5217, and VOA5436) require the contributing laboratory's derivation references or unpublished derivation descriptions and culture records. Exact-name literature searches did not establish a protocol for every model. Do not assign culture conditions from another histotype panel by analogy.

For all models, confirm medium supplier and formulation, supplements, temperature, oxygen, CO2, passage at each harvest, identity of the harvested stock, and whether the RNA and protein aliquots came from the same culture. Existing passage metadata is useful but does not answer all of these questions.

## RNA library preparation

**M02 — RNA extraction and integrity.** The detailed library-preparation paragraph supplied in Word comment 17 is incorporated into v6. It includes 250 ng total RNA, the NEBNext modules and kit, NovaSeq S4 Xp sequencing at 2 × 100 cycles, 1% PhiX, RTA v3, and bcl2fastq2 v2.20. Confirm the RNA extraction protocol, integrity metric and acceptance criterion, and any pre-library quality checks. Do not infer these from the library kit.

## Proteomics protocol draft

**M03 — Confirm the actual acquisition and processing protocol.** A directly relevant primary study involving Gian Luca Negri and Gregg Morin is [Orlando et al. 2020](https://doi.org/10.7554/eLife.59073), which used 11-plex TMT on cell pellets. A short draft for laboratory review follows:

> Proteins from frozen cell pellets were digested with trypsin. Peptides were labelled with 11-plex tandem mass tags, pooled, separated by high-pH reversed-phase chromatography into 48 fractions, and concatenated into 12 fractions. Desalted peptides were analysed with an Easy-nLC 1000 coupled to an Orbitrap Fusion operating in MS3 mode. Spectra were searched with Sequest HT in Proteome Discoverer 2.1.1.21 against the UniProt reference proteome dated 3 August 2018. Precursor and fragment tolerances were 20 ppm and 0.8 Da. Oxidation and N-terminal acetylation were variable modifications; carbamidomethylation and TMT labelling were fixed modifications. Percolator filtering retained peptide-spectrum matches at a false-discovery rate below 1%.

Every parameter in this draft must be checked against the present acquisition records. Its cited study retained proteins supported by at least two unique peptides, whereas the current source table contains lower-support features. The data must not be described as meeting that criterion without verification. Confirm lysis, reduction/alkylation, digestion amounts and timing, LC gradient, acquisition settings, protein-level false-discovery control, reporter interference thresholds, impurity correction, shared-peptide handling, protein summarisation, channel mapping, internal-standard composition, and all normalisation steps. Define exactly what the supplied CV column measures; bridge-derived repeatability is separately calculated.

Another relevant Morin/Negri paper is [Ji et al. 2022](https://doi.org/10.1002/path.6006), which profiles archival ovarian tumour tissue with SP3-CTP. It uses tissue-specific processing and is less directly applicable to cell pellets. Its accession PXD032355 belongs to that published study and must not be used as an accession for the present data.

## Exome run records

The original WES analysis is accessible on the cluster (author update, 5 September 2026). Use the [prioritized cluster retrieval checklist](wes_cluster_retrieval.md) and its recorded path/model tables for M04/M05. These are gaps in the local archive to retrieve or assess, not presumed irrecoverable data.

**M04 — Cluster run provenance.** Retrieve the nf-core/Sarek version and complete run configuration; sequencing instrument, read length, library kit, original capture kit and target BED; alignment-reference checksum; and existing sequencing/coverage/contamination QC. All 23 archived VCF headers establish GATK Mutect2 and FilterMutectCalls 4.5.0.0, VEP 113.0 on GRCh38.p14, 1000 Genomes phase 3, gnomAD exome/genome 4.1, ClinVar 202404, and dbSNP 156. These recovered versions are incorporated in the manuscript and no longer remain open requests. The recovered TOV3121D VCF resolves the earlier 22-versus-23 profile discrepancy locally; it does not remove the need to archive upstream pipeline provenance.

**M05 — Reference exomes and capture design.** The user states that cell-model and reference exomes used the same kit. Record the kit name and revision, the original vendor target files, and reference processing versions. The five accessions are recovered and included in v6: SRR4039087, SRR4039088, SRR4039089, SRR4039096, and SRR4039097 (PRJNA339046). The archive also contains identical CNVkit-derived target definitions with 290,475 bins; these do not establish the original vendor capture design. This is a documentation task, not a claim that incompatible capture designs were used. Retrieve per-exome target coverage and depth summaries to support genomic technical validation; positive-control alterations alone do not validate sensitivity throughout the exome.

## External references and stock records

**M06 — Cellosaurus release and retrieval documentation only.** DepMap provenance is resolved: [DepMap 24Q4 Public, version 1](https://doi.org/10.25452/figshare.plus.27993248.v1), published 10 December 2024. All four primary local inputs (`Model.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsSomaticMutationsMatrixHotspot.csv`, and `OmicsSomaticMutationsMatrixDamaging.csv`) matched the official version-1 metadata by MD5 and byte count. The exact filenames, MD5/SHA-256 values, and cached/derived reference inputs are recorded in `output/external_reference_provenance.csv` and its JSON summary; the official metadata response is archived in `data/reference/depmap_24Q4_figshare_v1_metadata.json`. The manuscript now cites the persistent dataset record.

The documented historical Cellosaurus snapshot is 23 July 2026, supported by `reports/review_code.md` and the dated external-data workflow. This is documentary evidence of the snapshot date, not an independently verified HTTP retrieval timestamp. Cached responses retain individual entry-version and last-updated fields, but neither the database release identifier nor original HTTP retrieval headers were saved. Confirm the database release and original retrieval documentation only if those records remain available. Do not infer a database release from entry versions or filesystem modification times.

**M09 — Stock authentication and contamination tests.** Obtain STR and mycoplasma results from each contributing laboratory. Record test method, date, result, and relationship to the aliquots profiled here. The manuscript no longer claims testing was not performed. The presence of a Cellosaurus reference STR profile is distinct from a comparison performed on a study stock.

## Ethics and author statements

**M07 — Ethics.** Confirm the source-tumour consent framework and approvals, the current institutional determination, and the permitted sharing of genomic data. Published derivation approvals are useful leads but should not be copied as approvals for this study without verification. Complete the journal's Human Data Checklist.

**M08 — AI assistance.** Confirm tool names, versions or access dates where appropriate, the scope of assistance, and the authors' final review and responsibility. Do not assert that the authors have performed checks that are still pending.

**A01–A06 — Authors and declarations.** Confirm author order and affiliations; previous publications using these exact datasets; individual contributions; competing interests; acknowledgements; and funding with grant identifiers. A02 is important for judging the new data contribution of a Data Descriptor.

## Public records and submission

**D01 — Data deposits.** Assign accessions and persistent identifiers for reads, spectra/search output, processed matrices, and metadata. Add licences and any legitimate access conditions. Check that deposited content exactly matches the versioned release manifest and checksums. Supply anonymous reviewer access appropriate to the repository.

**D02 — Browser.** The user-supplied address is [cooklab.ca/ovcan_viewer](https://www.cooklab.ca/ovcan_viewer). Confirm the hosted site's dataset matches the corrected final release and archive that version. The manuscript revision does not imply a deployment has occurred.

**D03 — Code.** Publish the intended repository and immutable revision, then archive the release with a DOI. Run the documented workflow from the archived inputs and check the generated manifest before submission.

The current [Scientific Data submission guidelines](https://www.nature.com/sdata/publish/submission-guidelines) request a descriptive title no longer than 110 characters, recommend an abstract within 170 words, and require the Data Descriptor section structure. Data and methods should be described for reuse, with validation tied to technical quality. v6 follows those conventions and retains four main figures. The [editorial policies](https://www.nature.com/sdata/policies/editorial-and-publishing-policies) also make previous use or disclosure of the same data relevant to suitability. Final repository and human-data arrangements must follow the [data policies](https://www.nature.com/sdata/policies/data-policies). Guidance checked 5 September 2026.
