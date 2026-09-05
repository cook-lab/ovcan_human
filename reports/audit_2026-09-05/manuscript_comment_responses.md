# Responses to the version 5 Word comments

All 34 comments were extracted from the Word document. Comment IDs 2 and 3 belong to tracked edits rather than comments; no comment was omitted. Exact anchors and original text are preserved in v5_comments_anchored.json. The manuscript was revised from the current Word text, and v5 was preserved.

## Comment 0

**Status:** Resolved

**Comment:** I'd like this to be reconciled. No reason to only have one analysis type for one of the models

**Response:** Recovered the local TOV3121D annotated VCF and added its candidate-variant profile. RNA/proteomics overlap remains explicitly described. The manuscript, Table 1, Table 2, and Figure 4 now use 23 models for both exome analyses.

## Comment 1

**Status:** Resolved

**Comment:** It's not just re-processing. We *generated* the RNA-seq and proteomics fresh. Reprocessing was only the WES

**Response:** Abstract and Background distinguish newly generated RNA/protein profiles from reprocessed existing exome data. The Word edit from reprocessed to processed was accepted as the baseline.

## Comment 4

**Status:** Resolved

**Comment:** This term isn't intuitive to me. Patient? Donor sample?

**Response:** Replaced donor-family with patient of origin and models derived from the same patient.

## Comment 5

**Status:** Resolved

**Comment:** Maybe vague, maybe fine. I'm luke warm about this statement

**Response:** Named abundance matrices, variant annotations, copy-number segments, and model metadata instead of generic readiness claims.

## Comment 6

**Status:** Resolved

**Comment:** Tense is weird. Should we describe it as past tense?

**Response:** Past tense describes data generation and validation; present tense describes the resource contents and reuse.

## Comment 7

**Status:** Resolved

**Comment:** I prefer oxford comma throughout

**Response:** Applied serial/Oxford commas throughout the revised scientific prose.

## Comment 8

**Status:** Resolved

**Comment:** Ensure all references are accurate. I will re-do them with a proper reference manager once drafting is complete

**Response:** Verified all 17 original DOI records through Crossref and Europe PMC. Restored omitted punctuation in four titles. Removed two unused references after deleting the unsupported external correlation benchmark comparison. Verification records are in reference_verification.json.

## Comment 9

**Status:** Resolved

**Comment:** What "measurements" is this referring to? Maybe we can generalize "use" or "data"?

**Response:** Replaced nonspecific measurements with RNA, protein, and genomic data.

## Comment 10

**Status:** Resolved

**Comment:** I see the point, but the use of "measurement" is kind of awkward to me. It could just be a me thing, but I feel like the more general "finding" "data" or something is clearer

**Response:** Rewrote the cross-study comparison in terms of datasets and experimental conditions.

## Comment 11

**Status:** Preserved

**Comment:** This is described clearly--great

**Response:** Retained the clear statement that 13 models are sublines from five patients.

## Comment 12

**Status:** Resolved

**Comment:** Again, I don't really know what this means when I read it. Can it be more clear/explicit?

**Response:** Explained common patient identifiers and selection of one model per patient in ordinary language.

## Comment 13

**Status:** Resolved

**Comment:** CRITICAL - this must be fixed. I would rather just change the metadata table so that there's not an awkward reference to a file with a required filter to understand the cohort.

**Response:** Cohort is described directly as 42 models. The release includes a separate 42-row model table without requiring an inclusion filter. Original source metadata remains an internal provenance input.

## Comment 14

**Status:** Resolved

**Comment:** Again, if they weren't included, let's remove them from the metadata and avoid the clunkiness of describing this altogether (eg. no justification for why it's excluded)

**Response:** Removed the unrelated low-grade serous cohort explanation from the manuscript and biological release metadata. The technical bridge VOA3993 is documented in proteomics Methods because it actually contributes to repeatability validation.

## Comment 15

**Status:** Resolved

**Comment:** See previous flags--I dont' know what this means. Will stop flagging this phrase

**Response:** Patient grouping and the representative selection rule are explained explicitly. Mutation/CNV arm frequencies use patient unions; FGA uses within-patient means, avoiding a blanket claim that all summaries use one representative.

## Comment 16

**Status:** Draft supplied; confirmation needed

**Comment:** TODO - you can identify the publications where these lines were generated and get the methods details to at least include a draft here for collaborators to confirm

**Response:** Located CRCHUM derivation papers from 2000, 2008, 2012, 2015, and 2020, plus BIN67/COV434 sources. author_confirmation.md contains model mappings and draft culture text. Reported 5% vs 7% oxygen and alternative BIN67 media require actual-harvest confirmation. VOA-specific derivation/culture records and TOV3121D derivation remain M01.

## Comment 17

**Status:** Incorporated

**Comment:** RNA-Seq libraries were generated from 250 ng of total RNA as follows: mRNA enrichment was performed using the NEBNext Poly(A) Magnetic Isolation Module (New England BioLabs). cDNA synthesis was achieved with the NEBNext RNA First Strand Synthesis and NEBNext Ultra Directional RNA Second Strand Synthesis Modules (New England BioLabs). The remaining steps of library preparation were done using and the NEBNext Ultra II DNA Library Prep Kit for Illumina (New England BioLabs). Adaptors and PCR primers were purchased from New England BioLabs. The libraries were normalised and pooled and then denatured in 0.02 N NaOH and neutralised using HT1 buffer. The pool was loaded at 200pM on an Illumina NovaSeq S4 lane using Xp protocol as per the manufacturer’s recommendations. The run was performed for 2 × 100 cycles (paired-end mode). A phiX library was used as a control and mixed with libraries at a 1% level. Base calling was performed with RTA v3. Program bcl2fastq2 v2.20 was then used to demultiplex samples and generate fastq reads.

**Response:** Inserted the supplied 250 ng NEBNext library preparation, NovaSeq S4 Xp 2x100 sequencing, PhiX, RTA3, and bcl2fastq2 details, with grammatical and unit formatting edits. RNA extraction/integrity remains M02.

## Comment 18

**Status:** Resolved

**Comment:** Is deposited where?? As in ensembl made this available? Remove if not relevant

**Response:** Clarified that the exact-reference transcript map and checksums accompany this resource. Removed the ambiguous claim of an unspecified deposit by Ensembl.

## Comment 19

**Status:** Resolved analytically

**Comment:** Is this saying that our setup wasn't ideal, missed transcripts, and that it would be appropriate to re-analyze? If so, remove this. I don't want to raise an obvious (fixable) flag

**Response:** Recovered the exact Ensembl93 cDNA reference, matching all 185299 versioned targets and lengths in every library. Replaced the mismatched release105 map and reran affected analyses. All transcript estimates are represented; the avoidable mismatch paragraph is removed.

## Comment 20

**Status:** Draft supplied; confirmation needed

**Comment:** See if you can pull appropriate methods from a publication by Gregg Morin and Gian Negri (UBC)

**Response:** Located Orlando et al. 2020, DOI10.7554/eLife.59073, with Negri and Morin and an 11-plex cell-pellet protocol. Provided a concrete, explicitly provisional methods draft for laboratory confirmation (M03), rather than declaring an unrelated protocol to be the one used here.

## Comment 21

**Status:** Author TODO M04

**Comment:** This info is on the cluster where the pipeline was run. Add it to a running to-do list and I can grab it

**Response:** Consolidated cluster provenance requests: run versions/configuration, WES experimental methods, capture design, references, annotation, parameters, and QC.

## Comment 22

**Status:** Resolved

**Comment:** Make sure this makes sense. I'm not very familiar with it and I'm not sure what it means from the text alone

**Response:** Explained reconstructed amino-acid consequence labels with H214fs as an example, retained non-canonical flags for incomplete strings, and named genomic coordinates and transcript identifiers as the unambiguous reference.

## Comment 23

**Status:** Resolved

**Comment:** Sentence structure feels clunky?

**Response:** Rewrote tumour-only candidate language in direct prose while retaining the necessary rare-germline limitation.

## Comment 24

**Status:** Resolved

**Comment:** This paragraph hedges a lot. It says more what it's not than what it is.

**Response:** Explained what each evidence tier prioritises. Corrected 14 genes to 19 selected genes. Removed repeated defensive descriptions; retained the specific distinction between functional evidence and somatic origin.

## Comment 25

**Status:** Resolved; author TODO M05

**Comment:** This should be flagged as a TODO, not as manuscript text. The manuscript should have no "work in progress" prose. If there are things we should do to finalize the analysis and draft, flag it explicitly as "TODO"

**Response:** Moved missing original capture-design documentation into an explicit TODO. Recovered normal accessions and 290475-bin CNVkit target footprint are reported as available records.

## Comment 26

**Status:** Resolved

**Comment:** This reads weird to me. I think I like using language more about the biology/interpretation. Eg. here I would refer to transcriptional patterns across models

**Response:** Described principal components as transcriptional patterns across models. Explained raw component-score variance partitioning in Methods.

## Comment 27

**Status:** Resolved

**Comment:** Perhaps unnecessary unless this comes up later. The previous sentence delivers the message

**Response:** Removed the repetitive histotype-associated versus histotype-specific sentence. Retained the concrete centre/histotype design limitation once in the relevant validation paragraph.

## Comment 28

**Status:** Resolved

**Comment:** Could the meaning of this be clearer in the prose? ie. what are we doing filtering down? The next sentence reveals what filtering corrected, but the "issue" is not really clear

**Response:** Explained caller-quality failure, population-frequency filtering, and consequence filtering before reporting the 582474 to16081 to6194 cascade.

## Comment 29

**Status:** Resolved

**Comment:** All references could just be referred to with in text author citation + DOI for me to manually enter

**Response:** Replaced superscript citation tokens with in-text author/year and DOI references for later reference-manager entry.

## Comment 30

**Status:** Assessed; proposed generalisation unsupported

**Comment:** Could a more general claim be made comparing HGSC to the other histotypes more broadly? They should all have relatively less aberrations

**Response:** The corrected copy-number profiles include TOV3392D clear cell with FGA0.671, above the HGSC model median0.635. Replaced the single low-grade comparison with a concise description of heterogeneous other models and limited non-serous sample sizes.

## Comment 31

**Status:** Resolved; documentation TODO

**Comment:** Is this true? The external normals were generated with the same kit

**Response:** Accepted the author statement that the same capture kit was used. Requested its name/version and original design records without asserting incompatible designs.

## Comment 32

**Status:** Resolved

**Comment:** This is a little cliche. I think it's stronger to directly state the interpretation. ie. why might there be discordance between the two layers?

**Response:** Removed the rhetorical BIN67 interchangeability example and mechanism implication. Cross-layer differences are described as compatible with several biological and technical processes without assigning a causal explanation.

## Comment 33

**Status:** Resolved; author TODO M09

**Comment:** The contributing labs may very well have STR and mycoplasma testing. Do not anywhere state that it wasn't performed the study. Add a flag for me to confirm

**Response:** Removed the unsupported assertion about absent testing. Requested STR/mycoplasma records from contributing laboratories and distinguished them from external reference availability.

## Comment 34

**Status:** Resolved

**Comment:** I don't like the idea of reference to specific data files. That doesn't feel normal. Can't assume file naming is conserved/static.I think I would state this more along the lines of "Metadata is available with the resource to provide context about each model's..."

**Response:** Usage Notes now describe metadata contents without tying the prose to a working filename. Release records supply the file manifest and dictionary.

## Comment 35

**Status:** Resolved; deployment verification TODO D02

**Comment:** Don't refer to the file. This is hosted at https://www.cooklab.ca/ovcan_viewer

**Response:** Used the supplied cooklab.ca/ovcan_viewer URL. Requested that the hosted browser be checked against the corrected release; no public deployment is implied.
