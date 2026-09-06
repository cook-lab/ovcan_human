# Decisions needed before the cluster pilots run (2026-09-06)

Answer inline or in a reply commit; the cluster agent will act on this file. Everything below is
prepared and staged on Nibi (`/project/6090753/active/ovcan_human/wes/molecular_extension_2026-09-06/`);
nothing has been submitted. See `INITIAL_RETURN.md`, `task_status.tsv` and `proposed_commands/`.

## Q1. Alignment generation for each pilot

Three GRCh38 alignments exist per model (details in `input_inventory.tsv`):

| Generation | Reads present | BQSR | Used by |
|---|---|---|---|
| A. `md.cram` (duplicate-marked) | all reads, genome-wide incl. off-target and flanks | no | MultiQC/mosdepth QC |
| B. published `recal.cram` | reads overlapping the 405 calling shards only | yes | named in retained VCF headers; 10 Oct-24 + 13 Oct-31 re-executions |
| C. Oct-24 `recal.cram` in `work/` | as B | yes | exact input converted to the manual-CNVkit BAMs |

Proposal: **A for MSI (MEX04) and read/locus review (MEX02/03/08)**; **C for MEX05/MEX06** so
SNP allele counts and the CNVkit log-ratios derive from identical reads.
Alternative: A for MEX05 as well (more flanking SNPs beyond ~100 bp of targets; un-recalibrated
base qualities; different read set from the CNR). **Choose: proposal / alternative.**

## Q2. Which jobs to authorise

| Script | Task | Models | Resources | Authorise? |
|---|---|---|---|---|
| `mex02_03_locus_review.sbatch` | MEX02/03/08 | 25 regions, 25 loci; 5 normals for bedcov | 1 job, 2 CPU, 8 GB, <2 h | yes / no |
| `mex04_msisensor2_pilot.sbatch` | MEX04 | TOV21G, TOV2414, TOV3392D, OV3331 | 4 × (4 CPU, 16 GB, ≤2.75 h) | yes / no |
| `mex05_mutect2_germline_sites.sbatch` | MEX05 | TOV81D, OV2085, TOV2929D, OV3331, OV1369-R2 | 5 × (4 CPU, 16 GB, ≤7 h) | yes / no |
| `mex05b_normals_mutect2_for_mappingbias.sbatch` | optional mapping bias | 5 public normals | 5 × (4 CPU, 16 GB, ≤7 h) | yes / no |
| `mex06_purecn_pilot.sbatch` | MEX06 (after MEX05) | same 5 models, 2 fits each | 5 × (4 CPU, 16 GB, ≤3 h) | yes / no |

If MEX04 site yield is interpretable, extend MSIsensor2 to all 23 (19 more tasks of the same
size)? **yes / no / decide after pilot.**

## Q3. PureCN fit settings (MEX06)

Drafted: `--fun-segmentation none` (import target-only CNR + corrected CNS as SEG), `--model
betabin --post-optimize --sex F --seed 20260906`; fit A purity 0.30–0.99, ploidy 1.5–6, max CN 12;
fit B purity 0.80–0.99, ploidy 1.5–8, max CN 16. Mapping bias only if Q2's optional normals run
is approved (labelled as process-matched public normals, not donor normals).
**Accept, or specify other grids/flags.**

## Q4. scarHRD (MEX07)

Not available on Nibi. Proposal: install from a pinned GitHub commit into `~/R` only after MEX06
fits pass QC; report LOH/TAI/LST components without a clinical threshold.
**Approve the install when the time comes? yes / no / name a preferred implementation.**

## Q5. MEX10 and MEX11 locations

RNA alignments/junctions and proteomics PSM/reporter records are not in the WES directories and
were not searched. **Provide paths (cluster or elsewhere) or mark these as out of scope for the
cluster.**

## Q6. Repository hygiene

`docs/cluster/recovery/2026-09-05/` (first recovery report, model status, QC metrics, inventory)
and `docs/cluster/recovery/2026-09-06/FOLLOWUP_INVENTORY.md` (+ TSVs: workflow-revision
verification, reference hashes, antitarget cause, CNVkit input chain, missing job logs) are pushed
with this directory. Several items still listed as outstanding in `2026-09-06/FOLLOWUP.md` are
answered there. **Confirm the local analysis should consume them, and whether any should be
trimmed before the repository stays public** (they contain cluster paths and lab usernames, no
credentials or clinical identifiers).

## Already settled, for reference

- TOV81D: WES sample id `TOV81D_P23` used as-is; the `TPV81D_23-pool` alias is a provider question.
- The 22 legacy exomes remain excluded (author decision 2026-09-05).
- No MSI controls and no matched donor normals exist in the searched roots; pilots return
  continuous scores and fit alternatives with uncertainty, not classifications.
