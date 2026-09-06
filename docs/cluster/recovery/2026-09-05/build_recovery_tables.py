#!/usr/bin/env python3
"""Build the WES recovery deliverables from the evidence copied into
/project/6090753/active/ovcan_human and from read-only stat() of the large
originals that were deliberately not copied.  Standard library only."""
import csv, glob, json, os, re, sys, datetime, collections

NEW = "/project/6090753/active/ovcan_human"
W = NEW + "/wes"
SRC = "/project/6090753/active/ovcan_gq_analysis"
RES = SRC + "/wes_sarek_results_20251022"
R00 = "/project/6090753/active/ovcan_gq/ammmasson_OvCAN_WES_bioinformatics_analysis_R004741"
REPO = os.path.expanduser("~/ovcan_human")
OUT = REPO + "/docs/cluster/recovery/2026-09-05"
os.makedirs(OUT, exist_ok=True)
os.makedirs(W + "/inventory", exist_ok=True)
NOW = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
RUN_VCF = "nf-core/sarek 3.5.1 run angry_allen (session 26cf85bf-0ac3-4237-8f48-c36b2ab20bd2, Slurm 3273686, completed 2025-10-23 23:48:40 EDT)"
RUN_OCT24 = "nf-core/sarek 3.5.1 run angry_kilby (session 26cf85bf, Slurm 3377445, 2025-10-24, failed at GETPILEUPSUMMARIES TOV3133D_P66)"
RUN_OCT31 = "nf-core/sarek 3.5.1 run disturbed_kilby (session 1d80816f-4d0e-4a07-8fa8-818e7b7f9db6, Slurm 3731446, 2025-10-31, killed by 12h59m walltime)"
MANUAL_CNV = "manual CNVkit 0.9.10 batch (apptainer image cnvkit-0.9.10--pyhdfd78af_0.img), outputs dated 2025-11-14 to 2025-12-03"

def tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if v is None else v for v in r])
    print("wrote", path, len(rows), "rows")

def st(path):
    try:
        s = os.stat(path)
        return s.st_size, datetime.datetime.fromtimestamp(s.st_mtime, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"), "accessible"
    except PermissionError:
        return None, None, "permission denied"
    except FileNotFoundError:
        if os.path.islink(path):
            return None, None, "broken symlink"
        return None, None, "not found in searched roots"

# ---------------------------------------------------------------- models
models = list(csv.DictReader(open(REPO + "/reports/audit_2026-09-05/wes_cluster_models.csv")))
assert len(models) == 23
SID = {m["cell_line"]: m["cnv_sample_id"] for m in models}
SIDS = sorted(SID.values(), key=len, reverse=True)
patients = {m["cell_line"]: m["patient_id"] for m in models}
assert len(set(patients.values())) == 16, set(patients.values())
hgs_pat = {m["patient_id"] for m in models if m["histotype"] == "HGS"}
assert len(hgs_pat) == 11, hgs_pat

manifest = {}
for line in open(NEW + "/MANIFEST.sha256"):
    h, p = line.rstrip("\n").split("  ", 1)
    manifest[p] = h
vcfcheck = {r["sample"]: r for r in csv.DictReader(open(
    "/tmp/claude-3055330/-home-dcook/b0c11c8f-b271-4cf7-9676-37f8757b0bda/scratchpad/vcf_hash_check.tsv"), delimiter="\t")}
repo_manifest = {r["cell_line"]: r for r in csv.DictReader(open(REPO + "/output/wes_input_manifest.csv"))}

# lanes from fastp report names
lane = {}
for p in glob.glob(W + "/sarek_3.5.1_run/reports/fastp/*/*.fastp.json"):
    m = re.match(r"(.+)-(\d+)\.fastp\.json$", os.path.basename(p))
    lane[m.group(1)] = m.group(2)

# ---------------------------------------------------------------- QC parsing
qc = []
def add(model, run, metric, value, unit, denom, target, src, status="found", notes=""):
    qc.append([model, run, metric, value, unit, denom, target, src, status, notes])

TGT_BED = "intervals_sorted.bed (242,421 hg38 intervals lifted from SeqCap EZ Exome v3; 63,709,951 bp summed, 63,514,049 bp merged)"
TGT_CNV = "intervals_sorted.target.bed (290,475 CNVkit autobin-split target bins, sha256 2177970f...; 63,514,049 bp)"

def rel(p): return os.path.relpath(p, NEW)

for s in sorted(SID.values()):
    run_md = RUN_VCF
    # fastp
    fj = glob.glob(f"{W}/sarek_3.5.1_run/reports/fastp/{s}/{s}-*.fastp.json")
    if fj:
        j = json.load(open(fj[0])); sm = j["summary"]
        add(s, run_md, "fastp_total_reads_before_filtering", sm["before_filtering"]["total_reads"], "reads (R1+R2)", "all reads in provider BAM converted to FASTQ", "", rel(fj[0]))
        add(s, run_md, "fastp_total_reads_after_filtering", sm["after_filtering"]["total_reads"], "reads (R1+R2)", "reads passing fastp --length_required 15 (adapter trimming disabled)", "", rel(fj[0]))
        add(s, run_md, "fastp_q30_rate_after_filtering", round(sm["after_filtering"]["q30_rate"], 4), "fraction of bases", "bases after filtering", "", rel(fj[0]))
        add(s, run_md, "fastp_gc_content_after_filtering", round(sm["after_filtering"]["gc_content"], 4), "fraction", "bases after filtering", "", rel(fj[0]))
        add(s, run_md, "fastp_duplication_rate", round(j["duplication"]["rate"], 4), "fraction of reads", "fastp sequence-level estimate before alignment", "", rel(fj[0]))
        add(s, run_md, "fastp_insert_size_peak", j["insert_size"]["peak"], "bp", "fastp overlap-based estimate", "", rel(fj[0]))
        add(s, run_md, "fastp_read_length_mean_R1", sm["after_filtering"]["read1_mean_length"], "bp", "", "", rel(fj[0]))
    else:
        add(s, run_md, "fastp_total_reads_before_filtering", "", "", "", "", "", "not found in searched roots")
    # markduplicates
    mf = f"{W}/sarek_3.5.1_run/reports/markduplicates/{s}.md.cram.metrics"
    if os.path.exists(mf):
        lines = open(mf).read().splitlines()
        for i, l in enumerate(lines):
            if l.startswith("## METRICS CLASS"):
                hdr = lines[i + 1].split("\t"); val = lines[i + 2].split("\t"); d = dict(zip(hdr, val)); break
        for k, u in [("READ_PAIRS_EXAMINED", "read pairs"), ("UNPAIRED_READS_EXAMINED", "reads"), ("UNMAPPED_READS", "reads"),
                     ("READ_PAIR_DUPLICATES", "read pairs"), ("READ_PAIR_OPTICAL_DUPLICATES", "read pairs"),
                     ("PERCENT_DUPLICATION", "fraction"), ("ESTIMATED_LIBRARY_SIZE", "molecules")]:
            add(s, run_md, "markduplicates_" + k, d.get(k, ""), u, "all mapped reads (genome-wide, md.cram)", "", rel(mf), notes="GATK 4.5.0.0 MarkDuplicates, REMOVE_DUPLICATES=false")
    # samtools stats md and recal
    for stage, denom in [("md", "all reads, genome-wide (duplicate-marked CRAM)"),
                         ("recal", "reads retained after per-interval ApplyBQSR (interval-restricted CRAM used for calling)")]:
        sf = f"{W}/sarek_3.5.1_run/reports/samtools/{s}/{s}.{stage}.cram.stats"
        if os.path.exists(sf):
            sn = {}
            for l in open(sf):
                if l.startswith("SN\t"):
                    _, k, v = l.rstrip("\n").split("\t")[:3]; sn[k.rstrip(":")] = v
            for k in ["raw total sequences", "reads mapped", "reads mapped and paired", "reads properly paired", "reads duplicated",
                      "reads MQ0", "reads unmapped", "error rate", "average length", "insert size average", "insert size standard deviation",
                      "percentage of properly paired reads (%)", "bases mapped (cigar)", "average quality"]:
                if k in sn:
                    add(s, run_md, f"samtools_stats_{stage}_" + k.replace(" ", "_").replace("(%)", "pct").replace("(cigar)", "cigar"), sn[k],
                        "percent" if "%" in k else ("bp" if "insert" in k or "length" in k or "bases" in k else ("fraction" if "rate" in k else "reads")),
                        denom, "", rel(sf))
    # mosdepth summaries (md = genome-wide + region; recal = interval-restricted)
    for stage in ["md", "recal"]:
        ms = f"{W}/sarek_3.5.1_run/reports/mosdepth/{s}/{s}.{stage}.mosdepth.summary.txt"
        if os.path.exists(ms):
            for l in open(ms):
                f = l.rstrip("\n").split("\t")
                if f[0] == "total":
                    add(s, run_md, f"mosdepth_{stage}_mean_depth_genome", f[3], "x", f"whole reference ({f[1]} bp)", "", rel(ms))
                if f[0] == "total_region":
                    add(s, run_md, f"mosdepth_{stage}_mean_depth_target", f[3], "x", f"target bases ({f[1]} bp, --by intervals_sorted.bed)", TGT_BED, rel(ms),
                        notes="mean of per-base depth over the target footprint; mosdepth counts duplicates unless flagged (Sarek default)")
        rd = f"{W}/sarek_3.5.1_run/reports/mosdepth/{s}/{s}.{stage}.mosdepth.region.dist.txt"
        if os.path.exists(rd):
            want = {"1", "10", "20", "30", "50", "100"}
            for l in open(rd):
                f = l.rstrip("\n").split("\t")
                if f[0] == "total" and f[1] in want:
                    add(s, run_md, f"mosdepth_{stage}_fraction_target_bases_ge_{f[1]}x", f[2], "fraction", "target bases", TGT_BED, rel(rd))
    # contamination
    cf = f"{W}/sarek_3.5.1_run/variant_calling/mutect2/{s}/{s}.mutect2.contamination.table"
    if os.path.exists(cf):
        r = open(cf).read().splitlines()[1].split("\t")
        add(s, run_md, "gatk_calculatecontamination_estimate", r[1], "fraction", "GetPileupSummaries on af-only-gnomad.hg38 sites within calling intervals; tumour-only", "", rel(cf), notes="cross-sample contamination estimate, error=" + r[2] + "; not a stock-authentication result")
    # on-target (from Oct-24 Sarek CNVKIT_BATCH log; recal BAM converted from recal CRAM)
    ce = glob.glob(f"{W}/sarek_3.5.1_run/task_evidence/oct24_CNVKIT_BATCH_{s}.*/.command.err")
    if ce:
        txt = open(ce[0]).read()
        m = re.search(r"#bins=242421, #reads=(\d+)", txt); m2 = re.search(r"Percent reads in regions: ([\d.]+) \(of (\d+) mapped\)", txt)
        if m and m2:
            add(s, RUN_OCT24, "cnvkit_reads_in_target_bins", m.group(1), "reads", "mapped reads in recal BAM (interval-restricted)", "cnvkit.reference.target-tmp.bed (242,421 bins)", rel(ce[0]))
            add(s, RUN_OCT24, "cnvkit_percent_mapped_reads_in_target_bins", m2.group(1), "percent", f"{m2.group(2)} mapped reads in recal BAM (interval-restricted, so not a genome-wide on-target fraction)", "cnvkit.reference.target-tmp.bed (242,421 bins)", rel(ce[0]),
                notes="Sarek-internal CNVkit batch (flat reference) task log; recal CRAM excludes off-interval reads, so genome-wide on-target fraction must use md CRAM stats")
    # manual cnvkit coverage
    cn = f"{W}/cnvkit_0.9.10_manual/{s}_new/{s}.targetcoverage.cnn"
    if os.path.exists(cn):
        L = D = n = z = 0
        for i, l in enumerate(open(cn)):
            if i == 0: continue
            f = l.split("\t"); ln = int(f[2]) - int(f[1]); L += ln; D += float(f[4]) * ln; n += 1; z += (float(f[4]) == 0)
        add(s, MANUAL_CNV, "cnvkit_target_mean_depth_length_weighted", round(D / L, 2), "x", f"{L} bp in {n} target bins (CNVkit depth column)", TGT_CNV, rel(cn), notes="derived here from the retained coverage file; not a pipeline-reported number")
        add(s, MANUAL_CNV, "cnvkit_target_bins_zero_depth", z, "bins", f"{n} target bins", TGT_CNV, rel(cn), notes="derived here")
        cs = f"{W}/cnvkit_0.9.10_manual/{s}_new/{s}.call.cns"
        if os.path.exists(cs):
            add(s, MANUAL_CNV, "cnvkit_call_segments", sum(1 for _ in open(cs)) - 1, "segments", "", TGT_CNV, rel(cs))
    # variant counts
    bs = f"{W}/sarek_3.5.1_run/multiqc/multiqc_data/multiqc_bcftools_stats.txt"
    for r in csv.DictReader(open(bs), delimiter="\t"):
        if r["Sample"] == f"{s}.mutect2.filtered":
            for k in ["number_of_records", "number_of_SNPs", "number_of_indels", "number_of_MNPs", "tstv"]:
                add(s, run_md, "bcftools_stats_filtered_vcf_" + k, r[k], "ratio" if k == "tstv" else "records", "all records in mutect2.filtered.vcf.gz (PASS and non-PASS)", "", rel(bs))
    fs = glob.glob(f"{W}/sarek_3.5.1_run/reports/vcftools/**/{s}.mutect2.filtered.FILTER.summary", recursive=True)
    if fs:
        for l in open(fs[0]):
            f = l.split()
            if f and f[0] == "PASS":
                add(s, run_md, "vcftools_filtered_vcf_PASS_records", f[1], "records", "records in mutect2.filtered.vcf.gz", "", rel(fs[0]))

# normals
NORM = ["SRR4039087", "SRR4039088", "SRR4039089", "SRR4039096", "SRR4039097"]
for n_ in NORM:
    for kind, tg in [("target", TGT_CNV), ("antitarget", "intervals.antitarget.bed (42,709 bins)")]:
        cn = f"{W}/cnvkit_0.9.10_manual/OV1369_R2_P66_new/{n_}.sorted.{kind}coverage.cnn"
        L = D = n = z = 0
        for i, l in enumerate(open(cn)):
            if i == 0: continue
            f = l.split("\t"); ln = int(f[2]) - int(f[1]); L += ln; D += float(f[4]) * ln; n += 1; z += (float(f[4]) == 0)
        add(n_, MANUAL_CNV, f"cnvkit_{kind}_mean_depth_length_weighted", round(D / L, 2), "x", f"{L} bp in {n} {kind} bins", tg, rel(cn), notes="public CNV-reference exome (PRJNA339046); derived here from retained coverage file; identical copies exist in all 23 model directories")
        add(n_, MANUAL_CNV, f"cnvkit_{kind}_bins_zero_depth", z, "bins", f"{n} bins", tg, rel(cn), notes="derived here")
    for met in ["fastp_total_reads", "markduplicates_PERCENT_DUPLICATION", "samtools_stats_reads_mapped", "mosdepth_mean_depth_target", "gatk_calculatecontamination_estimate"]:
        add(n_, "bwa-mem2 alignment via cnvkit/normal_samples/run_bwa.sh (Slurm 4401780/4439216/4459055, 2025-11-13/14)", met, "", "", "", "", "",
            "not produced", "normal exomes were aligned outside Sarek; no FastQC/fastp/MarkDuplicates/mosdepth/flagstat output exists for them in the searched roots")

tsv(OUT + "/qc_metrics.tsv", ["model_or_reference", "run_id", "metric", "value", "unit", "denominator", "target_file_id", "source_path", "status", "notes"], qc)

# ---------------------------------------------------------------- model_status
def q(model, metric):
    for r in qc:
        if r[0] == model and r[2] == metric: return r[3]
    return ""
rows = []
# CRAM producing run by mtime bucket
for m in models:
    s = m["cnv_sample_id"]; cl = m["cell_line"]
    prov = f"{R00}/bam/{cl.replace('-R2','(R2)').replace('-R','(R)')}_{m['recorded_wes_passage']}.bam" if cl.startswith("OV") and "-R" in cl else None
    # provider BAM naming: OV1369(R2)_P66.bam, TOV2295(R)_P57.bam ...
    base = re.sub(r"-(R2?)$", r"(\1)", cl)
    prov = f"{R00}/bam/{base}_{m['recorded_wes_passage']}.bam"
    psz, pmt, pst = st(prov)
    cram = f"{RES}/preprocessing/recalibrated/{s}/{s}.recal.cram"
    csz, cmt, cst = st(cram); _, _, ist = st(cram + ".crai")
    crun = RUN_OCT31 if cmt and cmt >= "2025-10-31" else RUN_OCT24
    md = f"{RES}/preprocessing/markduplicates/{s}/{s}.md.cram"; msz, mmt, mst = st(md)
    wb = re.search(rf"(\.\./work/\S+/{s}\.recal\.bam)", open(REPO + "/docs/cluster/evidence/cnvkit_commands.txt").read()).group(1)
    wsz, wmt, wst = st(SRC + "/" + wb[3:])
    vc = vcfcheck[s]
    callcns = os.path.exists(f"{W}/cnvkit_0.9.10_manual/{s}_new/{s}.call.cns")
    unresolved = []
    if pst != "accessible": unresolved.append("provider BAM " + pst)
    if not callcns: unresolved.append("manual CNVkit call.cns missing")
    unresolved.append("retained recal CRAM is a re-execution (" + ("Oct-31 partial run" if crun == RUN_OCT31 else "Oct-24 run") + "), not the byte-identical Oct-23 CRAM that fed the retained VCF")
    unresolved.append("sequencing run/flowcell/date and library kit lot: not in cluster records; see WES 2014/2017/2018 spreadsheets and USB QC report (not parsed here)")
    rows.append([cl, m["histotype"], m["patient_id"], s, m["recorded_wes_passage"], lane.get(s, ""),
                 f"@RG ID:{s}_{lane.get(s,'')} PU:{lane.get(s,'')} SM:{s}_{s} LB:{s} PL:ILLUMINA (404 read groups after FASTQ splitting)",
                 prov, psz, pst,
                 RUN_VCF, "yes", "yes", vc["match"] == "YES", vc["cluster_uncompressed_sha256"],
                 cram, csz, cmt, cst, ist, crun, md, msz, mst,
                 SRC + "/" + wb[3:], wsz, wst,
                 "yes" if callcns else "no", q(s, "gatk_calculatecontamination_estimate"), q(s, "mosdepth_md_mean_depth_target"), q(s, "mosdepth_md_fraction_target_bases_ge_30x"),
                 q(s, "markduplicates_PERCENT_DUPLICATION"), "yes (MultiQC 2025-10-23; fastp, FastQC, MarkDuplicates, samtools stats, mosdepth, bcftools/vcftools stats)",
                 "; ".join(unresolved)])
tsv(OUT + "/model_status.tsv",
    ["cell_line", "histotype", "patient_id", "wes_sample_id", "wes_passage", "flowcell_lane_PU", "read_group_evidence",
     "provider_bam_path", "provider_bam_bytes", "provider_bam_status",
     "run_producing_retained_vcfs", "mutect2_filtered_vcf_present", "vep_annotated_vcf_present", "annotated_vcf_sha256_matches_repo_manifest", "annotated_vcf_sha256_uncompressed",
     "recal_cram_path", "recal_cram_bytes", "recal_cram_mtime_utc", "recal_cram_status", "recal_crai_status", "run_producing_current_recal_cram",
     "md_cram_path", "md_cram_bytes", "md_cram_status",
     "cnvkit_input_bam_workdir_path", "cnvkit_input_bam_bytes", "cnvkit_input_bam_status",
     "manual_cnvkit_call_cns_present", "contamination_estimate", "mosdepth_md_mean_target_depth_x", "fraction_target_ge_30x", "duplication_fraction",
     "qc_coverage", "unresolved_items"], rows)

# ---------------------------------------------------------------- additional models (separate table)
extra = []
usb = W + "/mcgill_r004741_provenance/USB_lists/"
completed = set(open(usb + "samples.completed").read().split()); missed = set(open(usb + "samples_missed").read().split())
for sname in open(usb + "samples_list").read().split():
    canon = re.sub(r"_P\d+$", "", sname); canon2 = re.sub(r"(OV\d+)(R2?)$", r"\1-\2", canon).replace("TOV2295R", "TOV2295-R")
    in23 = any(SID[c].replace("_R2_", "R2_").replace("_R_", "R_") == sname for c in SID)
    bai = f"{R00}/bam/{re.sub(r'(R2?)_P', r'(\1)_P', sname)}.bam.bai"; bam = bai[:-4]
    bsz, _, bst = st(bam); isz, _, ist_ = st(bai)
    extra.append([sname, canon2, "baseline 23" if in23 else "additional (not in manuscript WES cohort)", bam, bsz, bst, ist_,
                  "yes" if sname in completed else ("no (listed in samples_missed/missing_reports)" if sname in missed else ""),
                  "yes" if os.path.exists(f"{W}/mcgill_r004741_provenance/USB_cna_grch37/{sname}.cnvkit.cna.tsv") else "no",
                  "McGill GenPipes DNAseq (GRCh37; HaplotypeCaller+VQSR allSamples.hc.vqsr.vcf.gz split per sample; CNVkit germline; PCGR 0.9.2 by pmarquis, 2022-04)",
                  R00 + "/USB/{split,Reports,cna}"])
for sname, kit in [("BIN67", "SureSelectHumanAllExonV7.Target.bed (216,573 targets, GRCh37)"), ("CSCOE", "SureSelectHumanAllExonV7.Target.bed")]:
    extra.append([sname, sname, "additional (SCCOHT lines; different capture kit)", f"{R00}/cell_line/{sname}.vcf", *st(f"{R00}/cell_line/{sname}.vcf")[0:1], "VCF only, no BAM in searched roots", "", "", "yes (cnvkit.cna.tsv)", "McGill GenPipes DNAseq GRCh37 + PCGR; capture " + kit, R00 + "/cell_line"])
for f in sorted(glob.glob(W + "/mcgill_r004741_provenance/cna_grch37/*.cnvkit.cna.tsv")):
    sname = os.path.basename(f).split(".")[0]
    extra.append([sname, re.sub(r"p\d+(-\d+)?$", "", sname), "additional (VOA/OV866-2 lines; VCF+CNA+PCGR only)", f"{R00}/input/{sname}.vcf.gz", *st(f"{R00}/input/{sname}.vcf.gz")[0:1], "VCF only, no BAM in searched roots", "", "", "yes", "McGill GenPipes DNAseq GRCh37 + PCGR 0.9.2 (grch37)", R00 + "/{input,cna,PCGR}"])
tsv(OUT + "/additional_wes_models.tsv", ["sample_label", "canonical_guess_unverified", "cohort_status", "primary_file", "bytes", "bam_status", "bai_status", "pcgr_report_completed_2022", "cna_tsv_present", "pipeline_evidence", "location"], extra)

# ---------------------------------------------------------------- file inventory
inv = []
def model_of(name):
    for sid in SIDS:
        if sid in name.replace("(R2)", "_R2").replace("(R)", "_R").replace("R2_P", "_R2_P").replace("R_P", "_R_P"): return sid
    for n_ in NORM:
        if n_ in name: return n_
    return ""
def role_of(p):
    for k, v in [("/pipeline_info/", "workflow provenance"), ("/logs/", "launch/scheduler log"), ("/launch/", "launch configuration"),
                 ("/task_evidence/", "executed task script/log"), ("/multiqc/", "QC summary"), ("/reports/mosdepth", "depth QC"),
                 ("/reports/", "QC report"), ("/variant_calling/mutect2/", "variant calling output"), ("/annotation/", "annotated VCF"),
                 ("/preprocessing/recal_table", "BQSR table"), ("/csv/", "Sarek sample CSV"), ("cnvkit_0.9.10_manual/normal_samples", "CNV reference build script/log"),
                 ("cnvkit_0.9.10_manual/", "CNV output"), ("scratch_archive_extract/references", "reference identity"), ("scratch_archive_extract/wes_new", "superseded earlier Sarek run"),
                 ("scratch_archive_extract/wes_sarek_results_20251024", "superseded CNVkit-only test run"), ("scratch_archive_extract/", "original scratch root member"),
                 ("mcgill_r004741_provenance/", "original McGill/provider provenance"), ("/inventory/", "inventory")]:
        if k in p: return v
    return "other"
def orig_of(p):
    r = p.split("/", 1)[1] if p.startswith("wes/") else p
    if r.startswith("sarek_3.5.1_run/launch/nextflow_history.tsv"): return SRC + "/.nextflow/history"
    if r.startswith("sarek_3.5.1_run/launch/"): return SRC + "/" + r.split("/", 2)[2]
    if r.startswith("sarek_3.5.1_run/logs/"): return SRC + "/" + r.split("/", 2)[2]
    if r.startswith("sarek_3.5.1_run/task_evidence/"):
        d, f = r.split("/")[2], "/".join(r.split("/")[3:]); h = d.rsplit(".", 1)[1]; return SRC + "/work/" + h.replace("_", "/", 1) + "/" + f
    if r.startswith("sarek_3.5.1_run/"): return RES + "/" + r.split("/", 1)[1]
    if r.startswith("cnvkit_0.9.10_manual/"): return SRC + "/cnvkit/" + r.split("/", 1)[1]
    if r.startswith("scratch_archive_extract/"): return SRC + "/ovcan.tar.gz::ovcan/" + r.split("/", 1)[1]
    if r.startswith("mcgill_r004741_provenance/USB_lists/listgd"): return R00 + "/USB/Reports/listgd"
    if r.startswith("mcgill_r004741_provenance/USB_lists/"): return R00 + "/USB/" + r.rsplit("/", 1)[1]
    if r.startswith("mcgill_r004741_provenance/USB_cna_grch37/"): return R00 + "/USB/cna/" + r.rsplit("/", 1)[1]
    if r.startswith("mcgill_r004741_provenance/cna_grch37/"): return R00 + "/cna/" + r.rsplit("/", 1)[1]
    if r.endswith("SureSelectHumanAllExonV7.Target.bed"): return R00 + "/cell_line/SureSelectHumanAllExonV7.Target.bed"
    if r.startswith("mcgill_r004741_provenance/"): return R00 + "/" + r.rsplit("/", 1)[1]
    if r.startswith("inventory/"): return "generated 2026-09-05 from tar tvf of " + SRC + "/ovcan.tar.gz"
    return ""
rid = 0
for p in sorted(manifest):
    rid += 1
    sz = os.path.getsize(NEW + "/" + p)
    o = orig_of(p)
    inv.append([f"C{rid:05d}", model_of(os.path.basename(p)), "", "", "", role_of(p), o, NEW + "/" + p, "copied", sz, "sha256", manifest[p],
                "computed 2026-09-05 on the copy (cp -p/rsync -t preserve mtime); original not re-hashed", p, "MANIFEST.sha256", NOW, ""])
# large / not copied
def big(model, sample, passage, run, role, path, status_note="", checksum_ev="no checksum available in searched roots", notes=""):
    global rid; rid += 1
    sz, mt, stt = st(path)
    inv.append([f"L{rid:05d}", model, sample, passage, run, role, path, os.path.realpath(path) if stt == "accessible" else "", stt + ("; " + status_note if status_note else ""), sz, "", "", checksum_ev, "", "stat() 2026-09-05", NOW, (notes + ("; mtime_utc=" + mt if mt else "")).strip("; ")])
for m in models:
    s = m["cnv_sample_id"]; cl = m["cell_line"]; base = re.sub(r"-(R2?)$", r"(\1)", cl)
    big(cl, s, m["recorded_wes_passage"], "provider delivery (McGill; copied to project 2025-09-23)", "provider aligned BAM (GRCh37; Sarek input via samplesheet)", f"{R00}/bam/{base}_{m['recorded_wes_passage']}.bam", notes="Sarek converted to FASTQ and re-aligned to GRCh38")
    big(cl, s, m["recorded_wes_passage"], "provider delivery", "provider BAM index", f"{R00}/bam/{base}_{m['recorded_wes_passage']}.bam.bai")
    big(cl, s, m["recorded_wes_passage"], RUN_VCF + " (overwritten by later runs)", "duplicate-marked CRAM", f"{RES}/preprocessing/markduplicates/{s}/{s}.md.cram", notes="mtime shows producing run")
    big(cl, s, m["recorded_wes_passage"], "Oct-24 or Oct-31 re-execution (see mtime)", "recalibrated CRAM (interval-restricted)", f"{RES}/preprocessing/recalibrated/{s}/{s}.recal.cram", notes="VCF headers name this filename as input, but the retained bytes are from a later re-execution; samtools quickcheck OK 2026-09-05")
    big(cl, s, m["recorded_wes_passage"], "", "recalibrated CRAM index", f"{RES}/preprocessing/recalibrated/{s}/{s}.recal.cram.crai")
    wb = re.search(rf"(\.\./work/\S+/{s}\.recal\.bam)", open(REPO + "/docs/cluster/evidence/cnvkit_commands.txt").read()).group(1)
    big(cl, s, m["recorded_wes_passage"], RUN_OCT24 + " CNVKIT_BATCH task", "BAM converted from Oct-24 recal CRAM; input to manual CNVkit", SRC + "/" + wb[3:], notes="path resolved relative to " + SRC + "/cnvkit (cwd of run_cnvkit.sh)")
    big(cl, s, m["recorded_wes_passage"], RUN_VCF, "mosdepth per-base depth (md)", f"{RES}/reports/mosdepth/{s}/{s}.md.per-base.bed.gz", notes="excluded from copy (17 GB total across stages)")
    big(cl, s, m["recorded_wes_passage"], RUN_VCF, "mosdepth per-base depth (recal)", f"{RES}/reports/mosdepth/{s}/{s}.recal.per-base.bed.gz")
    big(cl, s, m["recorded_wes_passage"], RUN_OCT24, "Sarek-internal CNVkit coverage (flat reference; superseded)", f"{RES}/variant_calling/cnvkit/{s}/{s}.targetcoverage.cnn")
for n_ in NORM:
    d = f"{SRC}/cnvkit/normal_samples/{n_}"
    for f, role in [(f"{n_}.sra", "SRA download"), (f"{n_}_1.fastq", "FASTQ R1 (uncompressed)"), (f"{n_}_2.fastq", "FASTQ R2 (uncompressed)"), (f"{n_}.sam", "bwa-mem2 SAM"), (f"{n_}.bam", "unsorted BAM"),
                    (f"{n_}.sorted.bam", "sorted BAM (CNV reference input)"), (f"{n_}.sorted.bam.bai", "BAM index")]:
        big(n_, n_, "", "cnvkit/normal_samples/run_bwa.sh (Slurm arrays 4401780, 4439216, 4459055; Nov 2025)", role, f"{d}/{f}", notes="public exome PRJNA339046" if f.endswith(".sra") else "")
    big(n_, n_, "", "", "symlink used by cnvkit -n", f"{SRC}/cnvkit/normal_samples/all_bams/{n_}.sorted.bam", notes="symlink -> ../%s/%s.sorted.bam" % (n_, n_))
big("", "", "", "created 2026-05-06 by asmab", "plain (uncompressed) tar of /scratch/asmab/ovcan, misnamed .tar.gz; TRUNCATED", f"{SRC}/ovcan.tar.gz",
    status_note="tar reports 'Unexpected EOF' after 5,663 members (last member ovcan/work/01/5ffbd9bfc6e783c2150fd8d41788a3/OV1369_R2_P66-1_1.merged.fastq.gz cut off)", notes="index saved to wes/inventory/ovcan_tar_gz_index_2026-09-05.txt; small members extracted to wes/scratch_archive_extract")
big("", "", "", "", "original scratch root named in run_sarek.sh and CNVkit commands", "/scratch/asmab/ovcan", notes="mode 700 owned by asmab; contents partially preserved in ovcan.tar.gz")
# in-tar members of interest
idx = [l.rstrip("\n").split(None, 5) for l in open(W + "/inventory/ovcan_tar_gz_index_2026-09-05.txt")]
for f in idx:
    if len(f) < 6: continue
    name = f[5]; sz = f[2]
    if re.search(r"references/.*/(WholeGenomeFasta/Homo_sapiens_assembly38\.fasta|BWAIndex/.*\.64\.bwt|GATKBundle/(1000g_pon|af-only-gnomad|dbsnp_146|Mills_and_1000G|Homo_sapiens_assembly38\.known_indels)[^/]*vcf\.gz)$", name) \
       or re.search(r"wes_new/preprocessing/recalibrated/.*\.recal\.cram$", name):
        rid += 1
        inv.append([f"T{rid:05d}", model_of(name), "", "", "", "reference resource" if "references" in name else "superseded wes_new recal CRAM (Oct 2025 earlier runs)", SRC + "/ovcan.tar.gz::" + name, "", "inside truncated tar; not extracted (listed only)", sz, "", "", "no checksum in searched roots", "", "tar tvf index 2026-09-05", NOW, "tar mtime " + f[3] + " " + f[4]])
tsv(OUT + "/file_inventory.tsv", ["record_id", "model_or_reference", "sample_id", "passage", "run_id", "role", "original_path", "resolved_path", "status", "bytes", "checksum_algorithm", "checksum", "checksum_evidence", "copied_relative_path", "source_evidence", "observed_at", "notes"], inv)
tsv(W + "/inventory/large_files_not_copied.tsv", ["record_id", "model_or_reference", "sample_id", "passage", "run_id", "role", "original_path", "resolved_path", "status", "bytes", "checksum_algorithm", "checksum", "checksum_evidence", "copied_relative_path", "source_evidence", "observed_at", "notes"], [r for r in inv if r[8] != "copied"])
print("done")
