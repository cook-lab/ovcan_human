#!/usr/bin/env python3
"""
build_report.py — render the canonical results markdown into the branded
Cook Lab HTML (and then PDF via headless Chrome).

Pipeline: reuse the proven Cook Lab <style> block, author a designed hero
(masthead + stat cards + key-findings callout), convert the markdown body
(F1 onward) with python-markdown, wrap tables for horizontal scroll, and
auto-insert each figure after the heading that first references its `f_*`
token. The markdown at reports/01_multiomic_characterization_results.md is the
single source of truth; this script is a pure presentation build.

Usage:  python3 scripts/build_report.py
Then :  (Chrome step is run separately by the caller for the PDF)
"""
import os, re, html, markdown

PROJ = os.environ.get("OVCAN_PROJ", "/Users/dpcook/Analysis/ovcan_human")
MD   = os.path.join(PROJ, "reports", "01_multiomic_characterization_results.md")
SRC_HTML = os.path.join(PROJ, "reports", "01_multiomic_characterization_results.html")
OUT_HTML = SRC_HTML  # redeploy in place
ASSETS = os.path.join(PROJ, "reports", "assets")

# ---- figure captions (Figure N assigned in document order) -----------------
CAPS = {
    "f_rna_qc": "RNA-seq QC — per-line genes detected and pseudoalignment rate (n=31). ~20,000 genes/line at >85% pseudoalignment; the site difference in alignment does not affect detection.",
    "f_prot_bridge": "TMT bridge-channel reproducibility across the 5-plex design — each inter-plex technical replicate vs its original (Pearson ≈0.99, ~7,300 shared proteins).",
    "f_prot_compression": "TMT ratio compression — per-gene cross-line dynamic range (protein vs RNA; protein ≈0.30× RNA), ADC-target ranges, and per-plex structural block-missingness.",
    "f_rna_pca_subtype": "RNA PCA coloured by subtype — HGSC, clear cell, and rare subtypes occupy distinct regions; TOV112D sits with SCCOHT.",
    "f_rna_pca_site": "The same PCA coloured by source site — lines do not segregate by lab; Mes-Masson and Huntsman clear-cell lines intermix (biology, not batch).",
    "f_variance_partition": "Genome-wide variance decomposition (subtype vs source site vs patient) for RNA and protein — subtype ≥ site; site's independent contribution is ~0.",
    "f_rna_markers": "Canonical subtype-marker expression (z-scored VST), lines grouped by subtype — lineage markers elevated in the expected subtype (read by program, not single-gene argmax).",
    "f_concordance": "RNA–protein concordance across 30 lines (8,212 shared genes) vs external benchmarks — per-gene median 0.40, a TMT compression-limited ceiling matching CPTAC/CCLE/ProCan.",
    "f_wes_oncoplot": "Canonical-driver oncoplot — columns are cell lines grouped by a patient-family track; a TMB top-track shows coding-candidate load (TOV21G's ~1,416-variant spike marks the MSI-high line); fill = somatic-confidence tier (Tier 3 excluded from headline frequencies); gene bars count patients (n=16). TP53 universal in HGSC; the four-line 3133 family reads as one patient.",
    "f_wes_cnv": "Autosome copy-number landscape (chrX excluded as a pooled-normal sex artifact) with patient-family and subtype tracks — canonical HGSC events recover; TOV81D (LGSC) quietest.",
    "f_wes_hypermutation": "TOV21G hypermutation / MMR-MSI — ranked coding load (indel-fraction fill), MMR/POLE gene status, GRCh38 SBS-96 spectrum, and COSMIC cosine (SBS6/44/15, not POLE).",
    "f_external_concordance": "External identity validation — RNA concordance vs DepMap for the 5 overlapping lines; each self-matches its namesake at rank 1 of 67 ovarian lines (Spearman 0.74–0.88).",
    "f_auth_swisnf": "Multi-omic SWI/SNF panel (SMARCA4/SMARCA2 across RNA, protein, WES) — TOV112D, COV434, BIN67 deficient; BIN67 shows retained mRNA but lost protein (post-transcriptional).",
    "f_auth_mucinous": "Mucinous authenticity — ovarian-vs-GI discriminators; TOV2414 reads ovarian (KRT7+/PAX8+/SATB2−), VOA8762/VOA8771 read intestinal.",
    "f_adc": "Subtype-resolved ADC-target expression atlas (RNA + protein) — MSLN→HGS, HER2→CC/MC recover; FOLR1 bimodal within HGSC. Lead shortlists with RNA (protein compressed 3–5×).",
    "f_hgs_het": "Within-HGSC pathway-activity strata (Hallmark singscore + PROGENy) — three descriptive subgroups; same-patient families mostly co-cluster. Model-selection example only.",
    "f_passage_check": "Passage-sensitivity — passage vs top RNA PCs; passage is 83% collinear with site and not an independent structure driver.",
    "f_passage_check_crossassay": "Within-line RNA-vs-WES passage mismatch (TOV112D Δ20, OV3331 Δ17; median |Δ|=4).",
    "f_prot_pca": "Proteomic PCA by subtype — HGS carries the separation; weaker than RNA, with a real-but-secondary plex/site component.",
}

# ---- hero (designed; revised headline numbers) -----------------------------
HERO = """
<header class="masthead">
  <p class="eyebrow">Cook Lab · Research / Ovarian Cancer</p>
  <h1>A multi-omic resource of human ovarian cancer cell models <em>recapitulates known subtype biology and corroborates its own identity</em></h1>
  <p class="lede">Uniform bulk RNA-seq, TMT proteomics, and whole-exome sequencing on the OvCAN cell-line panel — 42 models from 34 independent patients. The data are technically sound, reproduce canonical subtype biology, and — read as three layers together and against external references (CCLE/DepMap, Cellosaurus) — corroborate cell-line identity and flag discordances. Prepared as the results backbone for a <em>Scientific Data</em> Data Descriptor.</p>
  <div class="meta">
    <span><b>Prepared</b> 2026-07-24 (post-review revision)</span>
    <span><b>Assays</b> RNA-seq · TMT proteomics · WES</span>
    <span><b>Cohort</b> 42 models · 34 patients · 13 tri-omic</span>
  </div>
</header>

<section id="summary">
  <p class="section-label">Bottom line</p>
  <h2>What we found</h2>
  <div class="stats">
    <div class="stat"><div class="n">42 / 34</div><div class="l">models / independent patients<br>RNA 31 · protein 31 · WES 23</div></div>
    <div class="stat"><div class="n">100%</div><div class="l">TP53 in HGSC<br>11 / 11 patients · positive control</div></div>
    <div class="stat"><div class="n">0.40</div><div class="l">RNA–protein concordance<br>per-gene median · TMT-limited ceiling</div></div>
    <div class="stat"><div class="n">16 / 22</div><div class="l">canonical markers land<br>median marker AUC 0.69</div></div>
  </div>
  <div class="callout">
    <h3>Key findings</h3>
    <ul class="kf">
      <li><b>The data recapitulate known subtype biology — and it is biology, not batch.</b> On a joint variance model, source site adds ≤0.2% of PC1 variance beyond subtype; clear cell (the only subtype at both labs) shows 4–6% site effect. Markers (median AUC 0.69), GO terms, and RNA–protein concordance (per-gene ρ 0.40, a TMT compression ceiling matching CPTAC/CCLE) all recover.</li>
      <li><b>WES recovers canonical genetics at the patient level.</b> TP53 is mutated in 11/11 HGSC patients. Collapsing same-patient sublines corrects the apparent CDK12 enrichment (35% of lines → 18% of patients) and the textbook CNV landscape is robust to collapse (3q gain 82%, 17p loss 82%). Somatic-confidence tiers replace the uninformative germline-VAF flag; defensible somatic BRCA1/2 = 0.</li>
      <li><b>The resource corroborates its own identity, with external support.</b> The three layers recover two published reclassifications (COV434 & BIN67 → SCCOHT — now backed by DepMap SMARCA4-damaging calls and a Cellosaurus misclassification flag; TOV112D → dedifferentiated), and the five lines in DepMap each self-match at rank 1 of 67 ovarian lines. Two "mucinous" lines (VOA8762/VOA8771) look intestinal, not ovarian.</li>
      <li><b>A new reuse feature.</b> TOV21G is a candidate MSI-high / MMR-deficient clear-cell model — 6.9× the panel's mutation load, indel-rich, with a defective-MMR SBS-6/44/15 signature and no POLE lesion.</li>
      <li><b>Stated where the evidence stops.</b> WES is tumor-only (drivers tiered, germline caveat), genomic HRD is not computable, and TMT ratio compression / per-plex missingness are quantified and carried into interpretation. STR/mycoplasma were not done in-house; 30/42 lines have citable originator STR records.</li>
    </ul>
  </div>
  <p><b>How far the evidence goes:</b> technical-validation, patient-level genomics, and identity claims are directly measured, independently reproduced across the parallel re-analysis, and externally cross-checked where public data exist. Small rare-subtype groups (EC/MMMT/SCCOHT n=2; MC n=3) are reported descriptively, never as powered comparisons. <span class="evi med">Evidence: High for validation / identity / patient-level genomics; descriptive for rare-subtype biology</span></p>
</section>
"""

# print-only override: let wide tables WRAP instead of clipping at the page edge
# (the base sheet uses nowrap + horizontal scroll, which has no scroll in print)
EXTRA_CSS = ("<style>@media print{ table th,table td{white-space:normal} "
             ".tw{overflow:visible} }</style>")

def build():
    # CSS from the proven existing HTML head
    src = open(SRC_HTML, encoding="utf-8").read()
    css = re.search(r"<style>.*?</style>", src, re.S).group(0)

    # markdown body: from the first "## F1" onward (hero replaces title+summary)
    md = open(MD, encoding="utf-8").read()
    body_md = md[md.index("## F1"):]
    body = markdown.markdown(body_md, extensions=["tables", "fenced_code", "sane_lists", "attr_list"])

    # wrap tables for horizontal scroll
    body = body.replace("<table>", '<div class="tw"><table>').replace("</table>", "</table></div>")

    # auto-insert figures after the heading that first references each f_* token
    used, counter = set(), [0]
    def fig_block(tok):
        counter[0] += 1
        cap = CAPS.get(tok, "")
        return (f'<figure><img src="assets/{tok}.png" alt="{tok}">'
                f'<figcaption><b>Figure {counter[0]}</b> — {html.escape(cap)}</figcaption></figure>')
    def after_heading(m):
        head = m.group(0)
        if "Figure inventory" in head:              # don't fan figures out under the inventory list
            return head
        toks = [t for t in re.findall(r"f_[a-z0-9_]+", head) if t in CAPS and t not in used]
        for t in toks: used.add(t)
        return head + "".join(fig_block(t) for t in toks)
    body = re.sub(r"<h[23][^>]*>.*?</h[23]>", after_heading, body, flags=re.S)

    # supplementary gallery: any captioned figure not yet placed
    leftover = [t for t in CAPS if t not in used and os.path.exists(os.path.join(ASSETS, t + ".png"))]
    if leftover:
        gallery = '<section><p class="section-label">Supplementary</p><h2>Supplementary figures</h2>'
        gallery += "".join(fig_block(t) for t in leftover) + "</section>"
        body += gallery

    doc = (f'<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n'
           f'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
           f'<title>OvCAN multi-omic characterization — Cook Lab</title>\n{css}\n{EXTRA_CSS}\n</head>\n'
           f'<body>\n<div class="wrap">\n{HERO}\n{body}\n'
           f'<p class="note">Cook Lab · Ottawa Hospital Research Institute · generated from '
           f'<code>reports/01_multiomic_characterization_results.md</code> via '
           f'<code>scripts/build_report.py</code>. Figures in <code>reports/assets/</code>; '
           f'per-line supplement <code>output/supplement_per_line.csv</code>.</p>\n'
           f'</div>\n</body>\n</html>\n')
    open(OUT_HTML, "w", encoding="utf-8").write(doc)
    print(f"Wrote {OUT_HTML}")
    print(f"Figures embedded: {counter[0]} (inline {len(used)}, supplementary {len(leftover)})")

if __name__ == "__main__":
    build()
