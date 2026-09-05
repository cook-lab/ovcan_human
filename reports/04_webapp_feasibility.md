# OvCAN Gene-Expression Explorer — Feasibility Assessment & Prototype

**Prepared:** 2026-07-24 · **Author:** Cook Lab (analyst: Claude) · **Status:** Working prototype built + verified

**Verdict:** **Feasible and cheap.** A fully client-side static app is the right tool. The entire
resource (RNA + protein for the whole panel) fits in a **2.6 MB gzipped payload** embedded in the
page, so there is no server, no database, and no runtime network dependency. A working prototype is
built under `app/` and verified end-to-end in a real browser. Recommended host: **GitHub Pages**
(or Netlify / Cloudflare Pages — all equivalent for this).

---

## 1. What the app does

Type a gene symbol → see its **per-cell-line RNA (bulk RNA-seq, TPM)** and **protein (TMT, log2
relative abundance)** across the panel, as horizontal bars **coloured by histotype**, grouped and
sortable by subtype / RNA / protein / name. Missing data is shown honestly (`n.d.`, "no RNA",
"not in protein data"). Modelled on the branded, self-contained `CLDN6_report.html` aesthetic
(rust/navy, mono eyebrows, masthead).

## 2. Data sizes (the reason this is easy)

| Matrix | Source | Shape | Key | On disk (CSV) |
|---|---|---|---|---|
| RNA | `output/rna_tpm.csv` | 39,568 genes × 31 lines (TPM) | **Ensembl gene ID** (rel 105) | 11 MB |
| Protein | `output/prot_abundance_matrix.csv` | 8,430 proteins × 31 lines (log2) | **gene symbol** | 4.4 MB |

After mapping RNA to symbols and dropping unmapped/junk rows, the app-facing dataset is **28,901 RNA
symbols + 8,427 proteins** across a **32-line union** (31 RNA + 31 protein; `VOA6861` is RNA-only,
`VOA14993` protein-only). This is tiny by web standards — the whole thing lives in memory in the
browser and every lookup is an O(1) object read.

## 3. Recommended architecture — self-contained static app

**A single HTML file + one data file, no server.** Chosen because the data is small, the app is
read-only, and it must be hostable anywhere and Artifact-compatible (no external calls).

- **`app/index.html`** (~23 KB): all HTML + CSS + vanilla JS. No frameworks, no charting library
  (bars are CSS `<div>`s), no web fonts fetched (brand fonts are used if installed, else system-ui).
- **`app/data.js`** (~3.4 MB): the payload as **gzipped-then-base64 JSON**, assigned to a global.
  The page decompresses it in-browser with the native **`DecompressionStream('gzip')`** API — no
  library. Loaded via `<script src="data.js">` so it also works from `file://` (unlike `fetch`,
  which browsers block for local files).
- **`app/build_payload.py`**: regenerates `data.js` from `output/` + `metadata/` (documented, ~1 s).

**Payload-size estimate (measured, not guessed):**

| Representation | Size |
|---|---|
| JSON (minified) | 7.6 MB |
| **gzip (what ships over the wire)** | **2.6 MB** |
| base64 (as stored in `data.js` on disk) | 3.4 MB |
| `index.html` | 23 KB |

The payload includes, per gene, the 31 values for each assay plus: the chosen Ensembl gene ID
(traceability), non-protein-coding biotype tags, per-line subtype/site/plex metadata, and global
scale domains. Static hosts serve `data.js` with gzip/brotli transfer encoding, so the **effective
download is ~2.6 MB (less with brotli)** even though the file is 3.4 MB on disk.

*Why not a columnar/binary format (Parquet/Arrow, typed arrays)?* It would shrink the payload
further, but at 2.6 MB gzipped there is no need — JSON keeps the generator and app trivially simple
and debuggable. Noted as a future option only if the panel grows by an order of magnitude.

## 4. Hosting options compared

| Option | Cost | Fit for this app | Notes |
|---|---|---|---|
| **GitHub Pages** ⭐ | Free | **Excellent** | Natural home if the resource ships in the manuscript's GitHub repo — one versioned source of truth, cited under Code Availability. Serves gzip; 1 GB site / 100 MB file limits are irrelevant here. Public repo (or GitHub Pro for private). |
| **Netlify** | Free tier | **Excellent** | Drag-and-drop or git-connected deploy in seconds; brotli; custom domains; trivial rollbacks. Best "just works" option if not tied to a repo. |
| **Cloudflare Pages** | Free | **Excellent** | Fastest global CDN, brotli, generous free limits, git-connected. Interchangeable with Netlify. |
| **Claude Artifact** | — | **Good for internal preview** | Requires a single self-contained file (CSP blocks external hosts) — would need `data.js` inlined into `index.html` (~3.5 MB single file) or a size-reduced payload (protein-coding + in-both genes). Private-by-default, shareable within the org. Great for quick sharing; not a public resource URL. |
| **shinyapps.io** | Free tier limited | **Not recommended** | This is a static app — a Shiny server buys nothing here, adds a runtime, sleeps on the free tier, and caps usage-hours. Only relevant if server-side R computation were ever required (it isn't). |

**Recommendation:** **GitHub Pages** if the app is committed alongside the code/manuscript (most
natural for a published resource and a Data Descriptor's Code Availability). **Netlify or Cloudflare
Pages** are equally good and slightly easier if you want a standalone deploy. **Artifact** is a fine
zero-effort way to circulate a preview internally (inline `data.js` first). Deployment is the PI's
call — nothing here has been deployed.

## 5. Caveats to disclose (surfaced in the app UI + `app/README.md`)

1. **Gene-ID harmonisation.** RNA is Ensembl-gene-ID-keyed (rel 105) and mapped to HGNC symbols;
   protein is symbol-keyed. ~17% of Ensembl IDs have no symbol (dropped). **1,835 symbols map to
   >1 Ensembl gene** (mostly KIR/LILR immune paralogs) — the app keeps the highest-mean-TPM gene and
   displays its Ensembl ID for traceability. Verified: `CLDN6 = ENSG00000184697`, values match the
   prior CLDN6 report exactly (VOA8771 1072, OV3331 892, OV90 85 TPM).
2. **Not every gene is present in every assay.** 8,212 symbols are in both; some genes are absent
   from the RNA matrix altogether (not quantified — e.g. `MALAT1`); ~215 proteins have no RNA match
   (symbol-synonym differences). The UI states this rather than implying zero.
3. **Protein block-missingness (~9%).** Missing values are structured by **TMT plex** (5 plexes): a
   protein absent from a whole plex is missing for every line in it. Shown as `n.d.` — **absence of
   protein is weaker evidence than presence.** 3 junk protein rows (`NA`/blank symbols) were dropped.
4. **Plex / site / subtype are partially confounded** (e.g. all clear-cell lines come from one site).
   A pattern that tracks subtype may partly reflect batch. Cross-subtype contrasts are hypotheses.
5. **TMT ratio compression + units.** Protein is log2 *relative* abundance, not absolute; TMT
   compresses ratios. RNA (TPM) and protein are different assays/units — **compare by rank/pattern,
   not matching numbers.** Bars default to per-gene-relative scaling (numbers always shown); an
   "Absolute" toggle uses fixed cross-gene domains.
6. **One library/replicate per line** — single measurements, no per-line variance. A few lines carry
   subtype **label conflicts** (`OV3331`, `OV90` adeno-vs-HGS; `TOV112D` EC-vs-dedifferentiated),
   flagged with `‡` in the app.

## 6. What works vs what is stubbed

**Works now (verified in headless Chrome + node round-trip of the exact embedded bytes):**
- Masthead metadata, ~37k-symbol typeahead (prefix→substring, keyboard nav), quick-example chips.
- Gene render: Ensembl ID + biotype, detection summary, RNA + protein bars coloured by subtype,
  subtype grouping, sort (subtype/RNA/protein/A–Z), per-gene vs absolute scale toggle.
- Honest states: `n.d.`, "no RNA"/"no protein" per line, gene-absent-from-assay flags, not-found,
  case-insensitive lookup, `?gene=SYMBOL` deep links.
- End-to-end verified: `CLDN6`, `MSLN`, `PAX8`, `WT1` (values match source); `GATD3B` (all-NA
  protein → n.d.); `TARSL2` (protein-only); `ZZZFAKE` (not found).

**Deliberately deferred (listed as "next", not built — feasibility-first):**
- Multi-gene compare / small panel view; CSV/PNG download of the current view.
- **RNA-vs-protein scatter** per gene (concordance) — data already in payload.
- Clickable subtype legend to filter lines; WES/CNV or gene-set/signature overlays.
- For Artifact hosting: inline `data.js` into a single file, or ship a reduced (protein-coding /
  in-both) payload.

## 7. Files

```
app/
├── index.html        # the app (open locally; no server)
├── data.js           # generated payload (gzip+base64 JSON)
├── build_payload.py  # regenerates data.js from output/ + metadata/
└── README.md         # preview + regeneration instructions
reports/
└── 04_webapp_feasibility.md   # this document
```

Preview: `open app/index.html`. Regenerate data: `python3 app/build_payload.py`.
