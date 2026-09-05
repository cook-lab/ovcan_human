# OvCAN Gene Explorer — static web app

A small, self-contained web app to look up any gene and see its **per-cell-line RNA (TPM)
and protein (TMT log2 abundance)** across the OvCAN ovarian cancer cell-line panel, coloured
by histotype. No server, no build step, no runtime network calls — all data is embedded.

## Files

| File | What it is |
|---|---|
| `index.html` | The entire app (HTML + CSS + vanilla JS, ~23 KB). Open it to run. |
| `data.js` | Auto-generated data payload: a gzipped-then-base64 JSON blob the page gunzips in-browser via `DecompressionStream`. |
| `build_payload.py` | Regenerates `data.js` from the read-only matrices in `../output/` + `../metadata/`. |
| `build_single.py` | Inlines `data.js` into `index.html` → one self-contained file for deployment. |
| `ovcan_viewer_standalone.html` | Generated single-file build; this is the deployment artifact. |

## Preview

Just open the file — no server needed:

```sh
open app/index.html          # macOS
# or drag app/index.html into any modern browser (Chrome/Edge/Firefox/Safari)
```

`index.html` loads `data.js` as a sibling `<script>` (works under `file://`), decompresses it,
and renders. Try `CLDN6`, `MSLN`, `PAX8`, `WT1`. Deep-link a gene with `?gene=SYMBOL`
(e.g. `index.html?gene=MSLN`).

> Requires a browser with `DecompressionStream` (Chrome/Edge ≥80, Firefox ≥113, Safari ≥16.4 —
> i.e. anything from 2023 on). The app shows a clear message if the API is missing.

## Regenerate the embedded data

Run from the **project root** using Python with `pandas` and `numpy`:

```sh
python3 app/build_payload.py
```

Inputs (all read-only):
- `output/rna_tpm.csv` — RNA, Ensembl gene IDs × 31 lines (TPM)
- `output/prot_abundance_matrix.csv` — protein, gene symbols × 31 lines (TMT log2)
- `output/tx2gene_matched.csv` — Ensembl gene ID → symbol map (matched release 93)
- `metadata/samples.csv` — subtype, source site, TMT plex per line

It writes `app/data.js`. The script prints a summary and a CLDN6 sanity check.

## Data harmonisation & caveats (baked into the payload/UI)

- **Gene IDs:** RNA uses the reference whose versioned transcript IDs and lengths match all
  kallisto targets. Genes without a symbol cannot be queried by symbol. Where multiple Ensembl
  genes share a symbol, TPM values are summed, matching the analysis workflows; every contributing
  gene ID is shown. The reference includes alternative loci, recorded in the gene annotations.
- **Coverage:** RNA and protein each cover 31 models; their union is 32 models and their
  intersection is 30. Feature totals are reported by the builder and vary with the reference.
  Secondary protein rows sharing a symbol remain identified as `SYMBOL|UNIPROT`.
- **Protein missingness (~9%)** is block-structured by TMT plex; a protein missing from a whole plex
  is missing for every line in it. Shown as **n.d.** (no quantified value).
- **Confounding:** histotype, source site, and TMT plex are only partly separable. Related models
  are not independent patient replicates. Treat cross-subtype contrasts as hypotheses.
- **Protein scale:** values are supplied log2 abundances normalised using a pooled standard;
  the upstream transformation still requires laboratory confirmation. Observed cross-model
  spread can reflect biology and measurement. The shared display axis does not calibrate
  absolute abundance across different proteins.

The payload stores hashes of both expression matrices, the reference map, and model metadata.
The release builder refuses stale inputs or an outdated standalone build.

## Deploy (cooklab.ca/ovcan_viewer)

`cooklab.ca` is an Astro site (repo `dpcook/cooklab`, auto-deployed via Vercel). Astro serves
everything in `public/` verbatim at the site root, so the app lives at:

    <cooklab>/public/ovcan_viewer/index.html   →   https://cooklab.ca/ovcan_viewer

Deploy the **single-file** build (not the two-file version): it has no external references, so it
is immune to trailing-slash / relative-path issues and needs no server config.

```sh
python3 app/build_payload.py      # only if the data changed
python3 app/build_single.py       # -> app/ovcan_viewer_standalone.html
cp app/ovcan_viewer_standalone.html /Users/dpcook/Projects/cooklab/public/ovcan_viewer/index.html
# then from the cooklab repo: git add public/ovcan_viewer && git commit && git push  (Vercel deploys on push)
```

> Deploying makes the full RNA + protein matrices publicly downloadable (the embedded payload IS
> the dataset). Align going live with the deposition / embargo plan for the descriptor.

See `../reports/04_webapp_feasibility.md` for architecture, hosting options, and the full assessment.
