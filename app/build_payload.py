#!/usr/bin/env python3
# =============================================================================
# build_payload.py  —  OvCAN gene-expression explorer: embedded-data generator
# -----------------------------------------------------------------------------
# Reads the read-only resource matrices in output/ + metadata, and writes a
# single compact, gzipped+base64 payload (app/data.js) that the static app
# (app/index.html) decompresses in-browser (no server, no network).
#
#   RNA   : output/rna_tpm.csv            (Ensembl gene IDs x 31 lines, TPM)
#   Protein: output/prot_abundance_matrix.csv (gene symbols x 31 lines, log2)
#   Map   : output/tx2gene_matched.csv (Ensembl gene ID -> symbol)
#   Meta  : metadata/samples.csv          (subtype, source_site, tmt_plex ...)
#
# Gene-ID harmonisation matches the analysis: sum TPM over Ensembl gene IDs
# sharing a non-empty symbol, retaining every contributing ID for traceability.
# This includes mapped alternative loci; use deposited gene annotations for
# primary-assembly-only analyses. Protein duplicate rows retain their identifiers.
#
# Run:  python3 app/build_payload.py            (from project root)
#   or: python3 build_payload.py --root /path/to/ovcan_human
# =============================================================================
import argparse, base64, gzip, json, os, sys, datetime, re, hashlib
import pandas as pd
import numpy as np

# a symbol that is really "no symbol": blank, nan/null, or the NA / NA.1 / NA.2 pandas artifacts
_JUNK = re.compile(r"^(nan|null|na(\.\d+)?)$", re.IGNORECASE)
def is_real_symbol(s):
    s = str(s).strip()
    return bool(s) and not _JUNK.match(s)

def sig(x, n=4):
    """Round to n significant figures; keeps small TPMs meaningful, shrinks big ones."""
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return None
    if x == 0:
        return 0
    return float(f"{x:.{n}g}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=None, help="project root (default: parent of this script)")
    ap.add_argument("--out",  default=None, help="output data.js path (default: <root>/app/data.js)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    root = args.root or os.path.dirname(here)
    out  = args.out  or os.path.join(root, "app", "data.js")
    O    = os.path.join(root, "output")
    M    = os.path.join(root, "metadata")

    # ---- 1. Ensembl gene ID -> symbol / biotype ----------------------------
    t2g = pd.read_csv(os.path.join(O, "tx2gene_matched.csv"))
    g2 = t2g.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
    sym_of = g2["external_gene_name"]
    bt_of  = g2["gene_biotype"]

    # ---- 2. RNA (TPM) ------------------------------------------------------
    rna = pd.read_csv(os.path.join(O, "rna_tpm.csv"), index_col=0)
    rna_lines = list(rna.columns)
    symbols = sym_of.reindex(rna.index)
    keep = symbols.map(is_real_symbol)
    rna = rna.loc[keep]
    sym_kept = symbols.loc[keep].astype(str)
    biotype_kept = bt_of.reindex(rna.index)

    # Sum linear TPM before log transformation, as in the RNA/protein analysis.
    rna_sum = rna.groupby(sym_kept, sort=True).sum()
    assert np.allclose(rna_sum.sum(axis=0), rna.sum(axis=0), rtol=1e-10)
    ids_by_symbol = {}
    for gene_id, symbol in sym_kept.items():
        ids_by_symbol.setdefault(symbol, []).append(gene_id)
    rna_dict, rna_gene_id, rna_biotype_nonpc = {}, {}, {}
    for symbol, vals in rna_sum.iterrows():
        ids = ids_by_symbol[symbol]
        biotypes = sorted(set(biotype_kept.loc[ids].dropna()))
        rna_dict[symbol] = [sig(float(v), 4) for v in vals]
        rna_gene_id[symbol] = ", ".join(ids)
        if biotypes != ["protein_coding"]:
            rna_biotype_nonpc[symbol] = ", ".join(biotypes) or "unannotated"

    # ---- 3. Protein (supplied log2 normalised abundance; NaN = unquantified)
    prot = pd.read_csv(os.path.join(O, "prot_abundance_matrix.csv"), index_col=0)
    prot_lines = list(prot.columns)
    prot_dict = {}
    n_prot_dropped = 0
    for symbol, row in zip(prot.index, prot.values):
        if not is_real_symbol(symbol):          # drop 'nan'/'NA'/'NA.1' junk rows
            n_prot_dropped += 1
            continue
        prot_dict[str(symbol)] = [None if pd.isna(v) else round(float(v), 3) for v in row]

    # ---- 4. Per-line metadata (single source of truth = samples.csv) -------
    samp = pd.read_csv(os.path.join(M, "samples.csv"))
    samp = samp.set_index("cell_line")
    all_lines = list(dict.fromkeys(rna_lines + prot_lines))   # union, stable order
    lines_meta = {}
    for ln in all_lines:
        row = samp.loc[ln] if ln in samp.index else None
        subtype = str(row["subtype"]) if row is not None else "NA"
        site    = str(row["source_site"]) if row is not None else "NA"
        plex = None
        if row is not None and pd.notna(row.get("tmt_plex")) and str(row.get("tmt_plex")) not in ("", "nan", "-"):
            try:
                plex = str(int(float(row["tmt_plex"])))       # "4.0" -> "4"
            except (ValueError, TypeError):
                plex = str(row["tmt_plex"])
        # subtype_status flags a few label conflicts (OV3331/OV90 adeno-vs-HGS; TOV112D EC-vs-dediff)
        stat = str(row["subtype_status"]) if row is not None else ""
        conflict = stat.startswith("CONFLICT") or "reassign" in stat.lower() or "conflict" in stat.lower()
        lines_meta[ln] = {"subtype": subtype, "site": site, "plex": plex,
                          "has_rna": ln in rna_lines, "has_prot": ln in prot_lines,
                          "conflict": bool(conflict)}

    # ---- 5. Assemble + compress -------------------------------------------
    payload = {
        "generated": datetime.date.today().isoformat(),
        "source": "OvCAN human ovarian cancer cell-line multi-omic resource",
        "rna_unit": "TPM (kallisto+tximport, summed by gene symbol)",
        "prot_unit": "log2 protein abundance (pooled-standard normalised TMT)",
        "ensembl_release": 93,
        "rna_symbol_aggregation": "sum of TPM across mapped Ensembl gene IDs",
        "source_sha256": {name: hashlib.sha256(open(os.path.join(O, name), "rb").read()).hexdigest()
                          for name in ("rna_tpm.csv", "prot_abundance_matrix.csv", "tx2gene_matched.csv")},
        "metadata_sha256": hashlib.sha256(open(os.path.join(M, "samples.csv"), "rb").read()).hexdigest(),
        "rna_lines": rna_lines,
        "prot_lines": prot_lines,
        "all_lines": all_lines,
        "subtype_order": ["HGS", "LGS", "CC", "EC", "MC", "MMMT", "SCCOHT"],
        "lines": lines_meta,
        "rna": rna_dict,
        "rna_gene_id": rna_gene_id,
        "rna_biotype_nonpc": rna_biotype_nonpc,
        "prot": prot_dict,
    }
    # global scale domains (for the app's optional "absolute" bar scaling) --
    rna_all = np.concatenate([np.asarray(v, dtype=float) for v in rna_dict.values()])
    prot_all = np.array([v for row in prot_dict.values() for v in row if v is not None], dtype=float)
    payload["scale"] = {
        "rna_max": sig(float(np.nanmax(rna_all)), 4),            # top of log10(TPM+1) domain
        "prot_lo": round(float(np.nanpercentile(prot_all, 1)), 3),
        "prot_hi": round(float(np.nanpercentile(prot_all, 99)), 3),
    }

    js_json = json.dumps(payload, separators=(",", ":"))
    gz = gzip.compress(js_json.encode("utf-8"), 9, mtime=0)
    b64 = base64.b64encode(gz).decode("ascii")

    os.makedirs(os.path.dirname(out), exist_ok=True)
    header = ("// AUTO-GENERATED by app/build_payload.py — do not edit by hand.\n"
              "// gzipped+base64 JSON payload; app/index.html gunzips it in-browser.\n")
    with open(out, "w") as f:
        f.write(header)
        f.write('window.OVCAN_B64="' + b64 + '";\n')

    # ---- 6. Report ---------------------------------------------------------
    n_prot_only = len(set(prot_dict) - set(rna_dict))
    print("payload written:", out)
    print(f"  RNA symbols       : {len(rna_dict):,}  (from {len(rna):,} symboled Ensembl genes; TPM summed across duplicate symbols)")
    print(f"  Protein symbols   : {len(prot_dict):,}  ({n_prot_only} protein-only; {n_prot_dropped} junk rows dropped)")
    print(f"  RNA lines         : {len(rna_lines)}   Protein lines: {len(prot_lines)}   Union: {len(all_lines)}")
    print(f"  JSON              : {len(js_json)/1e6:.2f} MB")
    print(f"  gzip              : {len(gz)/1e6:.2f} MB")
    print(f"  base64 (on disk)  : {len(b64)/1e6:.2f} MB")
    # sanity: CLDN6
    if "CLDN6" in rna_dict:
        i = rna_lines.index("OV3331")
        print(f"  sanity CLDN6 RNA OV3331 = {rna_dict['CLDN6'][i]} TPM, Ensembl {rna_gene_id['CLDN6']}")

if __name__ == "__main__":
    main()
