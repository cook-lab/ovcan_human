#!/usr/bin/env python3
"""
build_report_full.py — render the consolidated results markdown into the branded
Cook Lab HTML, then PDF via headless Chrome.

The markdown at reports/01_multiomic_characterization_results.md is the single
source of truth. This script is a pure presentation build: it converts the body
with python-markdown and expands three custom constructs that keep the markdown
readable on its own while still producing a designed page.

  {{stats}} ... {{endstats}}          -> the stat-card row in the masthead block
      one card per line:  value || label (a <br> is inserted at the em dash)

  {{callout: Title}} ... {{endcallout}} -> the bordered key-findings callout

  {{figure: id | caption}}            -> a <figure> with reports/assets2/<id>.svg
      the label is derived from the id, so fig3 -> "Figure 3" and figs7 ->
      "Figure S7", matching the Fig. N references in the prose.

Figures come from reports/assets2/, which holds the manuscript figure set
regenerated with OVCAN_FIG_PLAIN=1 (in-panel methodological footnotes
suppressed; geometry, palettes and statistics identical).

FIGURE FORMAT.  Figures are embedded as SVG by default, converted from the
cairo_pdf originals with `pdftocairo -svg` (regenerated automatically when the
PDF is newer).  Chrome preserves SVG as vector through print-to-PDF, so the
figures stay resolution-independent and the report PDF is ~half the size of the
PNG build (3.4 MB vs 6.4 MB) with the same pagination.  Two consequences worth
knowing: pdftocairo converts glyphs to paths, so figure text is crisp and
font-independent but NOT selectable in the output PDF (body text still is); and
layers deliberately rasterised in R with ggrastr (the dense bridge and
concordance scatters) stay raster inside the SVG, which is correct.  Pass
--figs=png to fall back to the 400 dpi rasters.

Usage:  python3 scripts/build_report_full.py [--pdf] [--figs=svg|png]
"""
import os
import re
import sys
import html
import subprocess
import markdown

PROJ = os.environ.get("OVCAN_PROJ", "/Users/dpcook/Analysis/ovcan_human")
MD = os.path.join(PROJ, "reports", "01_multiomic_characterization_results.md")
OUT_HTML = os.path.join(PROJ, "reports", "01_multiomic_characterization_results.html")
OUT_PDF = os.path.join(PROJ, "reports", "01_multiomic_characterization_results.pdf")
ASSETS = os.path.join(PROJ, "reports", "assets2")
ASSET_REL = "assets2"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
PDFTOCAIRO = "pdftocairo"          # poppler; ships with pdfinfo/pdftoppm


def ensure_svg(fid):
    """Convert reports/assets2/<fid>.pdf -> .svg, refreshing a stale SVG.

    Returns True when an up-to-date SVG exists. pdftocairo is cairo re-emitting
    its own drawing model, so this is a lossless vector transform of the same
    cairo_pdf that produced the manuscript figure.
    """
    pdf = os.path.join(ASSETS, fid + ".pdf")
    svg = os.path.join(ASSETS, fid + ".svg")
    if not os.path.exists(pdf):
        return os.path.exists(svg)
    if os.path.exists(svg) and os.path.getmtime(svg) >= os.path.getmtime(pdf):
        return True
    try:
        subprocess.run([PDFTOCAIRO, "-svg", pdf, svg], check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError) as e:
        print("  pdftocairo failed for %s (%s)" % (fid, e))
        return False
    return True

# --------------------------------------------------------------------------- #
# Style sheet — the proven Cook Lab report sheet (rustNavy tokens), matched to
# scripts/00b_figure_theme.R so the page and the figures read as one system.
# --------------------------------------------------------------------------- #
CSS = """<style>
@import url('https://fonts.googleapis.com/css2?family=Manrope:wght@300;400;500;600;700&family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
:root{
  --rust:#C2410C; --rust-400:#EB6235; --rust-700:#7A2A09; --rust-50:#FEF3EE; --rust-100:#FDE2D1;
  --navy:#0F172A; --slate:#334155; --slate-400:#64748B; --slate-300:#94A3B8; --slate-200:#CBD5E1;
  --teal:#0D9488; --bg:#FFFFFF; --paper:#FFFFFF; --ink:#0F172A;
  --line:rgba(15,23,42,.10); --line-strong:rgba(15,23,42,.16);
  --font-display:'Manrope',system-ui,sans-serif;
  --font-body:'Inter',system-ui,sans-serif;
  --font-mono:'JetBrains Mono',ui-monospace,monospace;
}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--font-body);font-size:14px;line-height:1.55;font-weight:400;-webkit-font-smoothing:antialiased;}
.wrap{max-width:1000px;margin:0 auto;padding:0 26px 72px}
.masthead{padding:44px 0 20px;border-bottom:2px solid var(--rust);margin-bottom:8px}
.eyebrow{font-family:var(--font-mono);font-size:10.5px;font-weight:500;text-transform:uppercase;letter-spacing:.15em;color:var(--rust);margin:0 0 14px}
h1{font-family:var(--font-display);font-weight:300;letter-spacing:-.035em;line-height:1.08;font-size:2.2rem;margin:0 0 12px;color:var(--navy)}
h1 em{font-style:normal;color:var(--rust)}
.lede{font-size:1rem;color:var(--slate);max-width:72ch;margin:0 0 16px;font-weight:400;line-height:1.55}
.meta{font-family:var(--font-mono);font-size:10.5px;text-transform:uppercase;letter-spacing:.1em;color:var(--slate-400);display:flex;gap:18px;flex-wrap:wrap}
.meta b{color:var(--slate);font-weight:500}
section{margin-top:38px}
.section-label{font-family:var(--font-mono);font-size:10.5px;font-weight:500;text-transform:uppercase;letter-spacing:.14em;color:var(--rust);margin:0 0 8px}
h2{font-family:var(--font-display);font-weight:400;letter-spacing:-.025em;font-size:1.42rem;margin:38px 0 12px;color:var(--navy);padding-top:10px;border-top:1px solid var(--line)}
h2:first-of-type{border-top:none;padding-top:0}
h3{font-family:var(--font-display);font-weight:500;font-size:1.04rem;margin:24px 0 6px;color:var(--navy)}
p{margin:0 0 11px;max-width:82ch}
li{max-width:80ch}
a{color:var(--rust);text-decoration:none;border-bottom:1px solid var(--rust-100)}
strong{font-weight:600;color:var(--navy)}
em{font-style:italic}
code{font-family:var(--font-mono);font-size:.84em;background:var(--rust-50);color:var(--rust-700);padding:1px 5px;border-radius:4px}
ul,ol{padding-left:22px;margin:0 0 12px}
ul li,ol li{margin-bottom:6px}
.callout{background:var(--paper);border:1px solid var(--line);border-left:3px solid var(--rust);border-radius:10px;padding:18px 22px;margin:18px 0;box-shadow:0 1px 3px rgba(15,23,42,.06)}
.callout h3{margin-top:0}
.kf{list-style:none;padding:0;margin:0}
.kf li{position:relative;padding:9px 0 9px 26px;border-bottom:1px solid var(--line);max-width:none}
.kf li:last-child{border-bottom:none}
.kf li::before{content:"";position:absolute;left:5px;top:17px;width:6px;height:6px;border-radius:50%;background:var(--rust)}
.kf b{color:var(--navy)}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}
.stat{background:var(--paper);border:1px solid var(--line);border-radius:10px;padding:15px 14px;text-align:left}
.stat .n{font-family:var(--font-display);font-weight:300;font-size:1.85rem;color:var(--rust);line-height:1;letter-spacing:-.03em}
.stat .l{font-family:var(--font-mono);font-size:9.5px;text-transform:uppercase;letter-spacing:.06em;color:var(--slate-400);margin-top:8px;line-height:1.45}
figure{margin:22px 0;background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:14px;box-shadow:0 1px 3px rgba(15,23,42,.05)}
figure img{width:100%;height:auto;display:block;border-radius:4px}
figcaption{font-size:12px;color:var(--slate);margin-top:10px;padding-top:10px;border-top:1px solid var(--line);line-height:1.5;max-width:none}
figcaption .lab{font-family:var(--font-mono);font-size:10.5px;text-transform:uppercase;letter-spacing:.08em;color:var(--rust);font-weight:500}
.tw{overflow-x:auto;margin:16px 0;border:1px solid var(--line);border-radius:10px;background:var(--paper)}
table{border-collapse:collapse;width:100%;font-size:12px}
th,td{padding:6px 11px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
thead th{font-family:var(--font-mono);font-size:9.5px;text-transform:uppercase;letter-spacing:.06em;color:var(--slate-400);font-weight:500;background:var(--rust-50);white-space:nowrap}
tbody tr:last-child td{border-bottom:none}
.note{font-size:12px;color:var(--slate-400);border-top:1px solid var(--line);padding-top:14px;margin-top:32px;max-width:none}
hr{border:none;border-top:1px solid var(--line);margin:0}
@media (max-width:760px){h1{font-size:1.8rem}.stats{grid-template-columns:repeat(2,1fr)}.meta{gap:12px}}
@media print{
  @page{size:letter;margin:10mm}
  body{background:#fff;font-size:10.5px;line-height:1.5}
  .wrap{max-width:none;padding:0}
  figure,.callout,.stat,.tw{box-shadow:none}
  h1{font-size:1.8rem}
  h2{font-size:1.2rem;break-after:avoid;break-before:auto}
  h3{font-size:.98rem;break-after:avoid}
  p,li{max-width:none}
  figure,.callout,.tw,.stat{break-inside:avoid}
  /* The image cap must leave room for the longest caption, or the whole figure
     block exceeds one page, break-inside:avoid cannot place it, and Chrome emits
     a near-blank page. Text height on letter at 10 mm margins is ~259 mm; the
     longest caption here runs ~35 mm, so 205 mm is the hard ceiling and 200 mm the
     working value. Pagination is insensitive to this between 160 and 200 mm, so
     the larger value is kept: with vector figures there is no reason to shrink. */
  figure{padding:8px;margin:16px 0}
  figure img{max-height:200mm;width:auto;max-width:100%;margin:0 auto}
  tr{break-inside:avoid}
  thead{display:table-header-group}
  thead th{white-space:normal}
  table{font-size:9px}
  th,td{padding:4px 7px}
  .stats{break-inside:avoid}
  .masthead{padding-top:0}
  figcaption{font-size:9px;margin-top:8px;padding-top:8px}
}
</style>"""


def fig_label(fid):
    """fig3 -> Figure 3 ; figs7 -> Figure S7 (matches the Fig. N prose refs)."""
    m = re.fullmatch(r"figs(\d+)", fid)
    if m:
        return "Figure S%s" % m.group(1)
    m = re.fullmatch(r"fig(\d+)", fid)
    if m:
        return "Figure %s" % m.group(1)
    return fid


def inline(md_text):
    """Markdown-render a fragment without wrapping it in <p> when it is one line."""
    out = markdown.markdown(md_text.strip(), extensions=["sane_lists"])
    if out.startswith("<p>") and out.endswith("</p>") and out.count("<p>") == 1:
        out = out[3:-4]
    return out


def build(fmt="svg"):
    src = open(MD, encoding="utf-8").read()

    # ---- title -----------------------------------------------------------
    m = re.search(r"^# (.+)$", src, re.M)
    title_md = m.group(1).strip()
    title_html = inline(title_md)          # *— consolidated…* becomes <em>
    plain_title = re.sub(r"[*_]", "", title_md)

    # ---- lede + meta -----------------------------------------------------
    lede = re.search(r"^\*\*Lede:\*\*\s*(.+)$", src, re.M).group(1).strip()
    meta_line = re.search(r"^\*\*Prepared\*\*\s*(.+)$", src, re.M).group(0).strip()
    meta_items = [inline(x.strip()) for x in meta_line.split(" · ")]
    meta_html = "".join("<span>%s</span>" % x for x in meta_items)

    # ---- stat cards ------------------------------------------------------
    stats_block = re.search(r"\{\{stats\}\}(.*?)\{\{endstats\}\}", src, re.S).group(1)
    cards = []
    for line in stats_block.strip().splitlines():
        if "||" not in line:
            continue
        val, lab = [x.strip() for x in line.split("||", 1)]
        lab = lab.replace(" — ", "<br>")
        cards.append('<div class="stat"><div class="n">%s</div><div class="l">%s</div></div>'
                     % (inline(val), inline(lab)))
    stats_html = '<div class="stats">%s</div>' % "".join(cards)

    # ---- callout ---------------------------------------------------------
    cm = re.search(r"\{\{callout:\s*(.+?)\}\}(.*?)\{\{endcallout\}\}", src, re.S)
    co_items = []
    for chunk in re.findall(r"^- (.*(?:\n(?!- ).*)*)", cm.group(2), re.M):
        co_items.append("<li>%s</li>" % inline(chunk))
    callout_html = ('<div class="callout"><h3>%s</h3><ul class="kf">%s</ul></div>'
                    % (html.escape(cm.group(1)), "".join(co_items)))

    # ---- body: everything from the first "## " heading onward ------------
    body_md = src[src.index("\n## "):]

    # pull figures out before markdown conversion so captions can hold markdown
    figs = {}

    def stash(m):
        fid, cap = m.group(1).strip(), m.group(2).strip()
        key = "@@FIG%d@@" % len(figs)
        figs[key] = (fid, cap)
        return key

    body_md = re.sub(r"\{\{figure:\s*([A-Za-z0-9_]+)\s*\|\s*(.*?)\}\}", stash, body_md, flags=re.S)

    body = markdown.markdown(body_md,
                             extensions=["tables", "fenced_code", "sane_lists", "attr_list"])
    body = body.replace("<table>", '<div class="tw"><table>').replace("</table>", "</table></div>")

    missing, fell_back = [], []
    for key, (fid, cap) in figs.items():
        ext = "png"
        if fmt == "svg":
            if ensure_svg(fid):
                ext = "svg"
            else:
                fell_back.append(fid)
        if not os.path.exists(os.path.join(ASSETS, "%s.%s" % (fid, ext))):
            missing.append(fid)
        blk = ('<figure><img src="%s/%s.%s" alt="%s">'
               '<figcaption><span class="lab">%s</span> — %s</figcaption></figure>'
               % (ASSET_REL, fid, ext, fig_label(fid), fig_label(fid), inline(cap)))
        # the placeholder sits alone in a paragraph after conversion
        body = body.replace("<p>%s</p>" % key, blk).replace(key, blk)

    doc = ('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n'
           '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
           '<title>%s — Cook Lab</title>\n%s\n</head>\n<body>\n<div class="wrap">\n'
           '<header class="masthead">\n'
           '  <p class="eyebrow">Cook Lab · Research / Ovarian Cancer</p>\n'
           '  <h1>%s</h1>\n  <p class="lede">%s</p>\n  <div class="meta">%s</div>\n'
           '</header>\n'
           '<section id="summary">\n  <p class="section-label">At a glance</p>\n'
           '  %s\n  %s\n</section>\n'
           '%s\n'
           '<p class="note">Cook Lab · Ottawa Hospital Research Institute · generated from '
           '<code>reports/01_multiomic_characterization_results.md</code> via '
           '<code>scripts/build_report_full.py</code>. Figures are vector, in '
           '<code>reports/assets2/</code> (SVG as embedded, PDF and 400 dpi PNG alongside); '
           'per-model master table <code>output/supplement_per_line.csv</code>.</p>\n'
           '</div>\n</body>\n</html>\n'
           % (html.escape(plain_title), CSS, title_html, inline(lede), meta_html,
              stats_html, callout_html, body))

    open(OUT_HTML, "w", encoding="utf-8").write(doc)
    print("Wrote %s (%.0f KB)" % (OUT_HTML, len(doc) / 1024))
    print("Figures embedded: %d as %s" % (len(figs), fmt))
    if fell_back:
        print("  no SVG available, used PNG for: %s" % ", ".join(fell_back))
    if missing:
        print("MISSING figure files: %s" % ", ".join(missing))
    return len(figs)


def to_pdf():
    if not os.path.exists(CHROME):
        print("Chrome not found; skipping PDF")
        return
    subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
                    "--virtual-time-budget=20000",
                    "--print-to-pdf=%s" % OUT_PDF, "file://%s" % OUT_HTML],
                   check=True, capture_output=True)
    print("Wrote %s (%.1f MB)" % (OUT_PDF, os.path.getsize(OUT_PDF) / 1e6))


if __name__ == "__main__":
    fmt = "svg"
    for a in sys.argv[1:]:
        if a.startswith("--figs="):
            fmt = a.split("=", 1)[1]
    if fmt not in ("svg", "png"):
        sys.exit("--figs must be svg or png")
    build(fmt)
    if "--pdf" in sys.argv:
        to_pdf()
