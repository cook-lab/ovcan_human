#!/usr/bin/env python3
"""Package the 14 canonical figures after the R builders have completed.

Requires pypdf. The combined PDF has one figure per page, at its native dimensions;
legends are kept outside figure pages. No analysis data are modified.
"""
from pathlib import Path
import hashlib
import json
import zipfile
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parent
STEMS = [f'fig{i}' for i in range(1, 7)] + [f'figs{i}' for i in range(1, 9)]
LABELS = [f'Figure {i}' for i in range(1, 7)] + [f'Supplementary Figure {i}' for i in range(1, 9)]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    writer = PdfWriter()
    records = []
    for stem, label in zip(STEMS, LABELS):
        pdf = ROOT / f'{stem}.pdf'
        png = ROOT / f'{stem}.png'
        reader = PdfReader(pdf)
        assert len(reader.pages) == 1, f'{pdf}: expected exactly one figure page'
        page = reader.pages[0]
        assert not page.get('/Annots'), f'{pdf}: unexpected PDF annotation objects'
        writer.append(reader, outline_item=label)
        records.append({'figure': label, 'pdf': pdf.name, 'png': png.name,
                        'width_in': float(page.mediabox.width)/72,
                        'height_in': float(page.mediabox.height)/72,
                        'pdf_sha256': sha(pdf), 'png_sha256': sha(png)})
    writer.add_metadata({'/Title': 'OvCAN figures', '/Author': 'Cook Lab',
                         '/Subject': 'Clean figure exports; descriptions in figure_legends.md'})
    bundle = ROOT / 'OvCAN_figures_redesigned.pdf'
    with bundle.open('wb') as f:
        writer.write(f)
    manifest = {'figures': records, 'combined_pdf': {'path': bundle.name, 'sha256': sha(bundle)},
                'legends': {'path': 'figure_legends.md', 'sha256': sha(ROOT/'figure_legends.md')}}
    (ROOT/'figure_export_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    with zipfile.ZipFile(ROOT/'OvCAN_figure_exports.zip', 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for stem in STEMS:
            for ext in ['pdf', 'png']:
                p=ROOT/f'{stem}.{ext}'; z.write(p,p.name)
        for name in ['figure_legends.md', 'README.md', 'figure_export_manifest.json', 'build_figure_bundle.py']:
            z.write(ROOT/name,name)
    print(bundle)
    print(ROOT/'OvCAN_figure_exports.zip')

if __name__ == '__main__':
    main()
