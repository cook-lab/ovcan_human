from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sections={}
for name in ['legends_overview_adc.md','legends_genomics.md','legends_supp_expression.md']:
    s=(HERE/name).read_text()
    for m in re.finditer(r'^## ((?:Supplementary )?Figure (?:S)?\d+)\. ([^\n]+)\n+(.*?)(?=^## |\Z)',s,re.M|re.S):
        sections[m[1]]=(m[2],m[3].strip())
s=(HERE/'legends_expression.md').read_text()
for m in re.finditer(r'\*\*(Figure [23]) \| (.*?)\*\* (.*?)(?=^## |\Z)',s,re.M|re.S):
    sections[m[1]]=(m[2].rstrip('.'),m[3].strip())
keys=[f'Figure {i}' for i in range(1,7)]+[f'Supplementary Figure S{i}' for i in range(1,9)]
assert len(sections)==14 and set(sections)==set(keys),(len(sections),set(keys)-set(sections))
text='# OvCAN figure legends\n\nLegends for the September 2026 figure redesign. Figures 1–4 accompany the current Data Descriptor draft; Figures 5–6 and S1–S8 are additional analysis figures. The corresponding PDF pages contain the artwork alone.\n\n'
for k in keys:
    title,body=sections[k]
    if k=='Figure 3':
        body+=' Fifteen of 25 marker checks were recovered in the patient-representative assessment, compared with a mean of 8.46 under 20,000 joint label permutations within centre (P = 0.0101); MKI67 was excluded.'
    text+=f'## {k}. {title}\n\n{body}\n\n'
(ROOT/'docs/manuscript/figures/figure_legends.md').write_text(text)
print('Assembled all 14 legends.')
