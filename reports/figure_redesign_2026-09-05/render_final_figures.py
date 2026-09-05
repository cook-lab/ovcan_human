from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,json,hashlib
from PIL import Image, ImageOps, ImageDraw
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'final_pdf_render'; OUT.mkdir(exist_ok=True)
stems=[f'fig{i}' for i in range(1,7)]+[f'figs{i}' for i in range(1,9)]
def render(stem):
    p=ROOT/'docs/manuscript/figures'/f'{stem}.pdf'
    subprocess.run(['/opt/homebrew/bin/pdftoppm','-singlefile','-scale-to','2200','-png',str(p),str(OUT/stem)],check=True,capture_output=True)
    fonts=subprocess.check_output(['/opt/homebrew/bin/pdffonts',str(p)],text=True)
    (OUT/f'{stem}_fonts.txt').write_text(fonts)
    text=subprocess.check_output(['/opt/homebrew/bin/pdftotext','-layout',str(p),'-'],text=True)
    (OUT/f'{stem}.txt').write_text(text)
    rd=PdfReader(p); page=rd.pages[0]
    return {'figure':stem,'pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
            'pages':len(rd.pages),'annotations':len(page.get('/Annots',[])),
            'width_in':float(page.mediabox.width)/72,'height_in':float(page.mediabox.height)/72,
            'fonts':fonts,'extracted_text_chars':len(text)}
with ThreadPoolExecutor(max_workers=4) as pool: rows=list(pool.map(render,stems))
(OUT/'validation.json').write_text(json.dumps(rows,indent=2)+'\n')
for k in range(0,len(stems),6):
    canvas=Image.new('RGB',(1530,1240),'#E2E8F0'); d=ImageDraw.Draw(canvas)
    for j,stem in enumerate(stems[k:k+6]):
        im=Image.open(OUT/f'{stem}.png').convert('RGB'); thumb=ImageOps.contain(im,(500,590))
        x=(j%3)*510; y=(j//3)*620
        canvas.paste('white',(x,y,x+510,y+620)); canvas.paste(thumb,(x+5,y+25))
        d.text((x+8,y+5),stem,fill='#334155')
    canvas.save(OUT/f'contact_{k//6+1}.png')
print('Rendered and extracted all 14 final PDFs.')
