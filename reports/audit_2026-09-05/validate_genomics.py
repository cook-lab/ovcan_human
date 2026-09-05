"""Independent consistency checks for the 2026-09-05 genomics corrections."""
import csv, json, math
from pathlib import Path
from collections import Counter, defaultdict
root = Path(__file__).resolve().parents[2]
def rows(path): return list(csv.DictReader(open(root/path)))
variants=rows('output/wes_mutations_filtered.csv'); fam=rows('metadata/line_family_map.csv')
old=rows('reports/audit_2026-09-05/before_genomics/wes_mutations_filtered.csv')
identity=['cell_line','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','HGVSp_Short']
idkey=lambda r:tuple(r[k] for k in identity)
assert {idkey(r) for r in variants if r['cell_line']!='TOV3121D'} == {idkey(r) for r in old}
assert len(variants)==6194 and len(set(r['cell_line'] for r in variants))==23
patients={r['cell_line']:r['patient_id'] for r in fam}
hgsc={r['cell_line'] for r in variants if r['subtype']=='HGS'}
tp53={r['cell_line'] for r in variants if r['subtype']=='HGS' and r['Hugo_Symbol']=='TP53'}
assert hgsc==tp53 and len(hgsc)==18 and len({patients[x] for x in hgsc})==11
cascade=rows('output/wes_filter_cascade.csv')
assert [sum(int(r[k]) for r in cascade) for k in ['raw','pass','rare','coding']]==[582474,16081,15995,6194]
assert all(int(r['raw'])>=int(r['pass'])>=int(r['rare'])>=int(r['coding']) for r in cascade)
assert len(fam)==42 and len(set(patients.values()))==34
assert sum(r['has_wes_maf']=='TRUE' for r in fam)==23
segments=rows('output/wes_cnv_segments.csv');bounds=rows('output/wes_cnv_arm_boundaries.csv'); calls=rows('output/wes_cnv_arm_calls.csv')
assert len(calls)==23*39 and all(r['call']!='not assessed' for r in calls)
by_chrom=defaultdict(list)
for r in segments:by_chrom[r['chromosome']].append(r)
call_lookup={(r['cell_line'],r['arm']):r for r in calls}
# Recompute intersections without using the R arm-intersection output.
for b in bounds:
 if b['assessable']!='TRUE':continue
 length=float(b['arm_bp']); sums=defaultdict(lambda:[0.,0.,0.])
 for s in by_chrom[b['chromosome']]:
  width=max(0,min(float(s['end']),float(b['arm_end']))-max(float(s['start']),float(b['arm_start'])))
  v=float(s['log2c_auto']); a=sums[s['cell_line']];a[2]+=width
  if v>.2:a[0]+=width
  if v<-.2:a[1]+=width
 for line,a in sums.items():
  r=call_lookup[(line,b['arm'])]
  assert all(abs(float(r[k])-v/length)<1e-10 for k,v in zip(['gain_frac','loss_frac','assessed_fraction'],a))
  expected='gain' if a[0]/length>.5 else 'loss' if a[1]/length>.5 else 'neutral'
  assert r['call']==expected
freq=rows('output/wes_cnv_arm_freq_patient.csv')
for f in freq:
 c=[r for r in calls if r['arm']==f['arm'] and r['subtype']=='HGS']
 for direction in ['gain','loss']:
  assert int(f['n_patients_'+direction])==len({r['patient_id'] for r in c if r['call']==direction})
boot=rows('output/wes_signature_refit_bootstrap.csv')
assert len(boot)==23*(60+22) and len({(r['cell_line'],r['signature'],r['reference_set']) for r in boot})==len(boot)
assert all(int(r['n_boots'])==200 and 0<=float(r['boot_selected_frac'])<=1 and r['cache_sha256'] for r in boot)
ref=rows('output/wes_cosmic_target_normalized.csv')
assert len(ref)==96
for k in ref[0]:
 if k!='context':assert abs(sum(float(r[k]) for r in ref)-1)<1e-10
selfmatch=rows('output/external_selfmatch_margin.csv')
assert len(selfmatch)==5 and all(r['self_rank_of67']=='1' and r['reciprocal_best']=='TRUE' for r in selfmatch)
strs=rows('output/cellosaurus_str_status.csv');assert len(strs)==42
assert sum(r['str_profile_documented']=='TRUE' for r in strs)==30
assert all('pending author confirmation' in r['current_stock_str_status'] for r in strs)
tiers=rows('output/wes_driver_tiers.csv')
assert next(r for r in tiers if r['cell_line']=='OV90' and r['gene']=='CDKN2A')['tier']=='Tier3'
result={'status':'passed','variant_models':23,'coding_candidates':6194,'legacy_candidates_preserved':6036,'HGSC_TP53_models':18,'HGSC_TP53_patients':11,'independently_validated_arm_calls':len(calls),'bootstrap_rows':len(boot),'reciprocal_DepMap_matches':5,'Cellosaurus_reference_STR_profiles':30}
(root/'reports/audit_2026-09-05/genomics_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
