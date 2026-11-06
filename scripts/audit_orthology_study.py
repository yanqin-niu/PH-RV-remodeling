from pathlib import Path
import csv,json,hashlib,collections,math,gzip,os
# Resolve inputs relative to this script, independent of the shell working directory.
ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
import openpyxl
out=Path('outputs/Orthology_study_audit');out.mkdir(exist_ok=True)
def read_sheet(path,sn):
 w=openpyxl.load_workbook(path,read_only=True,data_only=True);s=w[sn];s.reset_dimensions();it=iter(s.values);h=next(it);result = [dict(zip(h,r)) for r in it]; w.close(); return result
def save(name,rs):
 if not rs:return
 with (out/name).open('w',newline='',encoding='utf-8-sig') as f:z=csv.DictWriter(f,fieldnames=list(rs[0]));z.writeheader();z.writerows(rs)
p=Path('data/reference/ORTHOLOGS_RAT.txt')
with p.open(encoding='utf-8') as f: reader=csv.reader((l for l in f if not l.startswith('#')),delimiter='\t'); header=next(reader); ort=[dict(RAT_GENE_SYMBOL=r[0], HUMAN_GENE_SYMBOL=r[4]) for r in reader if len(r)>=5]
rat=collections.defaultdict(set);human=collections.defaultdict(set)
for r in ort:
 a=r['RAT_GENE_SYMBOL'].strip().upper();b=r['HUMAN_GENE_SYMBOL'].strip().upper()
 if a and b:
  rat[a].add(b);human[b].add(a)
base='outputs/Fig5_outputs/Figure5_source_data.xlsx'
cons=read_sheet(base,'Consensus_genes');anchors=read_sheet(base,'Anchor_GSE266139_RV');shared=read_sheet(base,'Shared_input_genes')
shared_set={str(x['Gene']).strip().upper() for x in shared}
anchors=[r for r in anchors if r['Gene'] in shared_set and r['Significant']]
display=read_sheet(base,'Native_results')
keys=[(r['Dataset'],r['Gene']) for r in display]
duplicate_counts=collections.Counter(keys)
save('duplicate_gene_dataset.csv',[{'Dataset':d,'Gene':g,'n_rows':n} for (d,g),n in duplicate_counts.items() if n>1])
# Exact duplicate statistics can be deduplicated; conflicting candidate results must stop.
relevant={a['Gene'] for a in anchors}
relevant.update(t for g in list(relevant) for t in rat.get(g,set()))
groups=collections.defaultdict(list)
for r in display: groups[(r['Dataset'],r['Gene'])].append(r)
for (ds,g),rs in groups.items():
 if len(rs)>1 and g in relevant:
  stats={(r['log2FC'],r['padj']) for r in rs}
  if len(stats)>1: raise ValueError('Conflicting duplicate candidate results: '+ds+' / '+g)

lookup={(r['Dataset'],r['Gene']):r for r in display};datasets=sorted({r['Dataset'] for r in display});humans=[d for d in datasets if 'Human' in d];batches=[d for d in datasets if '240923' in d]
assert len(datasets)==8 and len(humans)==2 and len(batches)==2, 'Unexpected dataset structure'
final={r['Gene'] for r in cons};audit=[];support=[];long=[]
for a in anchors:
 g=a['Gene'];targets=rat.get(g,set());status='one_to_one' if len(targets)==1 and len(human[next(iter(targets))])==1 else 'unmapped' if not targets else 'ambiguous'
 target=next(iter(targets)) if status=='one_to_one' else ''
 votes={};available={}
 for ds in datasets:
  key=target if 'Human' in ds else g
  r=lookup.get((ds,key)) if key else None
  available[ds]=bool(r and isinstance(r['log2FC'],(int,float)) and math.isfinite(r['log2FC']) and isinstance(r['padj'],(int,float)) and math.isfinite(r['padj']))
  votes[ds]=bool(available[ds] and 0<=r['padj']<0.05 and abs(r['log2FC'])>=1 and (r['log2FC']>0)==(a['log2FC']>0))
  long.append({'Gene':g,'Dataset':ds,'query_symbol':key,'mapping_status':status,'discovery_log2FC':a['log2FC'],'log2FC':r['log2FC'] if r else None,'padj':r['padj'] if r else None,'test_available':available[ds],'same_direction_significant':votes[ds]})
 human_n=sum(votes[d] for d in humans);total=sum(votes.values());others=sum(v for d,v in votes.items() if d not in batches)
 audit.append({'Gene':g,'human_ortholog':target,'mapping_status':status,'all_targets':'|'.join(sorted(targets)),'uppercase_matches_ortholog':target==g,'in_primary_signature':g in final})
 row={'Gene':g,'mapping_status':status,'human_ortholog':target,'in_primary_signature':g in final,'n_human_tested':sum(available[d] for d in humans),'n_human_same':human_n,'n_all8_same':total,'mapped_6of8':total>=6,'study_any_6of7':others+any(votes[d] for d in batches)>=6,'study_both_6of7':others+all(votes[d] for d in batches)>=6,'study_any_5of7':others+any(votes[d] for d in batches)>=5,'omit_both_MCT_5of6':others>=5}
 row.update({d+'_same':votes[d] for d in datasets});support.append(row)
save('dataset_evidence.csv',long);save('orthology_audit_candidates.csv',audit);save('human_and_study_support.csv',support)
summary={'n_anchor_candidates':len(anchors),'n_primary_signature':len(final),'mapping_status_candidates':dict(collections.Counter(x['mapping_status'] for x in audit)),'mapping_status_primary_signature':dict(collections.Counter(x['mapping_status'] for x in audit if x['in_primary_signature'])),'primary_signature_human_support_counts':dict(collections.Counter(x['n_human_same'] for x in support if x['in_primary_signature']))}
for rule in ['mapped_6of8','study_any_6of7','study_both_6of7','study_any_5of7','omit_both_MCT_5of6']:
 selected={x['Gene'] for x in support if x[rule]};summary[rule]={'n':len(selected),'retained_original':len(final&selected),'lost_original':sorted(final-selected),'added':sorted(selected-final)}
w=openpyxl.load_workbook('data/raw/GSE240923_processed-data-mct.xlsx',read_only=True,data_only=True);lists=[]
for s in w.worksheets[1:]:s.reset_dimensions();lists.append(list(next(s.values))[1:])
summary['batch_sample_counts']=[len(x) for x in lists];summary['batch_shared_column_ids']=sorted(set(lists[0])&set(lists[1]));summary['orthology_file_sha256']=hashlib.sha256(p.read_bytes()).hexdigest();summary['source']='https://download.rgd.mcw.edu/pub/data_release/RAT/ORTHOLOGS_RAT.txt'
(out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))

# Explicit study units are accession-level except for the merged GSE240923 batches.
# This does not merge GSE242014/GSE240921 from the same publication or weight models equally.
rules=[]
for key in ['mapped_6of8','study_any_6of7','study_both_6of7','study_any_5of7','omit_both_MCT_5of6']:
 z=summary[key];rules.append({'Rule':key,'Selected':z['n'],'Retained_original':z['retained_original'],'Lost_original':'; '.join(z['lost_original']),'Added':'; '.join(z['added'])})
save('rule_summary.csv',rules)
save('primary_signature_human_support.csv',[r for r in support if r['in_primary_signature']])
save('human_supported_subset.csv',[r for r in support if r['in_primary_signature'] and r['n_human_same']>=1])
# Parse actual GEO batch annotations and compare against workbook columns.
series=Path('data/reference/GSE240923_series_matrix.txt.gz')
with gzip.open(series,'rt',encoding='utf-8') as f:
 metadata=[next(csv.reader([l],delimiter='\t')) for l in f if l.startswith('!Sample_')]
accessions=next(r[1:] for r in metadata if r[0]=='!Sample_geo_accession')
titles=next(r[1:] for r in metadata if r[0]=='!Sample_title')
batch=next(r[1:] for r in metadata if r[0]=='!Sample_characteristics_ch1' and r[1].startswith('batch:'))
sex=next(r[1:] for r in metadata if r[0]=='!Sample_characteristics_ch1' and r[1].startswith('Sex:'))
cols=next(r[1:] for r in metadata if r[0]=='!Sample_description' and r[1].startswith('count matrix column name:'))
samples=[]
for gsm,title,b,sx,col in zip(accessions,titles,batch,sex,cols):
 name=col.split(': ',1)[1]; bn=int(b.rsplit(' ',1)[1]);assert name in lists[bn-1]
 samples.append({'GSM':gsm,'Count_column':name,'Batch':bn,'Sex':sx.split(': ',1)[1],'Title':title,'Center':'Not specified in series matrix','Animal_overlap':'No shared sample IDs; alias-level animal reuse not established'})
assert len(set(accessions))==53 and len(samples)==53
save('GSE240923_sample_manifest.csv',samples)
meta=[{'Item':'Reference URL','Value':summary['source']},{'Item':'Orthology SHA256','Value':summary['orthology_file_sha256']},{'Item':'Input Figure5 SHA256','Value':hashlib.sha256(Path(base).read_bytes()).hexdigest()},{'Item':'GEO matrix SHA256','Value':hashlib.sha256(series.read_bytes()).hexdigest()},{'Item':'Mapping policy','Value':'Reciprocal uniqueness within RGD pairwise one-to-one export; ambiguous and unmapped human queries unavailable'},{'Item':'Threshold','Value':'Finite effect and adjusted P; padj >=0 and <0.05; absolute log2FC >=1; anchor direction required'},{'Item':'Missingness','Value':'Missing human mapping or unestimable test is unavailable, not tested nonsignificant; no vote, fixed denominator for selection'},{'Item':'Weighting scope','Value':'Collapse GSE240923 only; not publication-level or model-level weighting across all accessions'},{'Item':'Scope','Value':'Shared-anchor selection sensitivity; not a refit of Figure4 or raw-count models'}]
save('provenance.csv',meta)

# XLSX export uses openpyxl when the optional artifact-tool authoring package is unavailable.
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
wb=openpyxl.Workbook();wb.remove(wb.active)
sheets={'Overview':'rule_summary.csv','Human_supported':'human_supported_subset.csv','Primary_signature':'primary_signature_human_support.csv','Orthology_candidates':'orthology_audit_candidates.csv','Candidate_sensitivity':'human_and_study_support.csv','Dataset_evidence':'dataset_evidence.csv','MCT_sample_manifest':'GSE240923_sample_manifest.csv','Provenance':'provenance.csv'}
for name,filename in sheets.items():
 ws=wb.create_sheet(name)
 with (out/filename).open(encoding='utf-8-sig') as f:
  for i,row in enumerate(csv.reader(f)):
   values=[]
   for v in row:
    if i and v in ('True','False'):v=v=='True'
    elif i:
     try:v=float(v) if '.' in v or 'e-' in v else int(v)
     except ValueError:pass
    values.append(v)
   ws.append(values)
 ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
 for cell in ws[1]:cell.font=Font(name='Calibri',bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='17365D')
 for row in ws:
  for c in row:
   c.alignment=Alignment(vertical='center',wrap_text=True)
   c.border=Border(bottom=Side(style='hair',color='D9D9D9'))
   if c.row>1:c.font=Font(name='Calibri',size=11)
   if c.row>1 and c.row%2==0:c.fill=PatternFill('solid',fgColor='EDF3F8')
 for col in range(1,ws.max_column+1):ws.column_dimensions[get_column_letter(col)].width=23
 ws.row_dimensions[1].height=65
 for i in range(2,ws.max_row+1):ws.row_dimensions[i].height=45
 if name=='Overview':
  ws.column_dimensions['A'].width=29
  for c in ['E','F']:ws.column_dimensions[c].width=65
  for i in range(2,ws.max_row+1):ws.row_dimensions[i].height=105
 if name=='Provenance':
  ws.column_dimensions['A'].width=30;ws.column_dimensions['B'].width=110
  for i in range(2,ws.max_row+1):ws.row_dimensions[i].height=55
 ws.sheet_properties.pageSetUpPr.fitToPage=True
 ws.page_setup.orientation='landscape';ws.page_setup.paperSize=ws.PAPERSIZE_A3
 ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0
 ws.print_title_rows='1:1'
wb.save(out/'Supplementary_Table_orthology_human_support_study_weighting.xlsx')
print('Saved supplementary workbook')
