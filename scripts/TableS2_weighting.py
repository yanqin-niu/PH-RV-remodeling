from pathlib import Path
import csv,collections,itertools,json,math,hashlib,openpyxl,re
ROOT=Path(__file__).resolve().parent.parent
out=ROOT/'outputs/TableS2_outputs';out.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(open(p,encoding='utf-8-sig')))
def save(name,rows):
 if not rows:return
 with (out/name).open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
evidence=read(ROOT/'outputs/Orthology_study_audit/dataset_evidence.csv');primary={r['Gene'] for r in read(ROOT/'outputs/Orthology_study_audit/primary_signature_human_support.csv')}
ds=sorted({r['Dataset'] for r in evidence});genes=sorted({r['Gene'] for r in evidence});lookup={(r['Gene'],r['Dataset']):r for r in evidence}
def acc(d):return re.search('GSE\\d+',d)[0]
def publication(d):return 'Boucherat_39196250' if acc(d) in ['GSE240921','GSE240923','GSE242014'] else acc(d)
def model(d):return 'Human' if 'Human' in d else 'MCT' if 'MCT' in d else 'Hypoxia' if '133402' in d else 'SuHx' if '186989' in d else 'PAB'
units=[dict(Dataset=d,Accession=acc(d),Publication_unit=publication(d),Animal_model=model(d),Human_separate=model(d)=='Human') for d in ds];save('unit_definitions.csv',units)
long=[];results=[]
for level,key,subset in [('Cohort',lambda d:d,ds),('Accession',acc,ds),('Publication',publication,ds),('Animal_model',model,[d for d in ds if 'Human' not in d]),('Omit_both_GSE240923',lambda d:d,[d for d in ds if acc(d)!='GSE240923']),('Omit_GSE198618_normalized',lambda d:d,[d for d in ds if acc(d)!='GSE198618'])]:
 group=collections.defaultdict(list)
 for d in subset:group[key(d)].append(d)
 for aggregation in ['any','all'] if level in ['Accession','Publication','Animal_model'] else ['any']:
  by_gene=[]
  for g in genes:
   votes=[];avail=[]
   for unit,members in group.items():
    rr=[lookup[g,d] for d in members];aa=[r['test_available']=='True' for r in rr];vv=[r['same_direction_significant']=='True' for r in rr]
    tested=any(aa) if aggregation=='any' else all(aa);support=any(vv) if aggregation=='any' else all(vv)
    votes.append(support);avail.append(tested)
    long.append(dict(Gene=g,Primary=g in primary,Level=level,Aggregation=aggregation,Unit=unit,Members=';'.join(members),Evaluable=tested,Support=support))
   human=sum(lookup[g,d]['same_direction_significant']=='True' for d in ds if 'Human' in d)
   by_gene.append(dict(Gene=g,Primary=g in primary,Level=level,Aggregation=aggregation,Support=sum(votes),Evaluable=sum(avail),Total_units=len(group),Support_fraction_fixed=sum(votes)/len(group),Support_fraction_evaluable=sum(votes)/sum(avail) if sum(avail) else None,Human_support=human))
  results.extend(by_gene)
save('gene_unit_support.csv',results);save('unit_evidence.csv',long)
summary=[]
for level,agg in sorted({(r['Level'],r['Aggregation']) for r in results}):
 rr=[r for r in results if r['Level']==level and r['Aggregation']==agg];n=rr[0]['Total_units'];thresholds=sorted({math.ceil(.75*n),math.floor(.75*n)})
 for threshold in thresholds:
  ss={r['Gene'] for r in rr if r['Support']>=threshold}
  summary.append(dict(Level=level,Aggregation=agg,Threshold=f'>={threshold}/{n}',Interpretation='Descriptive model breadth, not replacement panel' if level=='Animal_model' else 'Fixed denominator selection sensitivity',Selected=len(ss),Retained_primary=len(ss&primary),Lost=';'.join(sorted(primary-ss)),Added=';'.join(sorted(ss-primary)),Genes=';'.join(sorted(ss))))
save('rule_summary.csv',summary)
# Inventory actual expression-input sample headers; no inference of animal identity from generic labels.
manifest=[]
def header(filename,sheet=0):
 w=openpyxl.load_workbook(ROOT/'data/raw'/filename,read_only=True,data_only=True);s=w.worksheets[sheet] if isinstance(sheet,int) else w[sheet];s.reset_dimensions();r=list(next(s.values));w.close();return r
items=[]
h=header('GSE266139_12967_2025_6792_MOESM2_ESM.xlsx');items.append(('GSE266139',0,[v for v in h if isinstance(v,str) and re.fullmatch('(LP|LM|RP|RM)_\\d+',v)]))
h=header('GSE133402_Rat_RVandLV_HypoxiaNormoxia_Processed.xlsx');items.append(('GSE133402',0,[h[i-2] for i,v in enumerate(h) if v=='counts']))
for accession,filename,sheets in [('GSE240921','GSE240921_processed-data-human.xlsx',[1]),('GSE240923','GSE240923_processed-data-mct.xlsx',[1,2]),('GSE242014','GSE242014_processed-data-pab.xlsx',[1])]:
 for sn in sheets:items.append((accession,sn,header(filename,sn)[1:]))
h=header('GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx');items.append(('GSE186989',0,[v for v in h if isinstance(v,str) and re.match(r'^(CTRL|SU/HX) - .* - Total gene reads$',v)]))
h=next(csv.reader(open(ROOT/'data/raw/GSE198618_Normalized_Counts_RV_ALL.csv')));items.append(('GSE198618',0,h[1:]))
gsmrows=read(ROOT/'outputs/Orthology_study_audit/GSE240923_sample_manifest.csv');gsmmap={r['Count_column']:r for r in gsmrows}
for accession,batch,names in items:
 for name in names:
  gm=gsmmap.get(name,{}) if accession=='GSE240923' else {}
  group='Ctrl' if re.search('CTRL|Control|Normal|RMNT|^[LR]P_',name,re.I) else 'Other_exposure_or_stage'
  manifest.append(dict(Accession=accession,Batch=batch,Sample_column=name,GSM=gm.get('GSM',''),Group_inventory=group,Inventory_scope='All input expression samples; includes LV and intermediate-stage libraries used in model fitting, not all enter disease contrast',Animal_identity='Matching LV/RV prefixes support pairing' if accession in ['GSE133402','GSE266139'] else 'No cross-study animal-level identity crosswalk available'))
save('sample_inventory.csv',manifest)
audit=[]
cohorts=collections.defaultdict(list)
for r in manifest:cohorts[r['Accession']+' / batch '+str(r['Batch'])].append(r)
for a,b in itertools.combinations(cohorts,2):
 aa,bb=cohorts[a],cohorts[b];names=set(r['Sample_column'] for r in aa)&set(r['Sample_column'] for r in bb);gsma={r['GSM'] for r in aa if r['GSM']};gsmb={r['GSM'] for r in bb if r['GSM']};common=gsma&gsmb
 audit.append(dict(Cohort_A=a,Cohort_B=b,Shared_column_labels=';'.join(sorted(names)),Shared_GSM=';'.join(sorted(common)),GSM_comparison_complete=bool(gsma and gsmb),Interpretation='Generic label collisions do not establish reuse; absent label overlap does not exclude aliases. Animal/control reuse not fully verifiable across studies.'))
save('pairwise_overlap.csv',audit)
b1=[r for r in manifest if r['Accession']=='GSE240923' and r['Batch']==1];b2=[r for r in manifest if r['Accession']=='GSE240923' and r['Batch']==2]
assert len(b1)==26 and len(b2)==27 and not ({r['GSM'] for r in b1}&{r['GSM'] for r in b2})
meta=[dict(Item='Input',Value='Current Figure5 results; curated reciprocal-unique RGD mapping for human queries. Existing DESeq2 fits unchanged.'),dict(Item='Selection',Value='Shared-anchor candidates fixed; significance FDR<0.05 and abs(log2FC)>=1 in discovery direction; missing tests unavailable, not nonsignificant. Fixed unit denominators used for selection; evaluable denominators also exported.'),dict(Item='Publication grouping',Value='GSE240921, GSE240923 and GSE242014 conservatively grouped as one publication-level unit; different species/models are NOT asserted to share animals. Other accessions separate.'),dict(Item='Model grouping',Value='Rat MCT, Hypoxia, SuHx and PAB each one unit; human evidence separately counted. Includes discovery MCT; descriptive breadth, not independent model replication.'),dict(Item='Aggregation',Value='any: at least one constituent cohort supports; all: every constituent cohort supports, unavailable constituents cannot provide support. Neither is a pooled-effect meta-analysis.'),dict(Item='Threshold',Value='ceil(75% * units) and floor(75% * units) both shown. With small denominators these differ materially; floor is a relaxed sensitivity, not equal stringency.'),dict(Item='Independence',Value='GSE240923 has 26 and 27 sample records, 53 distinct GSMs and no overlapping count-column IDs, including controls. Distinct donor/animal identity and absence of control reuse under aliases not independently established. No site assignment claimed.'),dict(Item='Audit coverage',Value='All expression input columns inventoried; GEO GSM crosswalk currently complete only for GSE240923. This is a documented identifier audit, not proof that no biological samples overlap across all studies.'),dict(Item='Source hashes',Value=json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'outputs/Fig5_outputs/Figure5_source_data.xlsx',ROOT/'data/reference/GSE240923_series_matrix.txt.gz']}))]
save('methods_and_limits.csv',meta)
wb=openpyxl.Workbook();wb.remove(wb.active)
from openpyxl.styles import Font,PatternFill,Alignment
for title,file in [('Summary','rule_summary.csv'),('Gene_support','gene_unit_support.csv'),('Unit_evidence','unit_evidence.csv'),('Unit_definitions','unit_definitions.csv'),('Samples','sample_inventory.csv'),('Pairwise_overlap','pairwise_overlap.csv'),('Methods_limits','methods_and_limits.csv')]:
 ws=wb.create_sheet(title)
 for r in csv.reader(open(out/file,encoding='utf-8-sig')):ws.append(r)
 ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
 for c in ws[1]:c.font=Font(bold=True,color='FFFFFF');c.fill=PatternFill('solid',fgColor='24476A')
 for col in ws.columns:
  ws.column_dimensions[col[0].column_letter].width=min(60,max(15,len(str(col[0].value))+3))
  for c in col:c.alignment=Alignment(vertical='top',wrap_text=True)
wb.save(out/'TableS2_independence_weighting.xlsx')
print(json.dumps(summary,indent=2));print('Input samples',len(manifest),'GSE240923 controls',sum(r['Group_inventory']=='Ctrl' for r in b1),sum(r['Group_inventory']=='Ctrl' for r in b2))
