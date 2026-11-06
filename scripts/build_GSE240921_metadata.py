"""Recover a documented technical stratum and GEO sex; do not infer clinical IDs."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import collections
import openpyxl

project = Path(__file__).resolve().parents[1]
analysis = project
source = project / 'data/raw/GSE240921_processed-data-human.xlsx'
geo = project / 'data/reference/GSE240921_series_matrix.txt.gz'
with gzip.open(geo, 'rt', encoding='utf-8') as stream:
    rows = list(csv.reader(stream, delimiter='\t'))
def field(key):
    return next(row[1:] for row in rows if row and row[0] == key)
titles = field('!Sample_title')
gsm = field('!Sample_geo_accession')
chars = [row[1:] for row in rows if row and row[0] == '!Sample_characteristics_ch1']
sex = next(row for row in chars if all(v.startswith('Sex: ') for v in row))
group = next(row for row in chars if all(v.startswith('primary disease_characterization: ') for v in row))
def count_name(title):
    for prefix in ['RV-Control_', 'RV-Compensated_', 'RV-Decompensated_']:
        if title.startswith(prefix):
            return title[3:]
    return title
annotations = {}
for title, accession, sx, grp in zip(titles, gsm, sex, group):
    sample = count_name(title)
    assert sample not in annotations
    annotations[sample] = dict(GEO_title=title, GSM=accession, Sex=sx.split(': ', 1)[1],
                               Group={'Normal RV':'Ctrl','Compensated RV':'cRV','Decompensated RV':'dRV'}[grp.split(': ',1)[1]])
book = openpyxl.load_workbook(source, read_only=True, data_only=True)
sheet = book['info table'];sheet.reset_dimensions();info = list(sheet.values)
header = next(i for i,row in enumerate(info) if row and row[0] == 'id' and len(row)>1 and row[1]=='cluster ids')
fastq = {row[0]: row for row in info[header+1:] if row and row[0] in annotations and len(row)>2 and isinstance(row[2],str) and '.fastq' in row[2]}
star_header = next(i for i,row in enumerate(info) if row and row[0]=='id' and len(row)>2 and row[1]=='parameters' and row[2]=='reads processed')
feature_header = next(i for i in range(star_header+1,len(info)) if info[i] and info[i][0]=='FEATURE-COUNT:')
star = {row[0]:str(row[1]) for row in info[star_header+1:feature_header] if row and row[0] in annotations}
counts = book['count matrix'];counts.reset_dimensions();samples = list(next(iter(counts.values)))[1:]
assert len(samples)==40 and set(samples)==set(annotations)==set(fastq)==set(star)
metadata = []
for sample in samples:
    row = fastq[sample]
    r1 = row[2];r2 = row[3] if len(row)>3 and row[3] else ''
    single = 'mode: single-end' in star[sample]
    if single:
        assert not r2 and r1.startswith('Prakash_')
        batch = 'SE15'
    else:
        assert r1.startswith('US-') and r2 and r1.replace('_R1','_R2')==r2
        batch = 'PE25'
    metadata.append(dict(Sample=sample,**annotations[sample],Batch=batch,FASTQ_R1=r1,FASTQ_R2=r2,
       Batch_evidence='info table per-sample FASTQ mapping and STAR parameters; GEO mapping-b2 identifies single-end protocol',
       STAR_parameters=star[sample]))
assert collections.Counter(row['Batch'] for row in metadata)=={'SE15':15,'PE25':25}
assert collections.Counter(row['Group'] for row in metadata)=={'Ctrl':13,'cRV':14,'dRV':13}
assert next(row for row in metadata if row['Sample']=='RV-Failing_5')['Group']=='cRV'
with (project/'data/metadata/GSE240921_sample_metadata.csv').open('w',newline='',encoding='utf-8-sig') as stream:
    writer=csv.DictWriter(stream,fieldnames=metadata[0].keys());writer.writeheader();writer.writerows(metadata)
provenance = {
 'count_workbook':source.name,'count_workbook_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'GEO_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240921/matrix/GSE240921_series_matrix.txt.gz',
 'GEO_sha256':hashlib.sha256(geo.read_bytes()).hexdigest(),
 'group_and_sex':'Directly matched GEO primary annotations and sex for all 40; documented RV- prefix normalization for first 15 titles.',
 'batch':'Operational technical batch PE25/SE15, supported by sample-level paired FASTQ versus explicitly single-end STAR records. Not inferred from expression clusters, donor country or clinical A/B labels.',
 'paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11358157/; Methods and Figure3 describe two human sequencing batches with paired-/single-end layouts.',
 'unavailable':'Sample-level age, ischemic time, cause of death and tissue-procurement details could not be reliably linked to all 40 expression samples. GEO provides sex. Sequencing layout was recovered separately from the per-sample FASTQ/STAR records.',
 'scope':'Supports the primary Batch+Sex+Group reanalysis and Figure S8 batch-stratum sensitivity. No SVA/RUV or adjustment for unavailable clinical covariates is claimed.'}
(project/'data/metadata/GSE240921_metadata_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
print('Recovered metadata for 40 samples; group x batch:',dict(collections.Counter((r['Group'],r['Batch']) for r in metadata)))
print('Sex x group:',dict(collections.Counter((r['Sex'],r['Group']) for r in metadata)))
