"""Validate current deliverables and create a local review package; never uploads."""
from pathlib import Path
import csv
import hashlib
import shutil
import zipfile
import json
import openpyxl

root = Path(__file__).resolve().parents[1]
destination = root.parent / 'review_materials_final'
destination.mkdir(exist_ok=True)
figures = [f'Fig{i}' for i in range(1, 7)] + [f'FigS{i}' for i in range(1, 9)]
sources = []
supplements = []
for name in figures:
    stem = name.replace('Fig', 'Figure', 1)
    folder = root / 'outputs' / f'{name}_outputs'
    for filename in [f'{stem}.pdf', f'{stem}_source_data.xlsx']:
        assert (folder / filename).is_file(), folder / filename
    sources.append(folder / f'{stem}_source_data.xlsx')
    if name.startswith('FigS'):
        supplements.append(folder / f'{stem}.pdf')

tables = [root / 'outputs' / f'TableS{i}_outputs' / name for i, name in enumerate([
    'TableS1_R2_5.xlsx', 'TableS2_independence_weighting.xlsx',
    'TableS3_consensus_stability.xlsx', 'TableS4_heldout_assessment.xlsx'], 1)]
# Locate the canonical Table S4 workbook if its descriptive suffix differs.
if not tables[-1].exists():
    candidates = list((root / 'outputs/TableS4_outputs').glob('TableS4*.xlsx'))
    assert len(candidates) == 1, candidates
    tables[-1] = candidates[0]
for path in tables + sources:
    assert path.is_file(), path
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    assert workbook.sheetnames, path
    workbook.close()

with zipfile.ZipFile(destination / 'Supplementary_Figures_S1-S8.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for path in supplements:
        archive.write(path, path.name)
with zipfile.ZipFile(destination / 'Supplementary_Tables_S1-S4.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for i, path in enumerate(tables, 1):
        archive.write(path, f'Table_S{i}.xlsx')
with zipfile.ZipFile(destination / 'Source_data_and_scripts.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for path in sources:
        archive.write(path, path.relative_to(root))
    for folder in ['scripts', 'data/metadata', 'data/reference', 'outputs/software_versions']:
        for path in sorted((root / folder).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.name not in {'.RDataTmp', '.Rhistory'}:
                archive.write(path, path.relative_to(root))
    for filename in ['README.md', 'VALIDATION_R3_3.md', 'output_provenance.json', 'FILE_MANIFEST.csv']:
        archive.write(root / filename, filename)
with zipfile.ZipFile(destination / 'Analysis_input_data.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for path in sorted((root / 'data/raw').glob('*')):
        if path.is_file():
            archive.write(path, path.relative_to(root))

report = {'figures': 14, 'supplementary_figures': 8, 'source_workbooks': len(sources),
          'supplementary_tables': len(tables), 'source_project': root.name,
          'raw_inputs_in_separate_archive': True,
          'publication_status': 'LOCAL ONLY; files must still be attached in the journal submission system and the revised repository must be uploaded.'}
(destination / 'PACKAGE_VALIDATION.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
(destination / 'READ_ME_FIRST.txt').write_text('''Review materials for PH-RV remodeling revision

Supplementary_Figures_S1-S8.zip: eight supplementary figure PDFs.
Supplementary_Tables_S1-S4.zip: four supplementary table workbooks.
Source_data_and_scripts.zip: all 14 figure source workbooks, executable scripts,
metadata/orthology references, README, and actual software-version records.
Analysis_input_data.zip: deposited input files used by the scripts.
Extract Source_data_and_scripts.zip and Analysis_input_data.zip into the SAME
project directory, preserving data/, scripts/ and outputs/. Install the R
dependencies listed in outputs/software_versions/R_package_versions.csv and
provide Python with openpyxl for TableS2.R. Run:
Rscript scripts/run_all_figures.R
Set PHRV_PYTHON to the Python executable if it is not available as python on PATH.
The run overwrites outputs within that extracted project directory.

These are local preparation files. They have not been uploaded to the journal
or to GitHub. Attach the supplement and table files to the revision; do not rely
only on the repository URL. Future analysis edits require rebuilding this bundle.
Software versions describe the execution environment, not a locked environment
or proof of reproduction on a clean machine. Shared data-scale and independence
limitations remain described in the manuscript and project README.
''', encoding='utf-8')
with (destination / 'CHECKSUMS.csv').open('w', newline='', encoding='utf-8-sig') as stream:
    writer = csv.writer(stream)
    writer.writerow(['File', 'Bytes', 'SHA256'])
    for path in sorted(destination.glob('*.zip')):
        with zipfile.ZipFile(path) as archive:
            assert archive.testzip() is None, path
        writer.writerow([path.name, path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest()])
print(json.dumps(report, indent=2))
