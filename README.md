# PH RV remodeling revision project

This folder is the reproducible analysis snapshot for the revised manuscript. It contains the analysis inputs, metadata and reference files, executable scripts, final figures, source-data workbooks, and supplementary tables. The previous working projects remain unchanged outside this revision package.

## Reproduce the outputs

The workflow was run with R 4.5.0. R package versions are recorded in `outputs/software_versions/R_package_versions.csv` and the Python/openpyxl versions are in `outputs/software_versions/Python_versions.txt`. Table S2 requires Python with `openpyxl`; set `PHRV_PYTHON` if Python is not on `PATH`.

From this folder in PowerShell:

```powershell
$env:PHRV_PYTHON = "C:/path/to/python.exe"
Rscript scripts/run_all_figures.R
```

The master runner executes Figures 1–6, Figures S1–S8, Tables S1–S4, and the software-version record. It runs each R analysis in a separate process, stops on error, and writes logs to `outputs/run_logs/`. Running it overwrites generated files under `outputs/` and the derived intermediate `data/processed/Figure1_source_data.xlsx`.

## Project structure

- `data/raw/`: downloaded GEO processed matrices and published supplementary workbooks used as analysis inputs.
- `data/reference/`: frozen GEO series-matrix snapshots and rat–human orthology reference.
- `data/metadata/`: sample annotation crosswalks, metadata provenance, preprocessing/model audit, and count-filter audit.
- `data/processed/`: derived intermediate data written by the workflow.
- `scripts/`: figure, table, helper, audit, and master-runner scripts.
- `outputs/`: figures, source-data workbooks, supplementary tables, run logs, and software records.

Dataset accessions, sample sizes, input scales, contrasts, exclusions, identifier handling, filters, and model formulas are documented in Table S1 under `outputs/TableS1_outputs/`. Public inputs are bundled in `data/raw/` and `data/reference/` so the workflow can run without downloading them again. See `data/README.md` for the data-folder guide.

## Key analysis details and interpretation limits

All 40 human RV samples from GSE240921 were retained: 13 controls, 14 compensated, and 13 decompensated. The primary DESeq2 model is `~Batch + Sex + Group`, contrasting decompensated RV with control while retaining compensated samples. `RV-Failing_5` follows its GEO compensated-RV annotation. The 15 single-end and 25 paired-end libraries define operational sequencing-layout strata (SE15 and PE25); they are not assigned donor-country or center labels and do not count as separate consensus votes. Figure S8 reports the joint and within-stratum analyses. Age, ischemic time, cause of death, and tissue-procurement details could not be reliably linked to all samples; no SVA/RUV adjustment is claimed.

The primary Shared-derived panel contains 35 genes, and the non-anchored RV-wide panel contains 97. The 35-gene panel is conditional on discovery eligibility and the prespecified support rule. External RV cohorts contribute to panel derivation, not unbiased validation of the finalized panel. Human support is reported separately in Figures 5 and S5 source data. Consensus stability is threshold-dependent: 14/35 genes remain in all eight conditional LODO runs and 21/35 in two runs; exactly 6/8-supported genes survive two runs by construction. Jaccard indices, omissions, threshold sensitivities, and exploratory effect synthesis are in Figure S6 and Table S3.

GSE198618 provides normalized counts; its rounded-count DESeq2 fit is an explicitly exploratory approximation. GSE133402 provides estimated counts without transcript-length offsets. GSE266139 source Sacrifice and Library factors could not be reliably recovered at sample level; direct LV–RV contrasts are exploratory where animal pairing is unverified. CTEPH effects were imported from published log2 fold-change and adjusted-P results and are assessed descriptively. See Table S1 and the manuscript for details and limitations.

## Deliverables

The project contains Figures 1–6 and S1–S8, their source-data workbooks, and Tables S1–S4. The sibling `submission_documents/` directory contains the tracked manuscript, reviewer response, and reviewer comments for the revision package. The `FILE_MANIFEST.csv` records relative paths, file sizes, and SHA256 hashes for the project snapshot.
