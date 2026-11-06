# Analysis input data

This folder contains the input files and frozen references used by the analysis scripts.

- `raw/` contains downloaded GEO processed matrices and published supplementary workbooks used as inputs.
- `reference/` contains the GEO series-matrix snapshots and curated rat–human orthology reference used for sample annotation and cross-species comparisons.
- `metadata/` contains sample crosswalks, metadata provenance, the preprocessing/model audit, and the common count-filter audit.
- `processed/` contains derived intermediate data written by the workflow.

Dataset accessions, sample sizes, input measurement scales, contrasts, preprocessing, and exclusions are documented in `outputs/TableS1_outputs/` and the main project README. The raw source files are included so the workflow can run without re-downloading these inputs. Source publications and GEO records remain the original data providers; see Table S1 and the reference list in the manuscript for attribution.
