# R3.3 analysis and package verification

Run date: 2026-11-04. The complete `scripts/run_all_figures.R` workflow finished successfully with R 4.5.0 and `PHRV_PYTHON` set to the Python runtime with `openpyxl`.

## Executed outputs

- Figures 1–6 and Figures S1–S8 were regenerated in the packaged project.
- Tables S1–S4 and software-version records were regenerated.
- Per-script logs are in `outputs/run_logs/`; the master-run log is `outputs/run_logs/R3_3_full_run.log`.
- Figure S8 shows the technical-stratum sensitivity, separate group-only and sex-adjusted models, the fixed 35-gene panel, and the all-sample batch-by-disease interaction. Source data include sample/GSM metadata, fitted effects, standard errors, adjusted P values, and design-rank checks.

## Verified GSE240921 results

All 40 samples were included in the primary model (`~Batch + Sex + Group`): 13 controls, 14 compensated RV, and 13 decompensated RV. GEO supplied the group and sex annotations. Sample-linked FASTQ/STAR records supported 15 single-end (SE15) and 25 paired-end (PE25) libraries. GSE240921 remains one cohort-level contribution to consensus selection.

Sex-adjusted within-stratum models yielded 4,734 strict DEGs in SE15 and three in PE25; the joint model yielded 855. For the fixed 35-gene panel, strict discovery-direction support was present for 22/34 evaluable genes in SE15, 0/26 in PE25, and 7/33 in the joint model. Among 32 panel genes with finite effects in both strata, 19 matched the discovery direction in both. Across 19,784 common finite transcriptome-wide effects, the sex-adjusted stratum estimates had Spearman rho=0.144.

The all-sample interaction model tested `(dRV-Ctrl)SE15 − (dRV-Ctrl)PE25` with sex adjustment. It identified 1,021 genes at transcriptome-wide FDR<0.05, including PLOD2, POSTN, PRG4, SCN3B, SRPX2, and STC2 from the fixed panel. These results indicate heterogeneous effects; they do not establish whether sequencing layout or clinical/donor background explains the heterogeneity.

The primary panel remained 35 genes; TUBB3 moved from 7/8 to 6/8 cohort support. The non-anchored RV-wide panel changed from 104 to 97 genes. Human-support, effect-matching, weighting, LODO, held-out, protein-annotation, and source-data outputs that depend on GSE240921 were refreshed by the master run.

## Interpretation limits

Age, ischemic time, cause of death, and tissue-procurement details were not reliably linked to all 40 expression samples. No SVA/RUV adjustment was performed. PCA and DEG-count differences were not used to exclude samples or label a stratum as poor quality. The two sequencing-layout strata were not counted as independent datasets. These analyses do not eliminate residual confounding.

The updated manuscript and response documents are outside this analysis directory. LibreOffice was unavailable in the environment, so the DOCX files could not be rendered for page-layout inspection.
