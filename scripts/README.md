# PH ventricular-remodeling figure scripts

## Standard file names

- `Fig1.R` to `Fig6.R`: main figures.
- `FigS1.R` to `FigS4.R`: supplementary figures.
- `00_common_helpers.R`: shared functions required by Fig4, Fig6 and FigS2-S4.

Keep all 11 files together in the project `scripts/` directory. The project root
must contain `data/raw/`, `data/processed/` and `outputs/`.

## Recommended run order

1. `Fig1.R`
2. `Fig2.R`
3. `Fig3.R`
4. `Fig5.R`
5. `Fig4.R`
6. `Fig6.R`
7. `FigS1.R`
8. `FigS2.R`
9. `FigS3.R`
10. `FigS4.R`

Fig4 and FigS4 depend on the final Figure 5 source workbook. Fig6 and FigS2
also use the final 36-gene Figure 5 signature. FigS1 deliberately reruns
`Fig5.R` before applying its independent, non-anchored 6-of-8 rule.

## Outputs

Each script saves a complete figure with a standard name (`Figure1.pdf`,
`FigureS1.pdf`, etc.) and also saves every constituent panel as a separate PDF.
Panel width and height are defined immediately beside the individual `ggsave()`
calls, so later layout adjustments do not require changing the analysis.

Fig5, FigS1 and FigS2 are single-panel figures; their panel aliases are
`Fig5A_bubble.pdf`, `FigS1A_consensus_bubble.pdf` and
`FigS2A_full_CTEPH_heatmap.pdf`.

## Figure 4C color scheme

Figure 4C uses three fixed colors:

- discovery absolute log2FC: orange (`#E68613`)
- LV-only versus Shared: blue (`#1F78B4`)
- RV-only versus Shared: green (`#33A02C`)

Both confidence intervals and filled points use these colors.

## Reproducibility parameters

- DEG threshold: adjusted P < 0.05 and absolute log2FC >= 1.
- Figure 5: Shared discovery gene, significant in GSE266139 RV, with at least
  5 of 7 external RV datasets significant in the anchor direction.
- Figure S1: at least 6 of all 8 RV datasets significant in one concordant
  direction; no forced anchor.
- Figure S3: leave-one-dataset-out proportional rule of at least 6 of 7.
- Figure S4: effect-size-matched resampling; default 10,000 draws. Change with
  environment variable `N_PERMUTATIONS` if needed.

