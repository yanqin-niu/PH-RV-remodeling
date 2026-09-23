# Figure 6 — three-panel CTEPH validation of the current Figure 5 signature
#
# Required upstream output:
#   outputs/Fig5_outputs/Figure5_source_data.xlsx
#
# Panels:
#   A. Signature activity across pre-PEA disease strata
#   B. Regional and post-PEA-associated signature summaries
#   C. Heatmap of the most recurrent signature genes across all 12 contrasts

suppressPackageStartupMessages({
  library(patchwork)
})

.args <- commandArgs(trailingOnly = FALSE)
.hit  <- grep("^--file=", .args)
.dir  <- if (length(.hit)) {
  dirname(normalizePath(sub("^--file=", "", .args[.hit[1]])))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
    nzchar(rstudioapi::getActiveDocumentContext()$path)
) {
  dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))
} else {
  normalizePath(getwd())
}

source(file.path(.dir, "00_common_helpers.R"))

# ============================================================
# 1. Paths and input data
# ============================================================
script_dir  <- get_script_dir()
project_dir <- find_project_dir(script_dir)
out_dir     <- file.path(project_dir, "outputs", "Fig6_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f5 <- read_fig5_objects(project_dir)
signature_genes <- unique(TOUP(f5$consensus$Gene))
signature_genes <- signature_genes[!is.na(signature_genes) & signature_genes != ""]

ct  <- read_cteph_data(project_dir, signature_genes)
dat <- ct$data
sm  <- signature_summary(dat)

contrast_order <- as.character(ct$spec$label)

# Consistent colors used in the previous Figure 6 version.
col_disease  <- "#D73027"
fill_disease <- "#FDD0A2"
col_regional <- "#4575B4"
col_postpea  <- "#8073AC"

# ============================================================
# 2. Panel A — disease-stratum signature activity
# ============================================================
disease_order <- c(
  "Moderate pre-PEA RV vs control RV",
  "Intermediate pre-PEA RV vs control RV",
  "Severe pre-PEA RV vs control RV"
)

disease_labels <- c(
  "Moderate pre-PEA RV vs control RV"     = "Moderate",
  "Intermediate pre-PEA RV vs control RV" = "Intermediate",
  "Severe pre-PEA RV vs control RV"       = "Severe"
)

disease <- dat %>%
  dplyr::filter(Contrast %in% disease_order) %>%
  dplyr::mutate(
    Stratum = factor(
      disease_labels[as.character(Contrast)],
      levels = c("Moderate", "Intermediate", "Severe")
    )
  )

p6A <- ggplot(disease, aes(x = Stratum, y = log2FC)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_violin(
    fill = fill_disease,
    color = col_disease,
    linewidth = 0.55,
    scale = "width",
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.17,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.42
  ) +
  geom_jitter(
    width = 0.075,
    height = 0,
    size = 0.85,
    alpha = 0.45,
    color = "#B2182B"
  ) +
  labs(
    title = "A. Signature activity in pre-PEA right ventricle",
    subtitle = paste0("Figure 5 signature; n = ", length(signature_genes), " genes"),
    x = NULL,
    y = "Gene-level log2 fold change"
  ) +
  theme_manuscript(10.5) +
  theme(
    plot.title = element_text(size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 9, hjust = 0, color = "grey30"),
    axis.text.x = element_text(size = 9.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    plot.margin = margin(6, 8, 5, 5)
  )

# ============================================================
# 3. Panel B — regional and post-PEA-associated summaries
# ============================================================
other_label_map <- c(
  "Pre-PEA septum vs control septum" =
    "Pre-PEA septum\nvs control septum",
  "Pre-PEA septum vs pre-PEA RV" =
    "Pre-PEA septum\nvs pre-PEA RV",
  "Post-PEA septum vs pre-PEA RV" =
    "Post-PEA septum\nvs pre-PEA RV",
  "Post-PEA moderate septum vs pre-PEA moderate RV" =
    "Post-PEA moderate septum\nvs pre-PEA moderate RV",
  "Post-PEA intermediate septum vs pre-PEA intermediate RV" =
    "Post-PEA intermediate septum\nvs pre-PEA intermediate RV",
  "Post-PEA severe septum vs pre-PEA severe RV" =
    "Post-PEA severe septum\nvs pre-PEA severe RV"
)

other_sm <- sm %>%
  dplyr::filter(
    Contrast_group %in% c("Regional comparison", "Post-PEA-associated")
  ) %>%
  dplyr::mutate(
    Contrast_short = unname(other_label_map[as.character(Contrast)]),
    Contrast_short = factor(Contrast_short, levels = rev(unname(other_label_map)))
  )

p6B <- ggplot(
  other_sm,
  aes(x = mean_log2FC, y = Contrast_short, color = Contrast_group)
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_segment(
    aes(x = q25, xend = q75, yend = Contrast_short),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      "Regional comparison" = col_regional,
      "Post-PEA-associated" = col_postpea
    ),
    breaks = c("Regional comparison", "Post-PEA-associated")
  ) +
  labs(
    title = "B. Regional and post-PEA-associated contrasts",
    subtitle = "Point: mean; line: interquartile range across signature genes",
    x = "Signature log2 fold change",
    y = NULL,
    color = NULL
  ) +
  theme_manuscript(9.5) +
  theme(
    plot.title = element_text(size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 8.8, hjust = 0, color = "grey30"),
    axis.text.y = element_text(size = 8.2, lineheight = 0.9),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 1, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    plot.margin = margin(6, 5, 5, 8)
  )

# ============================================================
# 4. Panel C — recurrent-gene heatmap across all contrasts
# ============================================================
gene_rank <- dat %>%
  dplyr::mutate(
    significant = !is.na(padj) & padj < 0.05 &
      !is.na(log2FC) & abs(log2FC) >= 1
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    n_significant = sum(significant, na.rm = TRUE),
    n_measured = sum(!is.na(log2FC)),
    mean_abs_log2FC = mean(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_abs_log2FC = dplyr::if_else(
      is.nan(mean_abs_log2FC), NA_real_, mean_abs_log2FC
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_significant),
    dplyr::desc(mean_abs_log2FC),
    Gene
  ) %>%
  dplyr::slice_head(n = 20)

contrast_short_map <- c(
  "Moderate vs severe pre-PEA RV" = "Moderate vs\nsevere RV",
  "Moderate vs intermediate pre-PEA RV" = "Moderate vs\nintermediate RV",
  "Intermediate vs severe pre-PEA RV" = "Intermediate vs\nsevere RV",
  "Moderate pre-PEA RV vs control RV" = "Moderate RV\nvs control",
  "Intermediate pre-PEA RV vs control RV" = "Intermediate RV\nvs control",
  "Severe pre-PEA RV vs control RV" = "Severe RV\nvs control",
  "Pre-PEA septum vs control septum" = "Pre-PEA septum\nvs control",
  "Pre-PEA septum vs pre-PEA RV" = "Pre-PEA septum\nvs pre-PEA RV",
  "Post-PEA septum vs pre-PEA RV" = "Post-PEA septum\nvs pre-PEA RV",
  "Post-PEA moderate septum vs pre-PEA moderate RV" = "Post moderate\nvs pre moderate",
  "Post-PEA intermediate septum vs pre-PEA intermediate RV" = "Post intermediate\nvs pre intermediate",
  "Post-PEA severe septum vs pre-PEA severe RV" = "Post severe\nvs pre severe"
)

heat <- ct$full %>%
  dplyr::filter(Gene %in% gene_rank$Gene) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(gene_rank$Gene)),
    Contrast_short = unname(contrast_short_map[as.character(Contrast)]),
    Contrast_short = factor(
      Contrast_short,
      levels = unname(contrast_short_map[contrast_order])
    ),
    missing_value = is.na(log2FC),
    significance_mark = dplyr::case_when(
      missing_value ~ "NA",
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01  ~ "**",
      !is.na(padj) & padj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

finite_abs_lfc <- abs(heat$log2FC[is.finite(heat$log2FC)])
fill_limit <- if (length(finite_abs_lfc)) {
  max(1, as.numeric(stats::quantile(finite_abs_lfc, 0.98, na.rm = TRUE)))
} else {
  1
}

heat <- heat %>%
  dplyr::mutate(
    log2FC_display = pmax(-fill_limit, pmin(fill_limit, log2FC))
  )

p6C <- ggplot(heat, aes(x = Contrast_short, y = Gene, fill = log2FC_display)) +
  geom_tile(color = "white", linewidth = 0.30) +
  geom_text(
    aes(label = significance_mark, color = missing_value),
    size = 2.55,
    lineheight = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "grey35")) +
  scale_fill_gradient2(
    low = "blue",
    mid = "#FFFF7F",
    high = "red",
    midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    na.value = "grey82"
  ) +
  labs(
    title = paste0(
      "C. Most recurrent CTEPH-supported genes (",
      nrow(gene_rank), "/", length(signature_genes), ")"
    ),
    subtitle = paste0(
      "All 12 contrasts; grey/NA = not reported or not estimable; ",
      "* padj < 0.05, ** padj < 0.01, *** padj < 0.001"
    ),
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  theme_manuscript(9.5) +
  theme(
    plot.title = element_text(size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 8.8, hjust = 0, color = "grey30"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 7.5,
      lineheight = 0.88
    ),
    axis.text.y = element_text(size = 8.3, face = "italic"),
    panel.border = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(15, "mm"),
    plot.margin = margin(7, 8, 5, 5)
  )

# ============================================================
# 5. Final three-panel composition
# ============================================================
top_row <- (p6A | p6B) +
  patchwork::plot_layout(widths = c(0.82, 1.18))

fig6 <- top_row / p6C +
  patchwork::plot_layout(heights = c(0.86, 1.58)) +
  patchwork::plot_annotation(
    theme = theme(plot.margin = margin(5, 7, 5, 5))
  )

# ============================================================
# 6. Export the combined figure
# ============================================================
ggsave(
  file.path(out_dir, "Figure6.pdf"),
  fig6,
  width = 13.5,
  height = 10.0,
  units = "in"
)

ggsave(
  file.path(out_dir, "Figure6.tiff"),
  fig6,
  width = 13.5,
  height = 10.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# 7. Export every panel separately
#    Change these width/height values if you later adjust layout.
# ============================================================
ggsave(
  file.path(out_dir, "Fig6A_disease_strata.pdf"),
  p6A,
  width = 5.3,
  height = 4.1,
  units = "in"
)

ggsave(
  file.path(out_dir, "Fig6B_regional_postPEA.pdf"),
  p6B,
  width = 7.0,
  height = 4.1,
  units = "in"
)

ggsave(
  file.path(out_dir, "Fig6C_recurrent_gene_heatmap.pdf"),
  p6C,
  width = 13.2,
  height = 6.3,
  units = "in"
)

# Optional high-resolution TIFF files for individual panels.
ggsave(
  file.path(out_dir, "Fig6A_disease_strata.tiff"),
  p6A,
  width = 5.3,
  height = 4.1,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(out_dir, "Fig6B_regional_postPEA.tiff"),
  p6B,
  width = 7.0,
  height = 4.1,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(out_dir, "Fig6C_recurrent_gene_heatmap.tiff"),
  p6C,
  width = 13.2,
  height = 6.3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# 8. Source-data workbook
# ============================================================
openxlsx::write.xlsx(
  list(
    Signature_genes = data.frame(Gene = signature_genes),
    All_CTEPH_data = dat,
    Contrast_summary = sm,
    PanelA_gene_values = disease,
    PanelB_summary = other_sm,
    PanelC_gene_ranking = gene_rank,
    PanelC_full_grid = heat,
    Meta = data.frame(
      Figure5_source = f5$path,
      Signature_size = length(signature_genes),
      Heatmap_gene_number = nrow(gene_rank),
      Note = paste(
        "Post-PEA-associated comparisons are descriptive because",
        "anatomical region and sampling time both differ."
      )
    )
  ),
  file.path(out_dir, "Figure6_source_data.xlsx"),
  overwrite = TRUE
)

print(fig6)
message("Figure 6 finished: ", out_dir)
