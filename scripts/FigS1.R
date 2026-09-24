# ================================================================
# Supplementary Figure S1
# Cross-ventricular behavior of threshold-defined LV-only and RV-only genes
#
# Purpose:
#   Demonstrate that "LV-only" and "RV-only" denote significance in only
#   one ventricle, rather than an absent or opposite effect in the other.
#
# Required upstream output:
#   outputs/Fig1_outputs/Figure1_source_data.xlsx
#
# Outputs:
#   outputs/FigS1_outputs/
#     FigureS1.pdf / FigureS1.tiff
#     FigS1A_only_gene_scatter.pdf
#     FigS1B_counterpart_effect.pdf
#     FigS1C_direction_proportion.pdf
#     FigureS1_source_data.xlsx
# ================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tibble)
})

# ---------- 0. Relative project paths ----------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args)
  if (length(hit)) {
    return(dirname(normalizePath(sub("^--file=", "", args[hit[1]]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  normalizePath(getwd())
}

find_project_dir <- function(script_dir = get_script_dir()) {
  candidates <- unique(normalizePath(
    c(file.path(script_dir, ".."), script_dir, file.path(script_dir, "../..")),
    mustWork = FALSE
  ))
  ok <- vapply(
    candidates,
    function(x) dir.exists(file.path(x, "data", "raw")),
    logical(1)
  )
  if (!any(ok)) {
    stop("Cannot locate project root containing data/raw. Put this script in project/scripts.")
  }
  candidates[which(ok)[1]]
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
fig1_source <- file.path(
  project_dir, "outputs", "Fig1_outputs", "Figure1_source_data.xlsx"
)
out_dir <- file.path(project_dir, "outputs", "FigS1_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(fig1_source)) {
  stop("Cannot find Figure 1 source data. Run Fig1.R first: ", fig1_source)
}

# ---------- 1. Parameters and style ----------
alpha <- 0.05
lfc_thr <- 1
set_cols <- c("LV-only" = "#1F78B4", "RV-only" = "#33A02C")
direction_cols <- c("Concordant" = "#D73027", "Discordant" = "#4575B4")

theme_s1 <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey25"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(7, 8, 7, 7)
    )
}

read_set <- function(sheet, set_name) {
  x <- openxlsx::read.xlsx(fig1_source, sheet = sheet) %>% as_tibble()
  gene_col <- intersect(c("Genename", "Gene", "SYMBOL", "GeneSymbol"), names(x))[1]
  required <- c("lfc_lv", "padj_lv", "lfc_rv", "padj_rv")
  if (is.na(gene_col) || any(!required %in% names(x))) {
    stop("Unexpected columns in Figure 1 sheet: ", sheet)
  }
  x %>%
    transmute(
      Gene = toupper(trimws(as.character(.data[[gene_col]]))),
      Set = set_name,
      lfc_lv = suppressWarnings(as.numeric(lfc_lv)),
      padj_lv = suppressWarnings(as.numeric(padj_lv)),
      lfc_rv = suppressWarnings(as.numeric(lfc_rv)),
      padj_rv = suppressWarnings(as.numeric(padj_rv))
    ) %>%
    filter(!is.na(Gene), Gene != "", is.finite(lfc_lv), is.finite(lfc_rv)) %>%
    distinct(Gene, .keep_all = TRUE)
}

only_genes <- bind_rows(
  read_set("LV_only_genes", "LV-only"),
  read_set("RV_only_genes", "RV-only")
) %>%
  mutate(
    Set = factor(Set, levels = c("LV-only", "RV-only")),
    significant_log2FC = if_else(Set == "LV-only", lfc_lv, lfc_rv),
    counterpart_log2FC = if_else(Set == "LV-only", lfc_rv, lfc_lv),
    counterpart_padj = if_else(Set == "LV-only", padj_rv, padj_lv),
    signed_counterpart_log2FC = counterpart_log2FC * sign(significant_log2FC),
    counterpart_direction = if_else(
      signed_counterpart_log2FC >= 0,
      "Concordant", "Discordant"
    ),
    counterpart_meets_full_threshold =
      !is.na(counterpart_padj) & counterpart_padj < alpha &
      abs(counterpart_log2FC) >= lfc_thr
  )

# Internal audit: by definition, no "only" gene may cross the full threshold
# in the other ventricle.
if (any(only_genes$counterpart_meets_full_threshold, na.rm = TRUE)) {
  stop("Set-definition audit failed: an only gene is significant in both ventricles.")
}

direction_summary <- only_genes %>%
  count(Set, counterpart_direction, name = "n") %>%
  group_by(Set) %>%
  mutate(
    total = sum(n),
    proportion = n / total,
    percent = 100 * proportion
  ) %>%
  ungroup()

annotation_df <- direction_summary %>%
  filter(counterpart_direction == "Concordant") %>%
  transmute(
    Set,
    label = paste0(round(percent, 1), "% concordant")
  )

# Label a small number of the most reproducibly large concordant effects.
label_df <- only_genes %>%
  filter(counterpart_direction == "Concordant") %>%
  mutate(label_score = abs(lfc_lv) + abs(lfc_rv)) %>%
  group_by(Set) %>%
  slice_max(label_score, n = 5, with_ties = FALSE) %>%
  ungroup()

# ---------- 2. Panel A: effect-size scatter restricted to only genes ----------
pA <- ggplot(only_genes, aes(lfc_lv, lfc_rv, color = Set)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_hline(
    yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
    color = "grey50", linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
    color = "grey50", linewidth = 0.4
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.55) +
  geom_point(size = 1.9, alpha = 0.72) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 2.8,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = set_cols, drop = FALSE) +
  labs(
    title = "A. Cross-ventricular effect sizes",
    subtitle = "Dashed lines indicate |log2FC| = 1; dotted line indicates equal LV and RV effects",
    x = "LV log2FC (MCT vs control)",
    y = "RV log2FC (MCT vs control)",
    color = NULL
  ) +
  theme_s1(10) +
  theme(legend.position = "top")

# ---------- 3. Panel B: counterpart effect aligned to discovery direction ----------
pB <- ggplot(
  only_genes,
  aes(Set, signed_counterpart_log2FC, fill = Set)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.42, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.75, alpha = 0.30) +
  geom_text(
    data = annotation_df,
    aes(x = Set, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.4,
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(values = set_cols, guide = "none") +
  labs(
    title = "B. Effect in the non-significant ventricle",
    subtitle = "Positive values indicate the same direction as the significant ventricle",
    x = NULL,
    y = "Direction-aligned counterpart log2FC"
  ) +
  theme_s1(10)

# ---------- 4. Panel C: direct direction summary ----------
pC <- ggplot(
  direction_summary,
  aes(Set, proportion, fill = counterpart_direction)
) +
  geom_col(width = 0.66, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = if_else(percent >= 5, paste0(round(percent, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(values = direction_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "C. Direction in the non-significant ventricle",
    x = NULL,
    y = "Proportion of genes",
    fill = NULL
  ) +
  theme_s1(10) +
  theme(legend.position = "top")

# ---------- 5. Combined figure and individual panels ----------
fig_s1 <- (pA | (pB / pC)) +
  patchwork::plot_layout(widths = c(1.35, 1))

ggsave(
  file.path(out_dir, "FigureS1.pdf"),
  fig_s1, width = 11.0, height = 6.9, units = "in"
)
ggsave(
  file.path(out_dir, "FigureS1.tiff"),
  fig_s1, width = 11.0, height = 6.9, units = "in",
  dpi = 600, compression = "lzw"
)
ggsave(
  file.path(out_dir, "FigS1A_only_gene_scatter.pdf"),
  pA, width = 6.0, height = 5.6, units = "in"
)
ggsave(
  file.path(out_dir, "FigS1B_counterpart_effect.pdf"),
  pB, width = 4.8, height = 3.4, units = "in"
)
ggsave(
  file.path(out_dir, "FigS1C_direction_proportion.pdf"),
  pC, width = 4.8, height = 3.4, units = "in"
)

openxlsx::write.xlsx(
  list(
    Only_gene_effects = only_genes,
    Direction_summary = direction_summary,
    Meta = data.frame(
      item = c("LV-only definition", "RV-only definition", "Interpretation"),
      value = c(
        "padj < 0.05 and |log2FC| >= 1 in LV, but not the full threshold in RV",
        "padj < 0.05 and |log2FC| >= 1 in RV, but not the full threshold in LV",
        "Only denotes threshold-defined significance, not biological absence in the other ventricle"
      )
    )
  ),
  file.path(out_dir, "FigureS1_source_data.xlsx"),
  overwrite = TRUE
)

message("Supplementary Figure S1 finished: ", out_dir)
