# Shared helpers for Figure 4, Figure 6 and Supplementary Figures S2-S4
# Run Figure 1 and the final Figure 5 script before these downstream scripts.

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", args[hit[1]]))))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    z <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(z) && nzchar(z$path)) return(dirname(normalizePath(z$path)))
  }
  normalizePath(getwd())
}

find_project_dir <- function(script_dir = get_script_dir()) {
  candidates <- unique(normalizePath(
    c(file.path(script_dir, ".."), script_dir, file.path(script_dir, "../..")),
    mustWork = FALSE
  ))
  ok <- vapply(candidates, function(x) dir.exists(file.path(x, "data", "raw")), logical(1))
  if (!any(ok)) stop("Cannot locate project root containing data/raw. Put scripts in project/scripts or set working directory inside the project.")
  candidates[which(ok)[1]]
}

TOUP <- function(x) toupper(trimws(as.character(x)))

theme_manuscript <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.45),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank()
    )
}

set_cols <- c("Shared" = "#D73027", "LV-only" = "#1F78B4", "RV-only" = "#33A02C")
direction_cols <- c("Up" = "red", "Down" = "blue")
species_cols <- c("Human" = "#984EA3", "Rat" = "#4DAF4A")

resolve_fig5_source <- function(project_dir) {
  from_env <- Sys.getenv("FIG5_SOURCE_FILE", "")
  candidates <- c(
    from_env,
    file.path(project_dir, "outputs", "Fig5_outputs", "Figure5_source_data.xlsx")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop(
    "Cannot find final Figure5_source_data.xlsx. Run the final Figure 5 script first, ",
    "or set FIG5_SOURCE_FILE to its full path."
  )
  normalizePath(hit[1])
}

read_fig1_sets <- function(project_dir) {
  path <- file.path(project_dir, "outputs", "Fig1_outputs", "Figure1_source_data.xlsx")
  if (!file.exists(path)) stop("Cannot find Figure 1 source data: ", path)
  read_one <- function(sheet, set_name) {
    x <- openxlsx::read.xlsx(path, sheet = sheet)
    g <- intersect(c("Genename", "Gene", "SYMBOL", "GeneSymbol"), names(x))[1]
    if (is.na(g)) stop("No gene-symbol column in Figure 1 sheet: ", sheet)
    lfc <- switch(set_name,
      "Shared" = rowMeans(cbind(as.numeric(x$lfc_lv), as.numeric(x$lfc_rv)), na.rm = TRUE),
      "LV-only" = as.numeric(x$lfc_lv),
      "RV-only" = as.numeric(x$lfc_rv)
    )
    tibble::tibble(
      Gene = TOUP(x[[g]]), Set = set_name,
      discovery_log2FC = lfc,
      discovery_abs_log2FC = abs(lfc),
      discovery_direction = dplyr::if_else(lfc >= 0, "Up", "Down")
    ) %>%
      dplyr::filter(!is.na(Gene), Gene != "", is.finite(discovery_log2FC)) %>%
      dplyr::distinct(Gene, .keep_all = TRUE)
  }
  dplyr::bind_rows(
    read_one("Shared_genes", "Shared"),
    read_one("LV_only_genes", "LV-only"),
    read_one("RV_only_genes", "RV-only")
  ) %>% dplyr::mutate(Set = factor(Set, levels = c("Shared", "LV-only", "RV-only")))
}

read_fig5_objects <- function(project_dir) {
  path <- resolve_fig5_source(project_dir)
  sheets <- openxlsx::getSheetNames(path)
  required <- c("Display_df_8datasets", "Consensus_genes")
  miss <- setdiff(required, sheets)
  if (length(miss)) stop("Figure 5 source workbook lacks sheet(s): ", paste(miss, collapse = ", "))
  display <- openxlsx::read.xlsx(path, sheet = "Display_df_8datasets") %>%
    dplyr::transmute(
      Gene = TOUP(Gene), Dataset = as.character(Dataset),
      log2FC = as.numeric(log2FC), padj = as.numeric(padj),
      Significant = as.logical(Significant), Direction = as.character(Direction)
    ) %>%
    dplyr::mutate(
      Significant = !is.na(padj) & padj < 0.05 & !is.na(log2FC) & abs(log2FC) >= 1,
      Direction = dplyr::case_when(Significant & log2FC > 0 ~ "Up", Significant & log2FC < 0 ~ "Down", TRUE ~ "NS"),
      Species = dplyr::if_else(stringr::str_detect(Dataset, stringr::regex("Human", ignore_case = TRUE)), "Human", "Rat")
    ) %>% dplyr::distinct(Gene, Dataset, .keep_all = TRUE)
  consensus <- openxlsx::read.xlsx(path, sheet = "Consensus_genes")
  consensus$Gene <- TOUP(consensus$Gene)
  list(path = path, display = display, consensus = consensus)
}

build_replication_grid <- function(fig1_sets, display, exclude_anchor = TRUE) {
  d <- display
  if (exclude_anchor) d <- d %>% dplyr::filter(!stringr::str_detect(Dataset, stringr::fixed("GSE266139")))
  datasets <- unique(d$Dataset)
  tidyr::expand_grid(fig1_sets, Dataset = datasets) %>%
    dplyr::left_join(d %>% dplyr::select(Gene, Dataset, log2FC, padj, Significant, Direction, Species), by = c("Gene", "Dataset")) %>%
    dplyr::mutate(
      measured = !is.na(log2FC) | !is.na(padj),
      same_direction = Significant & Direction == discovery_direction,
      opposite_direction = Significant & Direction != "NS" & Direction != discovery_direction,
      Species = dplyr::if_else(stringr::str_detect(Dataset, stringr::regex("Human", ignore_case = TRUE)), "Human", "Rat")
    )
}

wilson_ci <- function(x, n, conf = 0.95) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - half), min(1, ctr + half))
}

summarise_replication <- function(grid) {
  grid %>% dplyr::group_by(Dataset, Species, Set) %>%
    dplyr::summarise(n_total = dplyr::n(), n_measured = sum(measured), n_same = sum(same_direction, na.rm = TRUE),
              n_opposite = sum(opposite_direction, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rowwise() %>% dplyr::mutate(rate = n_same / n_measured, ci = list(wilson_ci(n_same, n_measured)),
                        lower = ci[[1]], upper = ci[[2]]) %>% dplyr::ungroup() %>% dplyr::select(-ci)
}

sheet_col <- function(df, patterns, required = TRUE) {
  z <- names(df)[Reduce(`|`, lapply(patterns, function(p) grepl(p, names(df), ignore.case = TRUE)))]
  if (!length(z) && required) stop("Cannot identify column matching: ", paste(patterns, collapse = " / "))
  if (length(z)) z[1] else NA_character_
}

extract_cteph_sheet <- function(df, key, label, group, signature_genes, reverse = FALSE) {
  names(df) <- trimws(names(df))
  gcol <- if ("Ensembl.gene" %in% names(df)) "Ensembl.gene" else sheet_col(df, c("gene", "symbol"))
  lcol <- sheet_col(df, c("log2", "logFC")); pcol <- sheet_col(df, c("padj", "FDR", "adj"))
  bcol <- sheet_col(df, c("baseMean"), required = FALSE)
  out <- tibble::tibble(
    Gene = TOUP(df[[gcol]]),
    baseMean = if (is.na(bcol)) NA_real_ else suppressWarnings(as.numeric(df[[bcol]])),
    log2FC = suppressWarnings(as.numeric(df[[lcol]])),
    padj = suppressWarnings(as.numeric(df[[pcol]])),
    Contrast_key = key, Contrast = label, Contrast_group = group
  ) %>% dplyr::filter(Gene %in% signature_genes, Gene != "")
  if (reverse) out$log2FC <- -out$log2FC
  out %>%
    dplyr::arrange(Gene, padj, dplyr::desc(abs(log2FC))) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
}

read_cteph_data <- function(project_dir, signature_genes) {
  raw_dir <- file.path(project_dir, "data", "raw")
  spec <- tibble::tribble(
    ~file, ~sheet, ~label, ~group, ~reverse,
    "44161_2025_672_MOESM4_ESM.xlsx", "Fig 2c", "Moderate vs severe pre-PEA RV", "Within pre-PEA RV", FALSE,
    "44161_2025_672_MOESM4_ESM.xlsx", "Fig 2d", "Moderate vs intermediate pre-PEA RV", "Within pre-PEA RV", FALSE,
    "44161_2025_672_MOESM4_ESM.xlsx", "Fig 2e", "Intermediate vs severe pre-PEA RV", "Within pre-PEA RV", FALSE,
    "44161_2025_672_MOESM5_ESM.xlsx", "Fig 3c", "Moderate pre-PEA RV vs control RV", "Disease vs control", FALSE,
    "44161_2025_672_MOESM5_ESM.xlsx", "Fig 3d", "Intermediate pre-PEA RV vs control RV", "Disease vs control", FALSE,
    "44161_2025_672_MOESM5_ESM.xlsx", "Fig 3e", "Severe pre-PEA RV vs control RV", "Disease vs control", FALSE,
    "44161_2025_672_MOESM6_ESM.xlsx", "Fig 4c", "Pre-PEA septum vs control septum", "Regional comparison", FALSE,
    "44161_2025_672_MOESM6_ESM.xlsx", "Fig 4e", "Pre-PEA septum vs pre-PEA RV", "Regional comparison", TRUE,
    "44161_2025_672_MOESM7_ESM.xlsx", "Fig 5c", "Post-PEA septum vs pre-PEA RV", "Post-PEA-associated", FALSE,
    "44161_2025_672_MOESM7_ESM.xlsx", "Fig 5d", "Post-PEA moderate septum vs pre-PEA moderate RV", "Post-PEA-associated", FALSE,
    "44161_2025_672_MOESM7_ESM.xlsx", "Fig 5e", "Post-PEA intermediate septum vs pre-PEA intermediate RV", "Post-PEA-associated", FALSE,
    "44161_2025_672_MOESM7_ESM.xlsx", "Fig 5f", "Post-PEA severe septum vs pre-PEA severe RV", "Post-PEA-associated", FALSE
  )
  paths <- file.path(raw_dir, unique(spec$file))
  if (any(!file.exists(paths))) stop("Missing CTEPH source file(s): ", paste(basename(paths[!file.exists(paths)]), collapse = ", "))
  out <- lapply(seq_len(nrow(spec)), function(i) {
    row <- spec[i, ]; path <- file.path(raw_dir, row$file)
    extract_cteph_sheet(openxlsx::read.xlsx(path, sheet = row$sheet), row$sheet, row$label,
                        row$group, signature_genes, row$reverse)
  }) %>% dplyr::bind_rows()
  full <- tidyr::expand_grid(Gene = signature_genes, Contrast = spec$label) %>%
    dplyr::left_join(out, by = c("Gene", "Contrast")) %>%
    dplyr::left_join(spec %>% dplyr::select(Contrast = label, Contrast_group = group, Contrast_key = sheet) %>% dplyr::distinct(),
              by = c("Contrast", "Contrast_group", "Contrast_key"))
  list(data = out, full = full, spec = spec)
}

signature_summary <- function(df) {
  df %>% dplyr::group_by(Contrast, Contrast_group) %>%
    dplyr::summarise(
      n_genes = sum(!is.na(log2FC)),
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      median_log2FC = median(log2FC, na.rm = TRUE),
      q25 = quantile(log2FC, .25, na.rm = TRUE), q75 = quantile(log2FC, .75, na.rm = TRUE),
      n_sig_up = sum(!is.na(padj) & padj < .05 & log2FC >= 1),
      n_sig_down = sum(!is.na(padj) & padj < .05 & log2FC <= -1),
      .groups = "drop"
    )
}
