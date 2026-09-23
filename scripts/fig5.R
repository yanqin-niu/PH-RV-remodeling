# ================================================================
# Figure 5 — Shared-derived conserved ventricular remodeling signature
#   Final version: ONLY Fig5A bubble plot
#
# Selection rule:
#   1) Gene must belong to the mutually exclusive Shared set from Figure 1
#   2) Gene must be significant in GSE266139_RV (anchor)
#   3) >=5 of the other 7 RV datasets must be significant in the
#      anchor direction; GSE240921 is included in selection
#   4) Thus, retained genes have 6/8, 7/8, or 8/8 support
#      (anchor + >=5 of the 7 external datasets)
#
# Significance definition:
#   - abs(log2FC) >= 1
#   - padj < 0.05
#
# Output:
#   outputs/Fig5_outputs/
#     - Figure5.pdf
#     - Figure5.png
#     - Fig5A_bubble.pdf
#     - Fig5A_bubble.png
#     - Figure5_source_data.xlsx
#
# ⚠️ NOTE:
# GSE186989 file is not included in the repository due to size limitation.
# Please download it from GEO and place it under:
# data/raw/GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx
# ================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(readr)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(forcats)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
})

set.seed(1)

# ------------------------------------------------
# 0. Auto-detect script directory
# ------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  
  if (length(idx) > 0) {
    path <- sub(file_arg, "", args[idx[1]])
    return(dirname(normalizePath(path)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  
  normalizePath(getwd())
}

script_dir  <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))
raw_dir     <- file.path(project_dir, "data", "raw")
output_root <- file.path(project_dir, "outputs")
out_dir     <- file.path(output_root, "Fig5_outputs")

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------
# 1. Global parameters
# ------------------------------------------------
alpha <- 0.05
lfc_thr <- 1
min_sig_other7 <- 5L
anchor_dataset_name <- "GSE266139_RV (Rat_MCT)"

fig_width  <- 8.8
fig_height <- 7.8

# ------------------------------------------------
# 2. Theme
# ------------------------------------------------
theme_bubble <- function(base_size = 9) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# ------------------------------------------------
# 3. Paths
# ------------------------------------------------
file_266139 <- file.path(raw_dir, "GSE266139_12967_2025_6792_MOESM2_ESM.xlsx")
if (!file.exists(file_266139)) {
  fallback_266139 <- file.path(raw_dir, "12967_2025_6792_MOESM2_ESM.xlsx")
  if (file.exists(fallback_266139)) file_266139 <- fallback_266139
}

file_198618 <- file.path(raw_dir, "GSE198618_Normalized_Counts_RV_ALL.csv")
file_242014 <- file.path(raw_dir, "GSE242014_processed-data-pab.xlsx")
file_240923 <- file.path(raw_dir, "GSE240923_processed-data-mct.xlsx")
file_186989 <- file.path(raw_dir, "GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx")
file_133402 <- file.path(raw_dir, "GSE133402_Rat_RVandLV_HypoxiaNormoxia_Processed.xlsx")
file_240921 <- file.path(raw_dir, "GSE240921_processed-data-human.xlsx")
file_prot   <- file.path(raw_dir, "44161_2022_113_MOESM3_ESM.xlsx")
file_olink  <- file.path(raw_dir, "44161_2022_113_MOESM4_ESM.xlsx")
file_fig1_source <- file.path(output_root, "Fig1_outputs", "Figure1_source_data.xlsx")

stopifnot(
  file.exists(file_266139),
  file.exists(file_198618),
  file.exists(file_242014),
  file.exists(file_240923),
  file.exists(file_186989),
  file.exists(file_133402),
  file.exists(file_240921),
  file.exists(file_prot),
  file.exists(file_olink),
  file.exists(file_fig1_source)
)

olink_sheets <- c("Olink CTRL vs PAH", "Olink CTRL vs dRV PAH")

# ------------------------------------------------
# 4. Helper functions
# ------------------------------------------------
TOUP <- function(x) toupper(trimws(as.character(x)))

clean_sym <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | x == "NA"] <- NA_character_
  x
}

safe_round_counts <- function(mat) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  round(mat, 0)
}

add_sig_col <- function(df, alpha = 0.05, lfc_thr = 1) {
  df %>%
    dplyr::mutate(
      Significant = !is.na(padj) & padj < alpha & abs(log2FC) >= lfc_thr,
      Direction = dplyr::case_when(
        Significant & log2FC > 0  ~ "Up",
        Significant & log2FC < 0  ~ "Down",
        TRUE ~ "NS"
      )
    )
}

map_rat_ensembl_to_symbol <- function(ids) {
  ids2 <- sub("\\..*$", "", ids)
  sy <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = ids2,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(sy[ids2])
}

map_human_ensembl_to_symbol <- function(ids) {
  ids2 <- sub("\\..*$", "", ids)
  sy <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ids2,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(sy[ids2])
}

collapse_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat2 <- mat[keep, , drop = FALSE]
  symbols2 <- TOUP(symbols[keep])
  rownames(mat2) <- symbols2
  
  if (anyDuplicated(rownames(mat2))) {
    mat2 <- rowsum(mat2, group = rownames(mat2)) /
      as.vector(table(rownames(mat2)))
  }
  mat2
}

run_deseq_simple <- function(counts, groups, contrast_vec, is_normalized = FALSE) {
  coldata <- data.frame(condition = factor(groups))
  rownames(coldata) <- colnames(counts)
  
  cmat <- counts
  if (is_normalized) cmat <- round(cmat, 0)
  cmat <- safe_round_counts(cmat)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cmat,
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- DESeq2::results(dds, contrast = contrast_vec)
  df <- as.data.frame(res)
  df$Gene <- rownames(df)
  df <- df[!is.na(df$padj), , drop = FALSE]
  
  df %>%
    dplyr::transmute(
      Gene = TOUP(Gene),
      log2FC = log2FoldChange,
      padj = padj
    ) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
}

split_first <- function(x) {
  sapply(x, function(xx) {
    if (is.na(xx)) return(NA_character_)
    parts <- unlist(strsplit(as.character(xx), "[;,]"))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts)) parts[1] else NA_character_
  })
}

.norm <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

.pick_gene_col <- function(df, prefer_overlap = character(0)) {
  cn <- colnames(df)
  cnn <- .norm(cn)
  patt <- c("gene", "symbol", "genename", "assay", "protein", "name")
  hit <- which(Reduce(`|`, lapply(patt, function(p) grepl(p, cnn))))
  if (length(hit)) return(cn[hit[1]])
  
  if (length(prefer_overlap) > 0) {
    score <- sapply(seq_along(cn), function(i) {
      vals <- TOUP(df[[i]])
      mean(vals %in% prefer_overlap, na.rm = TRUE)
    })
    if (max(score, na.rm = TRUE) > 0.05) return(cn[which.max(score)])
  }
  
  cn[1]
}

# ------------------------------------------------
# 5. GSE266139 -> RV genes only (ANCHOR)
# ------------------------------------------------
read_rv_266139 <- function(file) {
  expr_full <- openxlsx::read.xlsx(file, sheet = 1)
  stopifnot(all(c("Gene.stable.ID", "Gene.name") %in% colnames(expr_full)))
  
  gene_annot <- expr_full[, c("Gene.stable.ID", "Gene.name")]
  colnames(gene_annot) <- c("Geneid", "Genename")
  gene_annot$Geneid   <- clean_sym(gene_annot$Geneid)
  gene_annot$Genename <- TOUP(clean_sym(gene_annot$Genename))
  
  sample_cols <- grep("^(LP|LM|RP|RM)_\\d+$", colnames(expr_full), value = TRUE)
  expr <- expr_full[, sample_cols, drop = FALSE]
  expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))
  
  na_rows <- rowSums(is.na(expr)) > 0
  if (sum(na_rows) > 0) {
    expr <- expr[!na_rows, , drop = FALSE]
    gene_annot <- gene_annot[!na_rows, , drop = FALSE]
  }
  
  expr_mat <- safe_round_counts(expr)
  rownames(expr_mat) <- gene_annot$Geneid
  
  group <- ifelse(grepl("^LP_", colnames(expr_mat)), "LP",
                  ifelse(grepl("^RP_", colnames(expr_mat)), "RP",
                         ifelse(grepl("^LM_", colnames(expr_mat)), "LM",
                                ifelse(grepl("^RM_", colnames(expr_mat)), "RM", NA_character_))))
  stopifnot(!any(is.na(group)))
  
  keep <- rowSums(expr_mat) >= 10
  expr_mat <- expr_mat[keep, , drop = FALSE]
  gene_annot <- gene_annot[match(rownames(expr_mat), gene_annot$Geneid), , drop = FALSE]
  
  coldata <- data.frame(condition = factor(group))
  rownames(coldata) <- colnames(expr_mat)
  
  keep_samp <- coldata$condition %in% c("RM", "RP")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = expr_mat[, keep_samp, drop = FALSE],
    colData = droplevels(coldata[keep_samp, , drop = FALSE]),
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("condition", "RM", "RP"))
  df <- as.data.frame(res)
  df$Geneid <- rownames(df)
  df <- df[!is.na(df$padj), , drop = FALSE]
  df <- dplyr::left_join(df, gene_annot, by = "Geneid")
  
  df %>%
    dplyr::transmute(
      Gene = Genename,
      log2FC = log2FoldChange,
      padj = padj,
      Dataset = anchor_dataset_name
    ) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr) %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

# ------------------------------------------------
# 6. GSE133402 -> RV genes only
# ------------------------------------------------
read_rv_133402 <- function(file) {
  raw <- readxl::read_excel(file, sheet = 1, col_names = FALSE)
  
  hdr <- as.character(unlist(raw[1, ]))
  count_idx <- which(tolower(hdr) == "counts")
  sample_idx <- count_idx - 2
  valid <- sample_idx >= 1
  
  count_idx <- count_idx[valid]
  sample_idx <- sample_idx[valid]
  sample_names <- trimws(hdr[sample_idx])
  
  gene_col_idx <- sample_idx[1]
  row_start <- 2
  
  gene_ids <- as.character(raw[[gene_col_idx]][row_start:nrow(raw)])
  count_mat <- sapply(count_idx, function(j) as.numeric(raw[[j]][row_start:nrow(raw)]))
  count_mat <- as.matrix(count_mat)
  colnames(count_mat) <- sample_names
  rownames(count_mat) <- gene_ids
  
  keep <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
  count_mat <- count_mat[keep, , drop = FALSE]
  
  rn <- rownames(count_mat)
  sy <- if (mean(grepl("^ENS", rn)) > 0.5) map_rat_ensembl_to_symbol(rn) else rn
  mat_sym <- collapse_by_symbol(count_mat, sy)
  
  group <- ifelse(grepl("RMCH", colnames(mat_sym), ignore.case = TRUE), "Hypoxia",
                  ifelse(grepl("RMNT", colnames(mat_sym), ignore.case = TRUE), "Ctrl", NA_character_))
  chamber <- ifelse(grepl("LV", colnames(mat_sym), ignore.case = TRUE), "LV",
                    ifelse(grepl("RV", colnames(mat_sym), ignore.case = TRUE), "RV", NA_character_))
  
  if (any(is.na(group)) || any(is.na(chamber))) {
    stop("GSE133402: sample names could not be parsed")
  }
  
  deg <- run_deseq_simple(
    mat_sym[, chamber == "RV", drop = FALSE],
    group[chamber == "RV"],
    c("condition", "Hypoxia", "Ctrl"),
    is_normalized = FALSE
  )
  deg$Dataset <- "GSE133402_RV (Rat_Hyp)"
  deg
}

# ------------------------------------------------
# 7. Other datasets
# ------------------------------------------------
read_GSE198618 <- function(file) {
  expr_raw <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ens_raw  <- expr_raw[[1]]
  mat      <- as.matrix(expr_raw[, -1, drop = FALSE])
  rownames(mat) <- sub("\\..*$", "", ens_raw)
  mode(mat) <- "numeric"
  
  sy <- map_human_ensembl_to_symbol(rownames(mat))
  mat_sym <- collapse_by_symbol(mat, sy)
  
  samples <- colnames(mat_sym)
  group <- rep(NA_character_, length(samples))
  group[grepl("(?i)(^control|^ctrl$|normal|healthy|rv[-_ ]?normal)", samples, perl = TRUE)] <- "Ctrl"
  group[grepl("(?i)(compensated|^crv$|rv[-_ ]?comp)", samples, perl = TRUE)] <- "cRV"
  group[grepl("(?i)(decomp|^drv$|rv[-_ ]?fail|fail|failure)", samples, perl = TRUE)] <- "dRV"
  if (any(is.na(group))) stop("GSE198618: unmatched sample names")
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = TRUE)
  deg$Dataset <- "GSE198618 (Human)"
  deg
}

read_GSE242014 <- function(file) {
  expr_pab <- openxlsx::read.xlsx(file, sheet = 2)
  ens <- expr_pab[[1]]
  mat <- as.matrix(expr_pab[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- ens
  
  group <- ifelse(grepl("^CTRL", colnames(mat), ignore.case = TRUE), "Ctrl",
                  ifelse(grepl("^Comp", colnames(mat), ignore.case = TRUE), "cRV",
                         ifelse(grepl("^Decomp", colnames(mat), ignore.case = TRUE), "dRV", NA_character_)))
  if (any(is.na(group))) stop("GSE242014: unmatched sample names")
  
  sy <- map_rat_ensembl_to_symbol(rownames(mat))
  mat_sym <- collapse_by_symbol(mat, sy)
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$Dataset <- "GSE242014 (Rat_PAB)"
  deg
}

read_GSE240923_batch1 <- function(file) {
  count1 <- openxlsx::read.xlsx(file, sheet = 2, rowNames = TRUE)
  count1 <- as.matrix(count1)
  mode(count1) <- "numeric"
  
  sy <- map_rat_ensembl_to_symbol(rownames(count1))
  mat_sym <- collapse_by_symbol(count1, sy)
  
  samples <- colnames(mat_sym)
  group <- dplyr::case_when(
    grepl("^CTRL", samples) ~ "Ctrl",
    grepl("^COMP", samples) ~ "cRV",
    grepl("^DECOMP-FP", samples) ~ "dRV",
    grepl("^DECOMP-MCT69", samples) ~ "dRV",
    TRUE ~ NA_character_
  )
  if (any(is.na(group))) stop("GSE240923 batch1: unmatched sample names")
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$Dataset <- "GSE240923_batch1 (Rat_MCT)"
  deg
}

read_GSE240923_batch2 <- function(file) {
  count2 <- openxlsx::read.xlsx(file, sheet = 3, rowNames = TRUE)
  count2 <- as.matrix(count2)
  mode(count2) <- "numeric"
  
  sy <- map_rat_ensembl_to_symbol(rownames(count2))
  mat_sym <- collapse_by_symbol(count2, sy)
  
  samples <- colnames(mat_sym)
  group <- dplyr::case_when(
    grepl("^CTRL", samples, ignore.case = TRUE) ~ "Ctrl",
    grepl("^COMP", samples, ignore.case = TRUE) ~ "cRV",
    grepl("^DECOMP", samples, ignore.case = TRUE) ~ "dRV",
    TRUE ~ NA_character_
  )
  sex <- dplyr::case_when(
    grepl("-f(_|$)", samples, ignore.case = TRUE) ~ "F",
    grepl("-m(_|$)", samples, ignore.case = TRUE) ~ "M",
    TRUE ~ NA_character_
  )
  if (any(is.na(group)) || any(is.na(sex))) stop("GSE240923 batch2: unmatched sample names")
  
  coldata <- data.frame(
    Sex = factor(sex, levels = c("F", "M")),
    Group = factor(group, levels = c("Ctrl", "cRV", "dRV"))
  )
  rownames(coldata) <- colnames(mat_sym)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = safe_round_counts(mat_sym),
    colData = coldata,
    design = ~ Sex + Group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- DESeq2::results(dds, contrast = c("Group", "dRV", "Ctrl"))
  deg <- as.data.frame(res)
  deg$Gene <- rownames(deg)
  deg <- deg[!is.na(deg$padj), , drop = FALSE] %>%
    dplyr::transmute(
      Gene = TOUP(Gene),
      log2FC = log2FoldChange,
      padj = padj
    ) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
  
  deg$Dataset <- "GSE240923_batch2 (Rat_MCT)"
  deg
}

read_GSE186989 <- function(file) {
  df <- readxl::read_excel(file, sheet = 1)
  deg <- df %>%
    dplyr::transmute(
      Gene = TOUP(`Feature ID`),
      log2FC = as.numeric(`EDGE test: CTRL vs SU/HX, tagwise dispersions - Fold change`),
      padj   = as.numeric(`EDGE test: CTRL vs SU/HX, tagwise dispersions - FDR p-value correction`)
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", !is.na(log2FC), !is.na(padj)) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
  
  deg$Dataset <- "GSE186989 (Rat_SuHx)"
  deg
}

read_GSE240921 <- function(file) {
  message("\n========== GSE240921 RAW-DATA DIAGNOSTIC ==========")
  message("Input file: ", normalizePath(file, winslash = "/", mustWork = TRUE))
  message("Workbook sheets: ", paste(readxl::excel_sheets(file), collapse = " | "))

  df <- readxl::read_excel(file, sheet = 2)
  message("Sheet 2 dimensions (rows x columns): ", nrow(df), " x ", ncol(df))
  message("First column name: ", names(df)[1])
  message("Expression columns read from sheet 2:")
  print(names(df)[-1])
  
  expr <- as.matrix(df[, -1, drop = FALSE])
  rownames(expr) <- as.character(df[[1]])
  mode(expr) <- "numeric"
  
  samp <- colnames(expr)

  # GSE240921 contains 40 human RV samples from two sequencing batches.
  # Do NOT retain only columns beginning with "RV-": the first 15 columns
  # (Control_*, Compensated_*, Decompensated_*) are also RV tissue.
  # Primary clinical groups reported by GEO are n=13 normal, n=14
  # compensated, and n=13 decompensated. RV-Failing_5 is clinically
  # classified as compensated RV in the GEO sample metadata.
  group <- dplyr::case_when(
    grepl("^Control_|^RV-Normal_", samp, ignore.case = TRUE) ~ "Ctrl",
    grepl("^Compensated_|^RV-Compen_", samp, ignore.case = TRUE) ~ "cRV",
    samp == "RV-Failing_5" ~ "cRV",
    grepl("^Decompensated_|^RV-Failing_", samp, ignore.case = TRUE) ~ "dRV",
    TRUE ~ NA_character_
  )
  if (any(is.na(group))) stop("GSE240921: unmatched sample names")
  expected_group_n <- c(cRV = 14L, Ctrl = 13L, dRV = 13L)
  observed_group_n <- table(factor(group, levels = names(expected_group_n)))
  if (!identical(as.integer(observed_group_n), as.integer(expected_group_n))) {
    stop(
      "GSE240921: unexpected primary-group sample counts. Expected ",
      "cRV=14, Ctrl=13, dRV=13; observed ",
      paste(names(observed_group_n), as.integer(observed_group_n),
            sep = "=", collapse = ", "), "."
    )
  }
  message("All 40 human RV samples retained and assigned to primary clinical groups:")
  print(data.frame(Sample = samp, Group = group, check.names = FALSE), row.names = FALSE)
  message("Group counts:")
  print(table(group, useNA = "ifany"))

  finite_values <- as.numeric(expr[is.finite(expr)])
  integer_fraction <- if (length(finite_values)) {
    mean(abs(finite_values - round(finite_values)) < 1e-8)
  } else {
    NA_real_
  }
  message("Fraction of finite values that are integer-like: ",
          format(integer_fraction, digits = 5))
  if (length(finite_values)) {
    message("Expression-value range: ",
            paste(format(range(finite_values), digits = 6), collapse = " to "))
  }
  
  rn <- rownames(expr)
  sy <- if (mean(grepl("^ENS", rn)) > 0.5) map_human_ensembl_to_symbol(rn) else rn
  expr_sym <- collapse_by_symbol(expr, sy)

  if ("NPPA" %in% rownames(expr_sym)) {
    message("NPPA values entering DESeq2:")
    print(
      data.frame(
        Sample = colnames(expr_sym),
        Group = group,
        NPPA = as.numeric(expr_sym["NPPA", ]),
        check.names = FALSE
      ),
      row.names = FALSE
    )
  } else {
    message("WARNING: NPPA was not found after ID-to-symbol mapping/collapse.")
    nppa_like <- unique(c(
      grep("NPPA", rownames(expr), value = TRUE, ignore.case = TRUE),
      grep("NPPA", sy, value = TRUE, ignore.case = TRUE)
    ))
    message("NPPA-like raw IDs/symbols: ", paste(nppa_like, collapse = " | "))
  }
  
  deg <- run_deseq_simple(expr_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$Dataset <- "GSE240921 (Human)"
  message("Recomputed dRV vs Ctrl result for NPPA:")
  print(
    deg %>%
      dplyr::filter(Gene == "NPPA") %>%
      dplyr::select(Gene, log2FC, padj, Significant, Direction, Dataset) %>%
      tibble::as_tibble(),
    n = Inf
  )
  message("Expected Table S6 reference: log2FC = 4.19978138; padj = 0.0001986085")
  message("====================================================\n")
  deg
}

# ------------------------------------------------
# 8. Read all data
# ------------------------------------------------
deg_266139_rv  <- read_rv_266139(file_266139)
deg_133402_rv  <- read_rv_133402(file_133402)
deg_198618     <- read_GSE198618(file_198618)
deg_242014     <- read_GSE242014(file_242014)
deg_240923_b1  <- read_GSE240923_batch1(file_240923)
deg_240923_b2  <- read_GSE240923_batch2(file_240923)
deg_186989     <- read_GSE186989(file_186989)
deg_240921     <- read_GSE240921(file_240921)

# Positive-control audit from the independently exported GSE240921 Table S6:
# NPPA must be significantly upregulated in dRV vs Ctrl. Stop rather than
# silently producing inconsistent Figure 5/S1 results if raw parsing changes.
nppa_240921_check <- deg_240921 %>%
  dplyr::filter(Gene == "NPPA")
if (
  nrow(nppa_240921_check) != 1L ||
  !isTRUE(nppa_240921_check$Significant[1]) ||
  nppa_240921_check$Direction[1] != "Up"
) {
  stop(
    "GSE240921 raw-data audit failed: NPPA should be significantly upregulated ",
    "in dRV vs Ctrl (Table S6 reference: log2FC about 4.20, padj about 1.99e-4). ",
    "The diagnostic block printed immediately above identifies the exact file, ",
    "sheet, grouping, NPPA values, and recomputed result. Do not bypass this audit."
  )
}

display_df <- dplyr::bind_rows(
  deg_266139_rv,
  deg_133402_rv,
  deg_198618,
  deg_242014,
  deg_240923_b1,
  deg_240923_b2,
  deg_186989,
  deg_240921
)

other7_df <- dplyr::bind_rows(
  deg_133402_rv,
  deg_198618,
  deg_242014,
  deg_240923_b1,
  deg_240923_b2,
  deg_186989,
  deg_240921
)

dataset_deg_audit <- display_df %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    n_tested = sum(!is.na(padj)),
    n_significant = sum(Significant, na.rm = TRUE),
    n_up = sum(Significant & Direction == "Up", na.rm = TRUE),
    n_down = sum(Significant & Direction == "Down", na.rm = TRUE),
    .groups = "drop"
  )

message("Dataset-level DEG audit:")
print(dataset_deg_audit, n = Inf)
message("GSE240921 NPPA audit from raw-data reanalysis:")
print(
  deg_240921 %>%
    dplyr::filter(Gene == "NPPA") %>%
    dplyr::select(Gene, log2FC, padj, Significant, Direction, Dataset) %>%
    tibble::as_tibble(),
  n = Inf
)

# All seven external datasets, including GSE240921, are retained in
# both the visualization and consensus selection.
# All seven external RV datasets, including GSE240921, enter replication.
other7_informative_df <- other7_df

# Strictly restrict Figure 5 candidates to the mutually exclusive
# Shared set exported by Figure 1.
shared_sheet <- openxlsx::read.xlsx(file_fig1_source, sheet = "Shared_genes")
shared_gene_col <- if ("Genename" %in% names(shared_sheet)) {
  "Genename"
} else if ("Gene" %in% names(shared_sheet)) {
  "Gene"
} else {
  stop("Figure1 Shared_genes sheet has no Genename or Gene column")
}
shared_genes <- unique(TOUP(shared_sheet[[shared_gene_col]]))
shared_genes <- shared_genes[!is.na(shared_genes) & shared_genes != ""]

# ------------------------------------------------
# 9. GSE266139-anchor selection
# ------------------------------------------------
anchor_sig <- deg_266139_rv %>%
  dplyr::filter(Significant, Gene %in% shared_genes) %>%
  dplyr::select(
    Gene,
    anchor_log2FC = log2FC,
    anchor_padj = padj,
    anchor_dir = Direction
  ) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

other7_support <- other7_informative_df %>%
  dplyr::filter(Gene %in% anchor_sig$Gene) %>%
  dplyr::left_join(anchor_sig %>% dplyr::select(Gene, anchor_dir), by = "Gene") %>%
  dplyr::mutate(
    Sig_and_same_as_anchor = Significant & (Direction == anchor_dir)
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    n_other7_sig_same = sum(Sig_and_same_as_anchor, na.rm = TRUE),
    n_other7_up = sum(Significant & Direction == "Up", na.rm = TRUE),
    n_other7_down = sum(Significant & Direction == "Down", na.rm = TRUE),
    .groups = "drop"
  )

consensus_genes <- anchor_sig %>%
  dplyr::left_join(other7_support, by = "Gene") %>%
  dplyr::mutate(
    n_other7_sig_same = dplyr::coalesce(n_other7_sig_same, 0L),
    n_other7_up = dplyr::coalesce(n_other7_up, 0L),
    n_other7_down = dplyr::coalesce(n_other7_down, 0L),
    n_informative_total = 1L + n_other7_sig_same,
    support_class = dplyr::case_when(
      n_informative_total == 8L ~ "8/8 support",
      n_informative_total == 7L ~ "7/8 support",
      n_informative_total == 6L ~ "6/8 support",
      TRUE ~ "Below 6/8"
    ),
    support_label = paste0(Gene, " (", n_informative_total, "/8)"),
    ok = n_other7_sig_same >= min_sig_other7
  ) %>%
  dplyr::filter(ok) %>%
  dplyr::arrange(
    dplyr::desc(n_other7_sig_same),
    dplyr::desc(abs(anchor_log2FC)),
    Gene
  )

cat("Final Shared-derived conserved genes: ", nrow(consensus_genes), "\n", sep = "")
cat(
  "Rule: Shared + anchor significant + >= ", min_sig_other7,
  "/7 external datasets in the same direction (>=6/8 overall)\n",
  sep = ""
)
cat("  8/8 genes: ", sum(consensus_genes$n_informative_total == 8L), "\n", sep = "")
cat("  7/8 genes: ", sum(consensus_genes$n_informative_total == 7L), "\n", sep = "")
cat("  6/8 genes: ", sum(consensus_genes$n_informative_total == 6L), "\n", sep = "")

# ------------------------------------------------
# 10. Plot table
# ------------------------------------------------
plot_long <- display_df %>%
  dplyr::filter(Gene %in% consensus_genes$Gene, Direction != "NS") %>%
  dplyr::distinct(Gene, Dataset, .keep_all = TRUE)

dataset_order <- c(
  "GSE266139_RV (Rat_MCT)",
  "GSE133402_RV (Rat_Hyp)",
  "GSE198618 (Human)",
  "GSE242014 (Rat_PAB)",
  "GSE240923_batch1 (Rat_MCT)",
  "GSE240923_batch2 (Rat_MCT)",
  "GSE186989 (Rat_SuHx)",
  "GSE240921 (Human)"
)

# The first column is the Shared-gene anchor, rather than an ordinary
# independent validation cohort. Keep its data provenance explicit while
# visually distinguishing Figure 5 from the non-anchored all-RV analysis.
dataset_axis_labels <- c(
  "GSE266139_RV (Rat_MCT)" = "Shared genes\n(GSE266139 RV anchor)",
  "GSE133402_RV (Rat_Hyp)" = "GSE133402_RV (Rat_Hyp)",
  "GSE198618 (Human)" = "GSE198618 (Human)",
  "GSE242014 (Rat_PAB)" = "GSE242014 (Rat_PAB)",
  "GSE240923_batch1 (Rat_MCT)" = "GSE240923_batch1 (Rat_MCT)",
  "GSE240923_batch2 (Rat_MCT)" = "GSE240923_batch2 (Rat_MCT)",
  "GSE186989 (Rat_SuHx)" = "GSE186989 (Rat_SuHx)",
  "GSE240921 (Human)" = "GSE240921 (Human)"
)

# ------------------------------------------------
# 11. Protein / Olink support
# ------------------------------------------------
prot_df_raw <- readxl::read_excel(file_prot, sheet = "dRV vs CTRL", col_names = TRUE)
prot_sig_genes <- prot_df_raw %>%
  dplyr::transmute(
    Gene = TOUP(split_first(.data[["Gene.names"]])),
    Status = tolower(trimws(as.character(.data[["Sign_zscore_qval_Limma"]])))
  ) %>%
  dplyr::filter(!is.na(Gene), Gene != "", Status %in% c("up", "down")) %>%
  dplyr::distinct(Gene) %>%
  dplyr::pull(Gene)

prefer_overlap <- consensus_genes$Gene

.read_olink_one <- function(sh) {
  df <- readxl::read_excel(file_olink, sheet = sh, col_names = TRUE)
  gcol <- .pick_gene_col(df, prefer_overlap = prefer_overlap)
  lastc <- tail(colnames(df), 1)
  
  df %>%
    dplyr::transmute(
      Gene = TOUP(trimws(as.character(.data[[gcol]]))),
      status = TOUP(trimws(as.character(.data[[lastc]])))
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", status == "SIGNIFICANT") %>%
    dplyr::distinct(Gene)
}

olink_sig_genes <- dplyr::bind_rows(lapply(olink_sheets, .read_olink_one)) %>%
  dplyr::distinct(Gene) %>%
  dplyr::pull(Gene)

annot_support <- tibble::tibble(Gene = consensus_genes$Gene) %>%
  dplyr::mutate(
    in_prot = Gene %in% prot_sig_genes,
    in_olink = Gene %in% olink_sig_genes,
    support_cat = dplyr::case_when(
      in_prot & in_olink ~ "Protein + Olink + transcript",
      in_prot & !in_olink ~ "Protein + transcript",
      !in_prot & in_olink ~ "Olink + transcript",
      TRUE ~ "Transcript only"
    )
  )

plot_long <- plot_long %>%
  dplyr::left_join(annot_support, by = "Gene")

# ------------------------------------------------
# 12. Gene order
# ------------------------------------------------
gene_stats <- consensus_genes %>%
  dplyr::transmute(
    Gene,
    support_label,
    n_support_other7 = n_other7_sig_same,
    n_informative_total = n_informative_total,
    anchor_abs_lfc = abs(anchor_log2FC),
    anchor_dir = anchor_dir
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_support_other7),
    dplyr::desc(anchor_abs_lfc),
    Gene
  )

gene_order <- gene_stats$support_label

plot_long$Dataset <- factor(plot_long$Dataset, levels = dataset_order)
plot_long <- plot_long %>%
  dplyr::left_join(
    consensus_genes %>% dplyr::select(Gene, support_label, support_class),
    by = "Gene"
  )
plot_long$support_label <- factor(
  plot_long$support_label,
  levels = rev(gene_order)
)

# ------------------------------------------------
# 13. Fig5A — Bubble plot
# ------------------------------------------------
shape_map <- c(
  "Transcript only" = 16,
  "Protein + transcript" = 17,
  "Olink + transcript" = 15,
  "Protein + Olink + transcript" = 18
)

lfc_cap <- as.numeric(stats::quantile(abs(plot_long$log2FC), 0.9, na.rm = TRUE))
if (!is.finite(lfc_cap) || lfc_cap <= 0) lfc_cap <- lfc_thr
plot_long <- plot_long %>%
  dplyr::mutate(lfc_capd = pmin(abs(log2FC), lfc_cap))

pA <- ggplot2::ggplot(plot_long, ggplot2::aes(x = Dataset, y = support_label)) +
  ggplot2::geom_point(
    ggplot2::aes(size = lfc_capd, color = Direction, shape = support_cat),
    alpha = 0.95,
    stroke = 0
  ) +
  ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  ggplot2::scale_x_discrete(labels = dataset_axis_labels) +
  ggplot2::scale_shape_manual(values = shape_map, drop = TRUE) +
  ggplot2::scale_size_continuous(
    name = "|log2FC|",
    range = c(1.1, 3.2),
    limits = c(0, lfc_cap),
    breaks = round(seq(0, lfc_cap, length.out = 3), 2)
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    color = "Direction",
    shape = "External support",
    title = "Shared-gene anchored conserved ventricular remodeling signature"
  ) +
  theme_bubble(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 7.2),
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 8.5)
  )

# ------------------------------------------------
# 14. Save figure
# ------------------------------------------------
ggplot2::ggsave(
  file.path(out_dir, "Figure5.pdf"),
  pA,
  width = fig_width-2,
  height = fig_height,
  units = "in",
  limitsize = FALSE
)

ggplot2::ggsave(
  file.path(out_dir, "Figure5.png"),
  pA,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5A_bubble.pdf"),
  pA,
  width = fig_width,
  height = fig_height,
  units = "in",
  limitsize = FALSE
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5A_bubble.png"),
  pA,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

# ------------------------------------------------
# 15. Export source data
# ------------------------------------------------
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Anchor_GSE266139_RV")
openxlsx::writeData(wb, "Anchor_GSE266139_RV", deg_266139_rv)

openxlsx::addWorksheet(wb, "Other7_df")
openxlsx::writeData(wb, "Other7_df", other7_df)

openxlsx::addWorksheet(wb, "Other7_informative")
openxlsx::writeData(wb, "Other7_informative", other7_informative_df)

openxlsx::addWorksheet(wb, "Shared_input_genes")
openxlsx::writeData(wb, "Shared_input_genes", data.frame(Gene = shared_genes))

openxlsx::addWorksheet(wb, "Display_df_8datasets")
openxlsx::writeData(wb, "Display_df_8datasets", display_df)

openxlsx::addWorksheet(wb, "Dataset_DEG_audit")
openxlsx::writeData(wb, "Dataset_DEG_audit", dataset_deg_audit)

openxlsx::addWorksheet(wb, "Consensus_genes")
openxlsx::writeData(wb, "Consensus_genes", consensus_genes)

openxlsx::addWorksheet(wb, "Fig5A_long")
openxlsx::writeData(wb, "Fig5A_long", plot_long)

openxlsx::addWorksheet(wb, "Protein_support")
openxlsx::writeData(wb, "Protein_support", data.frame(Gene = prot_sig_genes))

openxlsx::addWorksheet(wb, "Olink_support")
openxlsx::writeData(wb, "Olink_support", data.frame(Gene = olink_sig_genes))

openxlsx::addWorksheet(wb, "Meta")
openxlsx::writeData(
  wb, "Meta",
  data.frame(
    Parameter = c(
      "alpha",
      "lfc_thr",
      "anchor_dataset",
      "selection_rule",
      "Figure5_width",
      "Figure5_height"
    ),
    Value = c(
      alpha,
      lfc_thr,
      anchor_dataset_name,
      paste0(
        "Gene must belong to Figure1 Shared genes, be significant in anchor (",
        anchor_dataset_name,
        "), and be significant in the same direction in >= ",
        min_sig_other7,
        "/7 external datasets; all eight datasets enter the analysis"
      ),
      fig_width,
      fig_height
    )
  )
)

openxlsx::addWorksheet(wb, "sessionInfo")
openxlsx::writeData(wb, "sessionInfo", capture.output(sessionInfo()))

openxlsx::saveWorkbook(
  wb,
  file.path(out_dir, "Figure5_source_data.xlsx"),
  overwrite = TRUE
)

message("✅ Figure 5 finished. Outputs saved to: ", out_dir)
