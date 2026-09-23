# Figure 4 — adjusted cross-dataset replication of discovery gene sets
suppressPackageStartupMessages({library(patchwork)})
.args <- commandArgs(trailingOnly = FALSE); .hit <- grep("^--file=", .args)
.dir <- if (length(.hit)) dirname(normalizePath(sub("^--file=", "", .args[.hit[1]]))) else if (requireNamespace("rstudioapi", quietly=TRUE) && nzchar(rstudioapi::getActiveDocumentContext()$path)) dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)) else normalizePath(getwd())
source(file.path(.dir, "00_common_helpers.R"))

script_dir <- get_script_dir(); project_dir <- find_project_dir(script_dir)
out_dir <- file.path(project_dir, "outputs", "Fig4_outputs"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig1_sets <- read_fig1_sets(project_dir); f5 <- read_fig5_objects(project_dir)
grid <- build_replication_grid(fig1_sets, f5$display, exclude_anchor = TRUE)
summary_df <- summarise_replication(grid)
dataset_order <- unique(grid$Dataset)

pA <- ggplot(summary_df, aes(color = Set)) +
  geom_segment(aes(x=lower,xend=upper,y=as.numeric(factor(Dataset, levels=rev(dataset_order))) + c("Shared"=-.2,"LV-only"=0,"RV-only"=.2)[as.character(Set)],
                   yend=as.numeric(factor(Dataset, levels=rev(dataset_order))) + c("Shared"=-.2,"LV-only"=0,"RV-only"=.2)[as.character(Set)]), linewidth=.45) +
  geom_point(aes(x=rate,y=as.numeric(factor(Dataset, levels=rev(dataset_order))) + c("Shared"=-.2,"LV-only"=0,"RV-only"=.2)[as.character(Set)]), size = 2.4) +
  scale_y_continuous(breaks=seq_along(rev(dataset_order)), labels=rev(dataset_order)) +
  scale_color_manual(values = set_cols) + scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "A. Same-direction replication by dataset", x = "Replication proportion (95% CI)", y = NULL, color = "Gene set") +
  theme_manuscript(10) + theme(legend.position = "top")

gene_support <- grid %>% group_by(Gene, Set, discovery_abs_log2FC) %>%
  summarise(n_same = sum(same_direction, na.rm = TRUE), n_tested = sum(measured), .groups = "drop")
pB <- ggplot(gene_support, aes(Set, n_same, fill = Set)) +
  geom_violin(scale = "width", trim = FALSE, alpha = .35, color = NA) +
  geom_boxplot(width = .18, outlier.shape = NA, linewidth = .4) +
  geom_jitter(width = .12, size = .55, alpha = .35) +
  scale_fill_manual(values = set_cols) + scale_y_continuous(breaks = 0:7, limits = c(0, 7)) +
  labs(title = "B. Per-gene replication breadth", x = NULL, y = "Datasets with same-direction replication") +
  theme_manuscript(10) + theme(legend.position = "none")

model_df <- grid %>% filter(measured) %>% mutate(Set = relevel(factor(Set), ref = "Shared"),
  z_abs_lfc = as.numeric(scale(discovery_abs_log2FC)), same_direction = as.integer(same_direction))
if (requireNamespace("lme4", quietly=TRUE)) {
  fit <- lme4::glmer(same_direction ~ Set + z_abs_lfc + (1 | Dataset) + (1 | Gene), data = model_df,
                     family = binomial(), control = lme4::glmerControl(optimizer = "bobyqa"))
  cf <- lme4::fixef(fit); se <- sqrt(diag(as.matrix(vcov(fit)))); model_name <- "Crossed random-intercept logistic model"
} else {
  warning("Package lme4 is unavailable; using dataset-adjusted logistic regression fallback.")
  fit <- glm(same_direction ~ Set + z_abs_lfc + Dataset, data=model_df, family=binomial())
  cf <- coef(fit); se <- sqrt(diag(vcov(fit))); model_name <- "Dataset-adjusted logistic regression"
}
or_df <- tibble(term=names(cf), estimate=exp(cf), conf.low=exp(cf-1.96*se), conf.high=exp(cf+1.96*se), Model=model_name) %>%
  filter(term %in% c("SetLV-only", "SetRV-only", "z_abs_lfc")) %>%
  mutate(term = recode(term, "SetLV-only" = "LV-only vs Shared", "SetRV-only" = "RV-only vs Shared",
                       "z_abs_lfc" = "Discovery |log2FC| (per SD)"))
pC_cols <- c(
  "Discovery |log2FC| (per SD)" = "#E68613",
  "LV-only vs Shared" = "#1F78B4",
  "RV-only vs Shared" = "#33A02C"
)
pC <- ggplot(or_df, aes(estimate, reorder(term, estimate), color=term)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
  geom_segment(
    aes(x=conf.low,xend=conf.high,y=reorder(term,estimate),yend=reorder(term,estimate)),
    linewidth=1.05,
    lineend="round"
  ) +
  geom_point(aes(fill=term), shape=21, size=4.0, stroke=.65, color="white") +
  scale_color_manual(values=pC_cols) + scale_x_log10() +
  scale_fill_manual(values=pC_cols) +
  labs(title = "C. Pooled external replication model (7 datasets)",
       subtitle="Dataset and gene random intercepts; adjusted for discovery effect size",
       x = "Odds ratio (95% CI; log scale)", y = NULL) +
  theme_manuscript(10) + theme(legend.position="none",plot.subtitle=element_text(hjust=.5,size=8.5))

species_df <- grid %>% filter(measured) %>% group_by(Species, Set) %>% summarise(n_same = sum(same_direction), n_total = n(), .groups = "drop") %>%
  rowwise() %>% mutate(rate = n_same/n_total, ci = list(wilson_ci(n_same,n_total)), lower=ci[[1]], upper=ci[[2]]) %>% ungroup()
pD <- ggplot(species_df, aes(Species, rate, fill = Set)) +
  geom_col(position = position_dodge(.75), width = .68, color = "black", linewidth = .25) +
  geom_errorbar(aes(ymin=lower,ymax=upper), position=position_dodge(.75), width=.18, linewidth=.4) +
  scale_fill_manual(values=set_cols) + scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
  labs(title="D. Replication stratified by species", x=NULL, y="Same-direction replication", fill="Gene set") +
  theme_manuscript(10) + theme(legend.position="top")

fig4 <- (pA | pB) / (pC | pD)
ggsave(file.path(out_dir,"Figure4.pdf"),fig4,width=11,height=8.2)
ggsave(file.path(out_dir,"Figure4.tiff"),fig4,width=11,height=8.2,dpi=600,compression="lzw")

# Individual panels: edit width/height here when assembling the final figure.
ggsave(file.path(out_dir,"Fig4A_replication_by_dataset.pdf"),pA,width=4.5,height=4)
ggsave(file.path(out_dir,"Fig4B_per_gene_replication.pdf"),pB,width=3.8,height=3.3)
ggsave(file.path(out_dir,"Fig4C_adjusted_model.pdf"),pC,width=4.3,height=3.2)
ggsave(file.path(out_dir,"Fig4D_species_summary.pdf"),pD,width=4.0,height=3.35)
write.xlsx(list(Replication_grid=grid, Dataset_summary=summary_df, Gene_support=gene_support,
                Mixed_model_OR=or_df, Species_summary=species_df,
                Meta=data.frame(Figure5_source=f5$path, Note="GSE266139 excluded; denominator is genes with an estimable external result; same-direction padj<0.05 and |log2FC|>=1")),
           file.path(out_dir,"Figure4_source_data.xlsx"), overwrite=TRUE)
message("Figure 4 finished: ", out_dir)
