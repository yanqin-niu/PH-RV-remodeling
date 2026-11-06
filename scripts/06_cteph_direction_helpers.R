# Directional assessment against the RV discovery effect, not the LV/RV mean.
cteph_direction_audit <- function(full_grid, consensus, alpha = 0.05, lfc_thr = 1) {
  required <- c("Gene", "anchor_log2FC")
  if (!all(required %in% names(consensus))) stop("Consensus lacks discovery RV effect")
  ref <- data.frame(Gene = toupper(trimws(consensus$Gene)),
                    discovery_RV_log2FC = as.numeric(consensus$anchor_log2FC))
  if (anyDuplicated(ref$Gene) || any(!is.finite(ref$discovery_RV_log2FC)) ||
      any(ref$discovery_RV_log2FC == 0)) stop("Invalid discovery direction reference")
  x <- as.data.frame(full_grid)
  x <- merge(x, ref, by = "Gene", all.x = TRUE, sort = FALSE)
  if (any(!is.finite(x$discovery_RV_log2FC))) stop("Missing discovery RV direction")
  if (anyDuplicated(x[c("Gene", "Contrast")])) stop("Duplicate gene/contrast")
  x$effect_available <- is.finite(x$log2FC)
  x$test_available <- x$effect_available & is.finite(x$padj)
  x$direction_product <- x$log2FC * sign(x$discovery_RV_log2FC)
  x$direction_status <- ifelse(!x$effect_available, "Unavailable",
    ifelse(x$direction_product > 0, "Concordant",
      ifelse(x$direction_product < 0, "Discordant", "Zero effect")))
  x$significant <- x$test_available & x$padj >= 0 & x$padj < alpha & abs(x$log2FC) >= lfc_thr
  summary <- do.call(rbind, lapply(split(x, x$Contrast), function(z) {
    n <- sum(z$effect_available)
    data.frame(Contrast = z$Contrast[1], n_signature = nrow(z), n_available = n,
      n_unavailable = sum(!z$effect_available), n_test_available = sum(z$test_available),
      n_positive = sum(z$log2FC > 0, na.rm = TRUE),
      median_log2FC = if(n) median(z$log2FC[z$effect_available]) else NA_real_,
      n_concordant = sum(z$direction_status == "Concordant"),
      n_discordant = sum(z$direction_status == "Discordant"),
      n_zero = sum(z$direction_status == "Zero effect"),
      concordance_rate = if(n) sum(z$direction_status == "Concordant") / n else NA_real_,
      n_sig_concordant = sum(z$significant & z$direction_status == "Concordant"),
      n_sig_discordant = sum(z$significant & z$direction_status == "Discordant"))
  }))
  rownames(summary) <- NULL
  list(gene_level = x, summary = summary, by_expected_direction = summarize_cteph_by_direction(x),
       definitions = data.frame(Definition = c(
         "Reference: GSE266139 RV disease-versus-control anchor_log2FC from Figure5 Consensus_genes",
         "Concordant: sign of CTEPH effect matches discovery RV sign; zero effects are separate",
         "Concordance denominator: finite CTEPH effect; missing effects are unavailable",
         "Significant concordance: finite adjusted P >=0 and <0.05 and absolute log2FC >=1",
         "Regional and post-PEA direction comparisons are descriptive and not disease replication")))
}
# Stratify descriptive CTEPH evidence by the discovery-expected direction.
summarize_cteph_by_direction <- function(x) {
 stopifnot(all(c('Gene','Contrast','discovery_RV_log2FC','effect_available','test_available','direction_status','significant')%in%names(x)))
 x$Expected_direction<-ifelse(x$discovery_RV_log2FC>0,'Up','Down')
 do.call(rbind,lapply(split(x,list(x$Contrast,x$Expected_direction),drop=TRUE),function(z){
  n<-sum(z$effect_available);nt<-sum(z$test_available)
  data.frame(Contrast=z$Contrast[1],Expected_direction=z$Expected_direction[1],N_signature=nrow(z),N_effect_available=n,N_test_available=nt,N_missing=nrow(z)-n,N_direction_concordant=sum(z$direction_status=='Concordant'),Direction_rate=if(n)sum(z$direction_status=='Concordant')/n else NA_real_,N_strict_concordant=sum(z$significant & z$direction_status=='Concordant'),N_strict_discordant=sum(z$significant & z$direction_status=='Discordant'),Concordant_genes=paste(z$Gene[z$direction_status=='Concordant'],collapse=';'),Discordant_genes=paste(z$Gene[z$direction_status=='Discordant'],collapse=';'))
 }))
}
