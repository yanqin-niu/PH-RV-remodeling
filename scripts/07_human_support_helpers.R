# Human support is descriptive evidence within datasets used for selection.
audit_human_support <- function(signature, results, direction_col, orthology_file) {
  o <- read.delim(orthology_file, comment.char="#", check.names=FALSE, quote="", stringsAsFactors=FALSE)
  pairs <- unique(data.frame(rat=toupper(trimws(o[[1]])), human=toupper(trimws(o[[5]]))))
  pairs <- pairs[nzchar(pairs$rat) & nzchar(pairs$human), ]
  nr <- table(pairs$rat); nh <- table(pairs$human)
  human_ds <- sort(unique(as.character(results$Dataset[grepl("Human", results$Dataset)])))
  if(length(human_ds)!=2L) stop("Expected two human datasets")
  detail <- list(); k <- 0L
  for(i in seq_len(nrow(signature))) {
    g <- toupper(trimws(signature$Gene[i])); expected <- as.character(signature[[direction_col]][i])
    if(!expected %in% c("Up","Down")) stop("Invalid reference direction")
    targets <- pairs$human[pairs$rat==g]
    status <- if(!length(targets)) "Unmapped" else if(length(targets)==1L && nh[targets]==1L) "One-to-one" else "Ambiguous"
    h <- if(status=="One-to-one") targets else NA_character_
    for(ds in human_ds) {
      r <- results[results$Dataset==ds & results$Gene==h & !is.na(results$Gene==h), ,drop=FALSE]
      if(nrow(r)>1L) stop("Duplicate human gene result: ",g," / ",ds)
      lfc <- if(nrow(r)) as.numeric(r$log2FC[1]) else NA_real_
      p <- if(nrow(r)) as.numeric(r$padj[1]) else NA_real_
      tested <- is.finite(lfc) && is.finite(p)
      sig <- tested && p>=0 && p<0.05 && abs(lfc)>=1
      same <- sig && ((expected=="Up" && lfc>0) || (expected=="Down" && lfc<0))
      k<-k+1L;detail[[k]]<-data.frame(Gene=g,Human_ortholog=h,Mapping=status,Dataset=ds,
        Reference_direction=expected,log2FC=lfc,padj=p,Test_available=tested,
        Significant=sig,Same_direction_support=same,
        Evidence=if(!tested) "Unavailable" else if(same) "Significant concordant" else if(sig) "Significant discordant" else "Tested not significant")
    }
  }
  detail<-do.call(rbind,detail)
  gene<-do.call(rbind,lapply(split(detail,detail$Gene),function(z)data.frame(
    Gene=z$Gene[1],Human_ortholog=z$Human_ortholog[1],Mapping=z$Mapping[1],Reference_direction=z$Reference_direction[1],
    N_human_tested=sum(z$Test_available),N_human_supported=sum(z$Same_direction_support),
    Supporting_datasets=paste(z$Dataset[z$Same_direction_support],collapse="; "))))
  summary<-data.frame(Human_datasets_supported=0:2,N_genes=vapply(0:2,function(n)sum(gene$N_human_supported==n),integer(1)))
  summary$N_signature<-nrow(gene);summary$Proportion<-summary$N_genes/nrow(gene)
  list(detail=detail,gene=gene,summary=summary,subset=gene[gene$N_human_supported>=1,,drop=FALSE])
}
