# Cross-modality annotations only: never used in transcriptomic selection.
audit_protein_support <- function(panel, direction_col, project_dir, out_dir) {
  if(!exists('rat_query',mode='function'))source(file.path(project_dir,'scripts/11_orthology_helpers.R'))
  query <- rat_query(panel$Gene, project_dir)
  query$Expected_direction <- panel[[direction_col]]
  rawdir <- file.path(project_dir,'data/raw')
  p <- readxl::read_excel(file.path(rawdir,'44161_2022_113_MOESM3_ESM.xlsx'),sheet='dRV vs CTRL')
  protein <- data.frame(Modality='Tissue',Contrast='dRV minus CTRL',
    Assay=as.character(p[['Gene.names']]),Assay_ID=as.character(p[['Protein.IDs']]),
    Effect=as.numeric(p[['log2Ratio']]),Adjusted_P=as.numeric(p[['qval_Limma']]),
    Z=as.numeric(p[['zscore']]),Source_call=as.character(p[['Sign_zscore_qval_Limma']]))
  protein$Query <- toupper(trimws(protein$Assay))
  # Multi-gene protein groups are retained in the assay audit but not assigned to one gene.
  protein$Mapping_rule <- ifelse(grepl('[;, ]',protein$Query),'Excluded_multi_gene','Exact_symbol')
  protein$Query[protein$Mapping_rule=='Excluded_multi_gene'] <- NA_character_
  protein$Significant <- is.finite(protein$Adjusted_P)&protein$Adjusted_P<.05&is.finite(protein$Z)&abs(protein$Z)>=2
  stopifnot(all(protein$Significant == (protein$Source_call %in% c('Up','Down'))))
  sheets <- c('Olink CTRL vs PAH','Olink CTRL vs dRV PAH')
  plasma <- do.call(rbind,lapply(sheets,function(sh) {
    x<-readxl::read_excel(file.path(rawdir,'44161_2022_113_MOESM4_ESM.xlsx'),sheet=sh,skip=1)
    x<-x[!is.na(x$Assay)&nzchar(trimws(x$Assay)),]
    disease<-if(sh==sheets[1])'PAH' else 'Decompensated'
    num<-function(z)suppressWarnings(as.numeric(z))
    effect<-num(x[[disease]])-num(x$Control)
    estimate<-num(x$estimate)
    stopifnot(all(abs(effect+estimate)<1e-8,na.rm=TRUE))
    q<-toupper(trimws(x$Assay));rule<-rep('Exact_symbol',length(q))
    # Source sheet 3B identifies NTproBNP as the NPPB-derived assay.
    rule[q=='NTPROBNP']<-'Explicit_NPPB_product';q[q=='NTPROBNP']<-'NPPB'
    a<-data.frame(Modality='Olink',Contrast=paste(disease,'minus CTRL'),Assay=x$Assay,
      Assay_ID=x$OlinkID,Effect=effect,Adjusted_P=num(x$Adjusted_pval),Z=NA_real_,
      Source_call=x$Threshold,Query=q,Mapping_rule=rule)
    a$Significant<-is.finite(a$Adjusted_P)&a$Adjusted_P<.05
    stopifnot(all(a$Significant==(a$Source_call=='Significant'),na.rm=TRUE))
    a
  }))
  assays<-rbind(protein,plasma)
  assays$Direction<-ifelse(assays$Effect>0,'Up',ifelse(assays$Effect<0,'Down','Zero'))
  long<-merge(query,assays,by='Query',all.x=TRUE,sort=FALSE)
  # Complete each panel gene x comparison grid, including assays absent from the input.
  comparisons<-unique(assays[c('Modality','Contrast')])
  grid<-merge(query,comparisons,by=NULL)
  long<-merge(grid,assays[!is.na(assays$Query),],by=c('Query','Modality','Contrast'),all.x=TRUE,sort=FALSE)
  long$Covered<-!is.na(long$Assay)
  long$Evaluable<-long$Covered&is.finite(long$Effect)&is.finite(long$Adjusted_P)
  long$Significant[is.na(long$Significant)]<-FALSE
  long$Concordant<-long$Evaluable&long$Direction==long$Expected_direction
  long$Concordant[is.na(long$Concordant)]<-FALSE
  long$Significant_concordant<-long$Significant&long$Concordant
  long$Significant_discordant<-long$Significant&long$Evaluable&!long$Concordant
  groups<-split(long,interaction(long$Gene,long$Modality,long$Contrast,drop=TRUE))
  gene<-do.call(rbind,lapply(groups,function(x) {
    y<-x[1,c('Gene','Query','Mapping','Expected_direction','Modality','Contrast')]
    y$N_assays<-sum(x$Covered);y$Covered<-any(x$Covered);y$Evaluable<-any(x$Evaluable)
    y$N_significant_concordant<-sum(x$Significant_concordant)
    y$N_significant_discordant<-sum(x$Significant_discordant)
    y$Conflicting_significant_assays<-y$N_significant_concordant>0 & y$N_significant_discordant>0
    y$Directional_support<-y$N_significant_concordant>0 & !y$Conflicting_significant_assays
    y
  }))
  annot<-data.frame(Gene=panel$Gene)
  supported<-function(mod)unique(gene$Gene[gene$Modality==mod & gene$Directional_support])
  annot$in_prot<-annot$Gene %in% supported('Tissue')
  annot$in_olink<-annot$Gene %in% supported('Olink')
  annot$support_cat<-ifelse(annot$in_prot & annot$in_olink,'Protein + Olink + transcript',
    ifelse(annot$in_prot,'Protein + transcript',ifelse(annot$in_olink,'Olink + transcript','Transcript only')))
  summary<-do.call(rbind,lapply(split(gene,interaction(gene$Modality,gene$Contrast,drop=TRUE)),function(x)
    data.frame(Modality=x$Modality[1],Contrast=x$Contrast[1],Panel_n=nrow(panel),Covered=sum(x$Covered),
      Evaluable=sum(x$Evaluable),Supported=sum(x$Directional_support),
      Significant_discordant=sum(x$N_significant_discordant>0),Conflicts=sum(x$Conflicting_significant_assays))))
  methods<-data.frame(Item=c('Source','Tissue','Olink','Mapping','Aggregation','Independence'),
    Definition=c('Boucherat et al 2022; DOI 10.1038/s44161-022-00113-w; imported statistical results, not refitted.',
    'MOESM3 dRV vs CTRL; 10 dRV and 10 CTRL tissue columns; supplied log2Ratio=dRV/CTRL; qval_Limma<0.05 and abs(zscore)>=2; no transcriptomic twofold cutoff.',
    'MOESM4 two sheets; actual header row 2; supplied adjusted P<0.05 from Welch tests; disease-minus-control NPX difference, opposite to supplied estimate; no extra FC cutoff or additional log transform. Adjustment method not re-inferred from summary data.',
    'Frozen reciprocal-unique rat-human RGD mapping, then exact human symbol; NTproBNP explicitly NPPB-derived; multi-gene tissue groups excluded; assay IDs retained.',
    'At least one significant concordant assay and no significant opposite assay within a contrast. Olink plot annotation is union of two contrast-specific support calls; contrasts are related, not independent replications. All assay rows retained.',
    'Tissue proteomics and GSE198618 originate from the same source study; participant crosswalk to transcriptomic cohorts and plasma not available in the files used. Overlap cannot be excluded. Cross-modality annotation, not independent-subject validation. Plasma cannot establish RV tissue origin.'))
  result<-list(Assays=assays,Gene_assay=long,Gene_contrast=gene,Summary=summary,Annotation=annot,Methods=methods)
  for(nm in names(result))write.csv(result[[nm]],file.path(out_dir,paste0('Protein_support_',nm,'.csv')),row.names=FALSE,na='')
  legend<-paste0('Bubble color and size show transcriptomic direction and absolute log2 fold change. Shapes annotate significant direction-consistent protein evidence, not independent-subject validation. Tissue: imported dRV-versus-CTRL log2Ratio, qval_Limma<0.05 and |zscore|>=2. Olink: supplied adjusted P<0.05 in PAH-versus-CTRL or decompensated-PAH-versus-CTRL, oriented as disease-minus-control NPX. No protein twofold cutoff is imposed. Significant opposite assays within a contrast prevent a support call. The Olink shape is the union of two related comparisons; their results are separately exported. Reference direction is discovery RV for Figure 5 and RV-wide consensus for Figure S5. Transcript only means no qualifying directional protein evidence, including unmeasured/unmapped genes, not a negative protein test. Frozen one-to-one orthology and exact human symbols are used; ambiguous tissue protein groups are excluded. Protein evidence does not enter transcriptomic selection. Participant overlap with transcriptomic cohorts cannot be excluded from available files. See Protein_Gene_assay, Protein_Gene_contrast, Protein_Summary and Protein_Methods in source data.')
  writeLines(legend,file.path(out_dir,'Protein_support_legend.txt'))
  result
}
add_protein_audit <- function(wb,audit) {
  for(nm in names(audit)) {sheet<-paste0('Protein_',nm);openxlsx::addWorksheet(wb,sheet);openxlsx::writeData(wb,sheet,audit[[nm]])}
}
