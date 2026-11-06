# R2.2 exploratory paired LV/RV external validation; leaves numbered figures unchanged.
suppressPackageStartupMessages({library(DESeq2);library(readxl);library(openxlsx);library(ggplot2);library(patchwork)})
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
source(file.path(sd,'00_common_helpers.R'));root<-find_project_dir(sd)
out<-file.path(root,'outputs/R2_2_biventricular');dir.create(out,recursive=TRUE,showWarnings=FALSE)
f<-file.path(root,'data/raw/GSE133402_Rat_RVandLV_HypoxiaNormoxia_Processed.xlsx')
x<-read_excel(f,col_names=FALSE);h<-as.character(unlist(x[1,]));ci<-which(tolower(h)=='counts');si<-ci-2
samples<-trimws(h[si]);ids<-lapply(si,function(j)as.character(x[[j]][-1]));stopifnot(all(vapply(ids,function(z)identical(z,ids[[1]]),logical(1))))
m<-sapply(ci,function(j)as.numeric(x[[j]][-1]));colnames(m)<-samples
genes<-toupper(trimws(ids[[1]]));k<-!is.na(genes)&nzchar(genes);m<-m[k,,drop=FALSE];genes<-genes[k]
stopifnot(all(is.finite(m)),all(m>=0));fractional<-mean(abs(m-round(m))>1e-6)
m<-rowsum(m,genes,reorder=FALSE);m<-prefilter_counts(round(m),"GSE133402 all12")
cd<-data.frame(Sample=samples,Animal=sub('(LV|RV)2$','',samples),Chamber=ifelse(grepl('LV',samples),'LV','RV'),condition=factor(ifelse(grepl('RMCH',samples),'Hypoxia','Ctrl'),levels=c('Ctrl','Hypoxia')),row.names=samples)
stopifnot(nrow(cd)==12,all(table(cd$Animal,cd$Chamber)==1),all(table(cd$condition,cd$Chamber)==3))
cd$GSM<-paste0('GSM',3908245:3908256)
expected<-c('RMCH1LV2','RMCH1RV2','RMCH2LV2','RMCH2RV2','RMCH4LV2','RMCH4RV2','RMNT1LV2','RMNT1RV2','RMNT2LV2','RMNT2RV2','RMNT3LV2','RMNT3RV2')
stopifnot(identical(samples,expected))
fit<-function(ch){s<-cd$Chamber==ch;d<-DESeq(DESeqDataSetFromMatrix(prefilter_counts(m[,s],paste0("GSE133402 ",ch)),cd[s,],~condition),quiet=TRUE);r<-as.data.frame(results(d,contrast=c('condition','Hypoxia','Ctrl'),alpha=.05));r$Gene<-rownames(r);r}
lv<-fit('LV');rv<-fit('RV');v<-merge(lv,rv,by='Gene',suffixes=c('_LV','_RV'))
# Descriptive effect comparisons require expression >=5 normalized counts in BOTH chambers.
v$eligible<-v$baseMean_LV>=5 & v$baseMean_RV>=5 & is.finite(v$log2FoldChange_LV)&is.finite(v$log2FoldChange_RV)
disc<-read.xlsx(file.path(root,'outputs/Fig1_outputs/Figure1_source_data.xlsx'),sheet='LV_RV_merge');disc<-disc[disc$set=='Shared',];disc$Gene<-toupper(disc$Genename)
disc<-disc[order(disc$Gene,disc$padj_rv,-abs(disc$lfc_rv),disc$Geneid),];dups<-disc[duplicated(disc$Gene)|duplicated(disc$Gene,fromLast=TRUE),];disc<-disc[!duplicated(disc$Gene),]
v$Discovery_direction<-sign(disc$lfc_rv[match(v$Gene,disc$Gene)]);v$Shared<-v$Gene%in%disc$Gene
v$LV_sig<-!is.na(v$padj_LV)&v$padj_LV<.05&abs(v$log2FoldChange_LV)>=1
v$RV_sig<-!is.na(v$padj_RV)&v$padj_RV<.05&abs(v$log2FoldChange_RV)>=1
v$Direction_class<-ifelse(!v$Shared|!v$eligible,NA,ifelse(sign(v$log2FoldChange_LV)==v$Discovery_direction & sign(v$log2FoldChange_RV)==v$Discovery_direction,'Both match discovery',ifelse(sign(v$log2FoldChange_LV)==-v$Discovery_direction & sign(v$log2FoldChange_RV)==-v$Discovery_direction,'Both oppose discovery','Mixed / zero')))
# Animal fixed effects absorb exposure main effect. RVHyp tests difference-in-differences.
cd$Animal<-factor(cd$Animal);cd$RV<-as.numeric(cd$Chamber=='RV');cd$RVHyp<-cd$RV*as.numeric(cd$condition=='Hypoxia')
mm<-model.matrix(~Animal+RV+RVHyp,cd);stopifnot(qr(mm)$rank==ncol(mm))
joint<-DESeq(DESeqDataSetFromMatrix(m,cd,~Animal+RV+RVHyp),quiet=TRUE)
inter<-as.data.frame(results(joint,name='RVHyp',alpha=.05));inter$Gene<-rownames(inter)
# Cross-check existing Fig5 RV effects: same input, rounding and contrast.
old<-read.xlsx(file.path(root,'outputs/Fig5_outputs/Figure5_source_data.xlsx'),sheet='Display_df_8datasets');old<-old[grepl('133402',old$Dataset),]
ix<-match(old$Gene,rv$Gene);delta<-max(abs(old$log2FC-rv$log2FoldChange[ix]),na.rm=TRUE);stopifnot(is.finite(delta),delta<1e-5)
stat<-function(z,label){z<-z[z$eligible,];data.frame(Subset=label,N=nrow(z),Spearman=cor(z$log2FoldChange_LV,z$log2FoldChange_RV,method='spearman'),Same_direction=sum(sign(z$log2FoldChange_LV)==sign(z$log2FoldChange_RV)),LV_DEG=sum(z$LV_sig),RV_DEG=sum(z$RV_sig),Both_DEG=sum(z$LV_sig&z$RV_sig),Both_DEG_same_direction=sum(z$LV_sig&z$RV_sig&sign(z$log2FoldChange_LV)==sign(z$log2FoldChange_RV)))}
stats<-rbind(stat(v,'All expression-eligible'),stat(v[v$Shared,],'Discovery Shared'))
s<-v[v$Shared&v$eligible,];lev<-c('Both match discovery','Mixed / zero','Both oppose discovery');tab<-data.frame(Class=factor(lev,levels=lev),N=as.numeric(table(factor(s$Direction_class,levels=lev))))
coverage<-data.frame(Gene=disc$Gene,Discovery_direction=sign(disc$lfc_rv),Present=disc$Gene%in%v$Gene,Eligible=v$eligible[match(disc$Gene,v$Gene)])
cols<-c('Both match discovery'='#D73027','Mixed / zero'='#E68613','Both oppose discovery'='#1F78B4')
base<-function(z,title,subtitle){ggplot(z,aes(log2FoldChange_LV,log2FoldChange_RV))+geom_hline(yintercept=0,color='grey65',linewidth=.3)+geom_vline(xintercept=0,color='grey65',linewidth=.3)+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+labs(title=title,subtitle=subtitle,x='LV log2 fold change (Hypoxia / Normoxia)',y='RV log2 fold change (Hypoxia / Normoxia)')+theme_manuscript(10)+theme(plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))}
pA<-base(v[v$eligible,],'A',sprintf('All eligible genes: n = %s; Spearman rho = %.2f',stats$N[1],stats$Spearman[1]))+geom_point(color='#4DAF4A',alpha=.18,size=.65)
pB<-base(s,'B',sprintf('Discovery Shared: n = %s; Spearman rho = %.2f',nrow(s),stats$Spearman[2]))+geom_point(aes(color=Direction_class),size=2,alpha=.8)+scale_color_manual(values=cols,drop=FALSE)+theme(legend.position='none')
pC<-ggplot(tab,aes(Class,N,fill=Class))+geom_col(width=.65)+geom_text(aes(label=N),vjust=-.4,size=4)+scale_fill_manual(values=cols,drop=FALSE)+scale_x_discrete(labels=c('Both match\ndiscovery','Mixed / zero','Both oppose\ndiscovery'))+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='C',subtitle='Direction relative to MCT discovery',x=NULL,y='Number of Shared genes')+theme_manuscript(10)+theme(legend.position='none',plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
fig<-pA+pB+pC+plot_layout(widths=c(1,1,.85));ggsave(file.path(out,'R2_2_biventricular.pdf'),fig,width=15,height=4.8);ggsave(file.path(out,'R2_2_biventricular.png'),fig,width=15,height=4.8,dpi=300)
meta<-data.frame(Item=c('Input','Rounding','Pairing','Contrast','Eligibility','Threshold','Interaction','Inference','Existing RV agreement'),Value=c('Published Kallisto estimated gene counts; NOT integer raw read counts. No transcript effective lengths available. Rounded DESeq2 reanalysis follows existing pipeline; exploratory inference.',paste('Fraction noninteger before symbol collapse:',fractional),'Six animal IDs each have LV and RV; original Methods: six rats, twelve RNA samples; Ctrl n=3 and Hypoxia n=3','Hypoxia / Normoxia in each chamber; positive = higher under hypoxia','baseMean >=5 in both separate-chamber models; finite effects in both; no significance preselection','padj <0.05 and absolute log2FC >=1; default independent filtering retained','~Animal + RV + RVHyp; RVHyp = RV hypoxia effect minus LV hypoxia effect; animal absorbs exposure main effect','Gene-level concordance summaries are descriptive; genes are not independent biological replicates. No binomial enrichment claim. Nonsignificant interaction is not equivalence.',paste('Maximum absolute log2FC deviation:',delta)))
write.xlsx(list(Samples=cd,LV_DESeq2=lv,RV_DESeq2=rv,Paired_interaction=inter,All_gene_effects=v,Shared_coverage=coverage,Shared_directions=s,Direction_counts=tab,Summary=stats,Duplicate_discovery_symbols=dups,Methods=meta,Session=data.frame(Text=capture.output(sessionInfo()))),file.path(out,'R2_2_source_data.xlsx'),overwrite=TRUE)
write.csv(stats,file.path(out,'Summary.csv'),row.names=FALSE);write.csv(s,file.path(out,'Shared_directions.csv'),row.names=FALSE)
writeLines(c(capture.output(stats),capture.output(tab),capture.output(table(coverage$Eligible,useNA='always')),paste('Interaction FDR <0.05:',sum(inter$padj<.05,na.rm=TRUE)),capture.output(meta)),file.path(out,'RESULTS.txt'))
message('Finished: ',out)

# Alternative external validation: transcriptome-wide FDR, no minimum effect magnitude.
# Original stringent outputs above are preserved as a sensitivity analysis.
fdrout<-file.path(out,'FDR_validation');dir.create(fdrout,showWarnings=FALSE)
s$LV_FDR_evaluable<-is.finite(s$padj_LV)
s$RV_FDR_evaluable<-is.finite(s$padj_RV)
s$LV_direction_support<-sign(s$log2FoldChange_LV)==s$Discovery_direction
s$RV_direction_support<-sign(s$log2FoldChange_RV)==s$Discovery_direction
s$LV_FDR_support<-s$LV_FDR_evaluable & s$padj_LV<.05 & s$LV_direction_support
s$RV_FDR_support<-s$RV_FDR_evaluable & s$padj_RV<.05 & s$RV_direction_support
s$FDR_class<-ifelse(!s$LV_FDR_evaluable|!s$RV_FDR_evaluable,'FDR unavailable',ifelse(s$LV_FDR_support & s$RV_FDR_support,'Both chambers',ifelse(s$LV_FDR_support,'LV only',ifelse(s$RV_FDR_support,'RV only','Neither chamber'))))
flev<-c('Both chambers','LV only','RV only','Neither chamber','FDR unavailable')
s$FDR_class<-factor(s$FDR_class,levels=flev)
ftab<-data.frame(Class=factor(flev,levels=flev),N=as.numeric(table(s$FDR_class)))
fcols<-c('Both chambers'='#D73027','LV only'='#1F78B4','RV only'='#33A02C','Neither chamber'='#B3B3B3','FDR unavailable'='#984EA3')
pBf<-base(s,'B',sprintf('Discovery Shared: n = %s; Spearman rho = %.2f',nrow(s),stats$Spearman[2]))+geom_point(aes(color=FDR_class),size=2,alpha=.85)+scale_color_manual(values=fcols,drop=FALSE,labels=c('Significant in both','Significant in LV alone','Significant in RV alone','Neither significant','FDR unavailable'))+labs(color='FDR < 0.05 + discovery direction')+theme(legend.position='bottom',legend.text=element_text(size=8),legend.title=element_text(size=9))+guides(color=guide_legend(nrow=2,byrow=TRUE))
pCf<-ggplot(ftab,aes(Class,N,fill=Class))+geom_col(width=.65)+geom_text(aes(label=N),vjust=-.4,size=4)+scale_fill_manual(values=fcols,drop=FALSE)+scale_x_discrete(labels=c('Both','LV alone','RV alone','Neither','FDR\nunavailable'))+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='C',subtitle='FDR < 0.05 + discovery direction; no FC cutoff',x=NULL,y='Number of Shared genes')+theme_manuscript(10)+theme(legend.position='none',plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
ffig<-pA+pBf+pCf+plot_layout(widths=c(1,1.1,1.05))
ggsave(file.path(fdrout,'R2_2_FDR_validation.pdf'),ffig,width=16,height=5.5)
ggsave(file.path(fdrout,'R2_2_FDR_validation.png'),ffig,width=16,height=5.5,dpi=300)
validation<-data.frame(Metric=c('Discovery Shared total','Present in external table','Expression-eligible in both','FDR evaluable in both among eligible','Both direction match (no FDR threshold)','LV FDR support','RV FDR support','Both FDR support','Both FDR support plus abs(log2FC)>=1','FDR significant opposite direction in LV','FDR significant opposite direction in RV'),N=c(nrow(disc),sum(coverage$Present),nrow(s),sum(s$LV_FDR_evaluable&s$RV_FDR_evaluable),sum(s$LV_direction_support&s$RV_direction_support),sum(s$LV_FDR_support),sum(s$RV_FDR_support),sum(s$LV_FDR_support&s$RV_FDR_support),sum(s$LV_FDR_support&s$RV_FDR_support&s$LV_sig&s$RV_sig),sum(s$LV_FDR_evaluable&s$padj_LV<.05&!s$LV_direction_support),sum(s$RV_FDR_evaluable&s$padj_RV<.05&!s$RV_direction_support)))
fmeta<-rbind(meta,data.frame(Item='FDR validation rule',Value='Primary external support: original transcriptome-wide DESeq2 padj <0.05 and discovery-consistent direction; NO minimum absolute FC. FDR is not recalculated within the 108 Shared genes. Expression eligibility remains unchanged. Missing padj in either chamber is explicitly separated from Neither chamber. Direction and strict-threshold sensitivity outputs are retained.'))
write.xlsx(list(Validation_summary=validation,FDR_categories=ftab,Shared_gene_results=s,Shared_coverage=coverage,LV_full=lv,RV_full=rv,Samples=cd,Methods=fmeta),file.path(fdrout,'R2_2_FDR_source_data.xlsx'),overwrite=TRUE)
write.csv(s,file.path(fdrout,'Shared_FDR_validation.csv'),row.names=FALSE)
writeLines(c(capture.output(validation),capture.output(ftab)),file.path(fdrout,'RESULTS.txt'))
message('FDR validation finished: ',fdrout)
