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

qout<-file.path(out,'QC');dir.create(qout,showWarnings=FALSE)
qc<-list();pcs<-list();cors<-list();sizes<-list();pca_plots<-list()
for(ch in c('LV','RV')){
 sel<-cd$Chamber==ch;d<-DESeq(DESeqDataSetFromMatrix(prefilter_counts(m[,sel],paste0("GSE133402 ",ch)),cd[sel,],~condition),quiet=TRUE)
 r<-as.data.frame(results(d,contrast=c('condition','Hypoxia','Ctrl'),alpha=.05));u<-as.data.frame(results(d,contrast=c('condition','Hypoxia','Ctrl'),alpha=.05,independentFiltering=FALSE))
 sig<-is.finite(r$padj)&r$padj<.05;eligible<-r$baseMean>=5
 qc[[ch]]<-data.frame(Chamber=ch,FDR05=sum(sig),FDR05_FC2=sum(sig&abs(r$log2FoldChange)>=1),FDR05_noIF=sum(u$padj<.05,na.rm=TRUE),P05=sum(r$pvalue<.05,na.rm=TRUE),P_available=sum(is.finite(r$pvalue)),FDR_available=sum(is.finite(r$padj)),Filtered_p_available=sum(is.finite(r$pvalue)&!is.finite(r$padj)),Median_abs_LFC=median(abs(r$log2FoldChange[eligible]),na.rm=TRUE),Median_SE=median(r$lfcSE[eligible],na.rm=TRUE),Median_dispersion=median(dispersions(d)[eligible],na.rm=TRUE))
 vt<-varianceStabilizingTransformation(d,blind=FALSE);a<-assay(vt);a<-a[r$baseMean>=5,,drop=FALSE]
 cc<-cor(a);cors[[ch]]<-data.frame(Sample=rownames(cc),cc,check.names=FALSE)
 pc<-prcomp(t(a),scale.=FALSE);pv<-100*pc$sdev^2/sum(pc$sdev^2)
 pp<-data.frame(Sample=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2]);pp$Condition<-cd[pp$Sample,'condition'];pp$Chamber<-ch;pp$PC1_percent<-pv[1];pp$PC2_percent<-pv[2];pcs[[ch]]<-pp
 sizes[[ch]]<-data.frame(Sample=colnames(d),Chamber=ch,Library_sum=colSums(counts(d)),Size_factor=sizeFactors(d))
 p<-ggplot(pp,aes(PC1,PC2,color=Condition,label=Sample))+geom_point(size=3)+geom_text(aes(label=sub('(LV|RV)2$','',Sample)),vjust=-.8,size=3,show.legend=FALSE)+scale_color_manual(values=c(Ctrl='#1F78B4',Hypoxia='#D73027'))+labs(title=paste(ch,'sample PCA'),x=sprintf('PC1 (%.1f%%)',pv[1]),y=sprintf('PC2 (%.1f%%)',pv[2]))+theme_manuscript(10)
 pca_plots[[ch]]<-p
 ggsave(file.path(qout,paste0(ch,'_PCA.png')),p,width=8,height=5,dpi=150)
}
qq<-do.call(rbind,qc);write.xlsx(c(list(Summary=qq,Sample_sizes=do.call(rbind,sizes),PCA=do.call(rbind,pcs)),cors),file.path(qout,'QC_source_data.xlsx'),overwrite=TRUE)
writeLines(c(capture.output(qq),capture.output(sizes),capture.output(cors)),file.path(qout,'QC_RESULTS.txt'))
