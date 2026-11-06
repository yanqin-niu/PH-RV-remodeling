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

# Joint PCA: common normalization and blind variance-stabilizing transformation.
# No batch correction or chamber/exposure regression is applied.
jout<-file.path(out,'Joint_PCA');dir.create(jout,showWarnings=FALSE)
jd<-DESeqDataSetFromMatrix(m,cd,~condition+Chamber)
jd<-estimateSizeFactors(jd)
jv<-varianceStabilizingTransformation(jd,blind=TRUE)
keep_joint<-rowMeans(counts(jd,normalized=TRUE))>=5
ja<-assay(jv)[keep_joint,,drop=FALSE]
jpc<-prcomp(t(ja),center=TRUE,scale.=FALSE);jpct<-100*jpc$sdev^2/sum(jpc$sdev^2)
jp<-data.frame(Sample=rownames(jpc$x),PC1=jpc$x[,1],PC2=jpc$x[,2],PC3=jpc$x[,3]);jp<-cbind(jp,cd[jp$Sample,c('Animal','Chamber','condition')]);jp$Label<-paste0(jp$Animal,'-',jp$Chamber)
pJoint<-ggplot(jp,aes(PC1,PC2))+geom_line(aes(group=Animal),color='grey72',linewidth=.45)+geom_point(aes(color=condition,shape=Chamber),size=3.4)+ggrepel::geom_text_repel(aes(label=Label,color=condition),size=2.5,seed=20261102,max.overlaps=Inf,box.padding=.45,min.segment.length=0,show.legend=FALSE)+scale_color_manual(values=c(Ctrl='#1F78B4',Hypoxia='#D73027'))+scale_shape_manual(values=c(LV=16,RV=17))+scale_x_continuous(expand=expansion(mult=.2))+scale_y_continuous(expand=expansion(mult=.18))+labs(title='A',subtitle='Paired LV/RV samples: six rats, twelve tissues',x=sprintf('PC1 (%.1f%%)',jpct[1]),y=sprintf('PC2 (%.1f%%)',jpct[2]),color='Exposure',shape='Chamber')+theme_manuscript(10)+theme(plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9),legend.position='bottom')
ggsave(file.path(jout,'Joint_12_samples_PCA.pdf'),pJoint,width=9,height=6)
ggsave(file.path(jout,'Joint_12_samples_PCA.png'),pJoint,width=9,height=6,dpi=200)
jvar<-data.frame(PC=paste0('PC',seq_along(jpct)),Variance_percent=jpct)
jmeta<-data.frame(Item=c('Normalization','Transformation','Genes','Scaling','Pairing','Interpretation'),Value=c('DESeq2 size factors estimated jointly across all twelve samples','Variance stabilizing transformation with blind=TRUE; no exposure/chamber residualization',paste(sum(keep_joint),'genes with joint normalized mean >=5'),'Gene-centered PCA; no unit-variance scaling','Grey lines connect the LV and RV from the same animal','Descriptive structure, not a formal test of relative chamber/exposure effect size. The first two PCs may omit other signals.'))
write.xlsx(list(Coordinates=jp,Variance=jvar,Methods=jmeta,Loadings=data.frame(Gene=rownames(jpc$rotation),jpc$rotation[,1:3])),file.path(jout,'Joint_PCA_source_data.xlsx'),overwrite=TRUE)
writeLines(c(capture.output(jvar[1:3,]),capture.output(jp)),file.path(jout,'RESULTS.txt'))
