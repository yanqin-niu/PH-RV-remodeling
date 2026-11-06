# Known technical stratum and GEO sex adjustment; GSE240921 remains one vote.
human240921_metadata <- function(samples,project_dir) {
 m<-read.csv(file.path(project_dir,'data/metadata/GSE240921_sample_metadata.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE)
 stopifnot(nrow(m)==40,!anyDuplicated(m$Sample),setequal(m$Sample,samples))
 m<-m[match(samples,m$Sample),];stopifnot(identical(m$Sample,samples),all(m$Sex%in%c('F','M')),all(m$Batch%in%c('PE25','SE15')))
 m$Group<-factor(m$Group);m$Sex<-factor(m$Sex,levels=c('F','M'));m$Batch<-factor(m$Batch,levels=c('PE25','SE15'));rownames(m)<-m$Sample
 stopifnot(all(table(m$Group)==c(cRV=14,Ctrl=13,dRV=13)),all(table(m$Batch)==c(PE25=25,SE15=15)))
 stopifnot(m$Group[m$Sample=='RV-Failing_5']=='cRV')
 mm<-model.matrix(~Batch+Sex+Group,m);stopifnot(qr(mm)$rank==ncol(mm))
 m
}
fit_human240921 <- function(counts,project_dir,out_dir) {
 m<-human240921_metadata(colnames(counts),project_dir)
 counts<-prefilter_counts(round(counts),'GSE240921 all40; Batch + Sex + Group')
 dds<-DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(counts,m,~Batch+Sex+Group),quiet=TRUE)
 z<-as.data.frame(DESeq2::results(dds,contrast=c('Group','dRV','Ctrl')));z$Gene<-rownames(z)
 write.csv(m,file.path(out_dir,'GSE240921_model_samples.csv'),row.names=FALSE)
 write.csv(z,file.path(out_dir,'GSE240921_batch_sex_DESeq2.csv'),row.names=FALSE,na='')
 data.frame(Gene=z$Gene,log2FC=z$log2FoldChange,lfcSE=z$lfcSE,padj=z$padj)
}
