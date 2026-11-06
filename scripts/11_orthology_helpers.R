# Curated reciprocal-unique mapping. Native human results are kept separately.
read_curated_orthology <- function(project_dir) {
 o<-read.delim(file.path(project_dir,'data/reference/ORTHOLOGS_RAT.txt'),comment.char='#',check.names=FALSE,quote='')
 p<-unique(data.frame(rat=toupper(trimws(o[[1]])),human=toupper(trimws(o[[5]]))))
 p<-p[!is.na(p$rat)&!is.na(p$human)&nzchar(p$rat)&nzchar(p$human),]
 nr<-table(p$rat);nh<-table(p$human);good<-p[nr[p$rat]==1 & nh[p$human]==1,]
 list(pairs=p,good=good)
}
rat_query <- function(genes,project_dir) {
 o<-read_curated_orthology(project_dir);q<-o$good$human[match(genes,o$good$rat)]
 data.frame(Gene=genes,Query=q,Mapping=ifelse(!is.na(q),'One_to_one',ifelse(genes%in%o$pairs$rat,'Ambiguous','Unmapped')))
}
harmonize_native_results <- function(x,project_dir) {
 o<-read_curated_orthology(project_dir);x$NativeGene<-x$Gene;h<-grepl('Human',x$Dataset)
 x$Gene[h]<-o$good$rat[match(x$NativeGene[h],o$good$human)]
 x$Mapping<-ifelse(h,'One_to_one','Same_species')
 x<-x[!is.na(x$Gene)&nzchar(x$Gene),]
 if(anyDuplicated(paste(x$Gene,x$Dataset)))stop('Duplicate mapped gene-dataset result')
 x
}
