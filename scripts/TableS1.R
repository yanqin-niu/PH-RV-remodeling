# Dataset registry and exact preprocessing documentation; no DE fitting.
suppressPackageStartupMessages(library(openxlsx))
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
root<-normalizePath(file.path(sd,'..'));out<-file.path(root,'outputs','TableS1_outputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-read.csv(file.path(root,'data/metadata/dataset_preprocessing_models.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
w<-loadWorkbook(file.path(root,'data/metadata/TableS1_template.xlsx'))
removeWorksheet(w,'Analysis_details');addWorksheet(w,'Analysis_details')
audit<-data.frame(Dataset=x$Dataset,Input_scale=paste(x$Exact_input,x$Scale,sep='; '),Model_filter=paste(x$Model_contrast,x$Prefilter,sep='; '),Samples=x$Samples_exclusions,Limitations=paste(x$Identifiers,x$Limitations,sep='; '))
writeData(w,'Analysis_details',audit);addWorksheet(w,'Preprocessing_R2_5');writeData(w,'Preprocessing_R2_5',x)
for(s in c('Analysis_details','Preprocessing_R2_5')){z<-if(s=='Analysis_details')audit else x;freezePane(w,s,firstRow=TRUE);setColWidths(w,s,cols=1:ncol(z),widths=40);addStyle(w,s,createStyle(wrapText=TRUE,valign='top'),rows=1:(nrow(z)+1),cols=1:ncol(z),gridExpand=TRUE);addStyle(w,s,createStyle(textDecoration='bold',fgFill='#24476A',fontColour='white',wrapText=TRUE),rows=1,cols=1:ncol(z),gridExpand=TRUE,stack=TRUE);setRowHeights(w,s,rows=2:(nrow(z)+1),heights=120);setRowHeights(w,s,rows=1,heights=30)}
# Empty window coordinates inherited from the supplied template are not valid integers.
w$workbook$bookViews <- gsub(' (xWindow|yWindow|windowWidth|windowHeight)=""','',w$workbook$bookViews)
modifyBaseFont(w,fontColour='black')
# Correct the human-cohort role and expose the full sample-level model metadata.
main<-read.xlsx(w,'Table S1',colNames=FALSE,skipEmptyRows=FALSE)
rr<-which(main[[1]]=='GSE240921')
if(length(rr)==1)writeData(w,'Table S1','Human RV cohort for panel derivation; Batch+Sex+Group; sensitivity Figure S8',startRow=rr,startCol=8,colNames=FALSE)
human<-read.csv(file.path(root,'data/metadata/GSE240921_sample_metadata.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE)
addWorksheet(w,'Human_samples');writeData(w,'Human_samples',human[,c('Sample','GEO_title','GSM','Group','Sex','Batch','FASTQ_R1','FASTQ_R2')])
freezePane(w,'Human_samples',firstRow=TRUE);setColWidths(w,'Human_samples',cols=1:8,widths=23)
addStyle(w,'Human_samples',createStyle(wrapText=TRUE,valign='top'),rows=1:41,cols=1:8,gridExpand=TRUE)
addStyle(w,'Human_samples',createStyle(textDecoration='bold',fgFill='#24476A',fontColour='white'),rows=1,cols=1:8,gridExpand=TRUE,stack=TRUE)
saveWorkbook(w,file.path(out,'TableS1_R2_5.xlsx'),overwrite=TRUE)
