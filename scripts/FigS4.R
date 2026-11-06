# Figure S4: independent paired biventricular validation in GSE133402. Prerequisites: existing Fig1 and Fig5 source data.
# Rebuilds separate-chamber effects, paired interaction, validation and QC.
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
source(file.path(sd,'08_biventricular_helpers.R'))
source(file.path(sd,'09_biventricular_QC.R'))
finalout<-file.path(root,'outputs/FigS4_outputs');dir.create(finalout,showWarnings=FALSE)
pa<-pca_plots[['LV']]+labs(title='A',subtitle='LV: three rats per condition')+scale_x_continuous(expand=expansion(mult=.18))+scale_y_continuous(expand=expansion(mult=.12))+theme(plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
pb<-pca_plots[['RV']]+labs(title='B',subtitle='RV: three rats per condition')+scale_x_continuous(expand=expansion(mult=.18))+scale_y_continuous(expand=expansion(mult=.12))+theme(plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
degcounts<-rbind(data.frame(Chamber=qq$Chamber,Rule='FDR < 0.05',N=qq$FDR05),data.frame(Chamber=qq$Chamber,Rule='FDR < 0.05 and |log2FC| >= 1',N=qq$FDR05_FC2))
pcbar<-ggplot(degcounts,aes(Chamber,N,fill=Rule))+geom_col(position=position_dodge(.75),width=.65)+geom_text(aes(label=N),position=position_dodge(.75),vjust=-.35,size=3)+scale_fill_manual(values=c('FDR < 0.05'='#E68613','FDR < 0.05 and |log2FC| >= 1'='#984EA3'))+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='C',subtitle='Transcriptome-wide differential expression',x=NULL,y='Number of genes',fill=NULL)+theme_manuscript(10)+theme(legend.position='bottom',legend.text=element_text(size=8),plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))+guides(fill=guide_legend(ncol=1))
pd<-pBf+labs(title='D')
pe<-pC+labs(title='E')
pf<-pCf+labs(title='F')
source(file.path(sd,'10_joint_PCA.R'))
figfinal<-(pJoint+(pa+labs(title='B'))+(pb+labs(title='C')))/((pcbar+labs(title='D'))+(pd+labs(title='E'))+pf)
ggsave(file.path(finalout,'FigureS4.pdf'),figfinal,width=16,height=10.5)
ggsave(file.path(finalout,'FigureS4.png'),figfinal,width=16,height=10.5,dpi=300)
# Collect all underlying data without recomputing adjusted P values.
basebook<-file.path(out,'R2_2_source_data.xlsx');fbook<-file.path(out,'FDR_validation/R2_2_FDR_source_data.xlsx');qbook<-file.path(out,'QC/QC_source_data.xlsx')
book<-list()
for(sh in getSheetNames(basebook))book[[sh]]<-read.xlsx(basebook,sheet=sh)
for(sh in c('Validation_summary','FDR_categories','Shared_gene_results'))book[[sh]]<-read.xlsx(fbook,sheet=sh)
for(sh in getSheetNames(qbook))book[[paste0('QC_',sh)]]<-read.xlsx(qbook,sheet=sh)
book[['DEG_plot_counts']]<-degcounts
book[['Joint_PCA_coordinates']]<-jp
book[['Joint_PCA_variance']]<-jvar
book[['Joint_PCA_methods']]<-jmeta
write.xlsx(book,file.path(finalout,'FigureS4_source_data.xlsx'),overwrite=TRUE)
writeLines(c(capture.output(validation),capture.output(qq),capture.output(ftab)),file.path(finalout,'RESULTS.txt'))
message('Figure S4 finished: ',finalout)
