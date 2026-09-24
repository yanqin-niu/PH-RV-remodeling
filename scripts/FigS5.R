# Supplementary Figure S5 — full current Figure 5 signature across all CTEPH contrasts
.args<-commandArgs(trailingOnly=FALSE);.hit<-grep("^--file=",.args);.dir<-if(length(.hit))dirname(normalizePath(sub("^--file=","",.args[.hit[1]])))else if(requireNamespace("rstudioapi",quietly=TRUE)&&nzchar(rstudioapi::getActiveDocumentContext()$path))dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))else normalizePath(getwd())
source(file.path(.dir,"00_common_helpers.R"))
project_dir<-find_project_dir(get_script_dir());out_dir<-file.path(project_dir,"outputs","FigS5_outputs");dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
f5<-read_fig5_objects(project_dir); genes<-unique(TOUP(f5$consensus$Gene));ct<-read_cteph_data(project_dir,genes);dat<-ct$data
rank<-dat%>%dplyr::mutate(sig=!is.na(padj)&padj<.05&abs(log2FC)>=1)%>%dplyr::group_by(Gene)%>%dplyr::summarise(n_sig=sum(sig),mean_abs=mean(abs(log2FC),na.rm=TRUE),.groups="drop")%>%dplyr::arrange(dplyr::desc(n_sig),dplyr::desc(mean_abs),Gene)
plot_df<-ct$full%>%dplyr::mutate(Gene=factor(Gene,levels=rev(rank$Gene)),Contrast=factor(Contrast,levels=ct$spec$label),missing_value=is.na(log2FC),
 mark=dplyr::case_when(missing_value~"NA",!is.na(padj)&padj<.001~"***",!is.na(padj)&padj<.01~"**",!is.na(padj)&padj<.05~"*",TRUE~""))
lim<-max(1,quantile(abs(plot_df$log2FC),.98,na.rm=TRUE))
p<-ggplot(plot_df,aes(Contrast,Gene,fill=pmax(-lim,pmin(lim,log2FC))))+geom_tile(color="white",linewidth=.16)+geom_text(aes(label=mark,color=missing_value),size=1.8,show.legend=FALSE)+
 scale_color_manual(values=c("FALSE"="black","TRUE"="grey35"))+scale_fill_gradient2(low="blue",mid="#FFFF7F",high="red",midpoint=0,limits=c(-lim,lim),na.value="grey82")+
 labs(title=paste0("Supplementary Figure S5. CTEPH validation of the current ",length(genes),"-gene signature"),subtitle="* adjusted P < 0.05, ** < 0.01, *** < 0.001; grey/NA = log2FC unavailable",x=NULL,y=NULL,fill="log2FC")+
 theme_manuscript(8)+theme(axis.text.x=element_text(angle=60,hjust=1),panel.border=element_blank(),plot.subtitle=element_text(hjust=.5))
h<-max(7,2.5+.18*length(genes));ggsave(file.path(out_dir,"FigureS5.pdf"),p,width=12,height=h,limitsize=FALSE);ggsave(file.path(out_dir,"FigureS5.tiff"),p,width=12,height=h,dpi=600,compression="lzw",limitsize=FALSE)
ggsave(file.path(out_dir,"FigS5A_full_CTEPH_heatmap.pdf"),p,width=12,height=h,limitsize=FALSE)
write.xlsx(list(Full_data=dat,Full_grid=plot_df,Gene_order=rank,Meta=data.frame(Figure5_source=f5$path)),file.path(out_dir,"FigureS5_source_data.xlsx"),overwrite=TRUE)
message("Figure S5 finished: ",out_dir)
