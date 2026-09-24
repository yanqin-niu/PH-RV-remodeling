# Supplementary Figure S4 — leave-one-dataset-out stability of Shared-derived consensus genes
suppressPackageStartupMessages({library(patchwork)})
.args<-commandArgs(trailingOnly=FALSE);.hit<-grep("^--file=",.args);.dir<-if(length(.hit))dirname(normalizePath(sub("^--file=","",.args[.hit[1]])))else if(requireNamespace("rstudioapi",quietly=TRUE)&&nzchar(rstudioapi::getActiveDocumentContext()$path))dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))else normalizePath(getwd())
source(file.path(.dir,"00_common_helpers.R"))
project_dir<-find_project_dir(get_script_dir());out_dir<-file.path(project_dir,"outputs","FigS4_outputs");dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
f5<-read_fig5_objects(project_dir);sets<-read_fig1_sets(project_dir);shared<-sets%>%dplyr::filter(Set=="Shared")%>%dplyr::pull(Gene);d<-f5$display%>%dplyr::filter(Gene%in%shared);datasets<-unique(d$Dataset);final<-unique(TOUP(f5$consensus$Gene))

# The full rule requires 6/8 concordant datasets (75%). After omitting one cohort,
# the proportional analogue is ceiling(0.75 * 7) = 6, therefore LODO uses >=6/7.
# Direction is determined by the majority among retained cohorts; no dataset is forced as anchor.
select_one<-function(left_out){d%>%dplyr::filter(Dataset!=left_out)%>%dplyr::group_by(Gene)%>%dplyr::summarise(n_up=sum(Significant&Direction=="Up"),n_down=sum(Significant&Direction=="Down"),.groups="drop")%>%dplyr::mutate(n_concordant=pmax(n_up,n_down),Direction=dplyr::if_else(n_up>=n_down,"Up","Down"))%>%dplyr::filter(n_concordant>=6)%>%dplyr::mutate(Left_out=left_out)}
lodo<-dplyr::bind_rows(lapply(datasets,select_one));selected_sets<-split(lodo$Gene,lodo$Left_out)
run_summary<-tibble::tibble(Left_out=datasets)%>%dplyr::rowwise()%>%dplyr::mutate(n_selected=length(unique(selected_sets[[Left_out]])),n_overlap=length(intersect(final,selected_sets[[Left_out]])),n_union=length(union(final,selected_sets[[Left_out]])),Jaccard=n_overlap/n_union)%>%dplyr::ungroup()
gene_stability<-tibble::tibble(Gene=sort(unique(c(final,lodo$Gene))))%>%dplyr::mutate(n_LODO=vapply(Gene,function(g)sum(vapply(selected_sets,function(x)g%in%x,logical(1))),integer(1)),In_final=Gene%in%final)%>%dplyr::arrange(dplyr::desc(In_final),dplyr::desc(n_LODO),Gene)
pA<-ggplot(run_summary,aes(Jaccard,factor(Left_out,levels=rev(datasets))))+geom_segment(aes(x=0,xend=Jaccard,y=Left_out,yend=Left_out),color="grey70",linewidth=.7)+geom_point(size=2.8,color="#D73027")+scale_x_continuous(limits=c(0,1))+
 labs(title="A. Agreement with the full Figure 5 set",x="Jaccard index",y="Dataset omitted")+theme_manuscript(10)
pB<-ggplot(gene_stability%>%dplyr::filter(In_final),aes(n_LODO,reorder(Gene,n_LODO)))+geom_col(fill="#D73027",width=.72)+scale_x_continuous(breaks=0:length(datasets),limits=c(0,length(datasets)))+
 labs(title="B. Stability of final consensus genes",x="Selected in LODO analyses",y=NULL)+theme_manuscript(9)
fig<-pA|pB;h<-max(6,2.2+.14*nrow(dplyr::filter(gene_stability,In_final)));ggsave(file.path(out_dir,"FigureS4.pdf"),fig,width=11,height=h,limitsize=FALSE);ggsave(file.path(out_dir,"FigureS4.tiff"),fig,width=11,height=h,dpi=600,compression="lzw",limitsize=FALSE)
ggsave(file.path(out_dir,"FigS4A_Jaccard_by_omitted_dataset.pdf"),pA,width=5.4,height=4.8)
ggsave(file.path(out_dir,"FigS4B_gene_stability.pdf"),pB,width=5.4,height=h)
write.xlsx(list(LODO_selected=lodo,Run_summary=run_summary,Gene_stability=gene_stability,Final_Fig5_genes=data.frame(Gene=final),Meta=data.frame(Rule=">=6/7 retained datasets significant in the same majority direction; no dataset is forced as anchor")),file.path(out_dir,"FigureS4_source_data.xlsx"),overwrite=TRUE)
message("Figure S4 finished: ",out_dir)
