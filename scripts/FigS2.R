# Supplementary Figure S2 — effect-size-matched random resampling analysis
suppressPackageStartupMessages({library(patchwork)})
.args<-commandArgs(trailingOnly=FALSE);.hit<-grep("^--file=",.args);.dir<-if(length(.hit))dirname(normalizePath(sub("^--file=","",.args[.hit[1]])))else if(requireNamespace("rstudioapi",quietly=TRUE)&&nzchar(rstudioapi::getActiveDocumentContext()$path))dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))else normalizePath(getwd())
source(file.path(.dir,"00_common_helpers.R"))
set.seed(20260921);B<-as.integer(Sys.getenv("N_RESAMPLING_ITERATIONS",Sys.getenv("N_PERMUTATIONS","10000")));project_dir<-find_project_dir(get_script_dir());out_dir<-file.path(project_dir,"outputs","FigS2_outputs");dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
f5<-read_fig5_objects(project_dir);sets<-read_fig1_sets(project_dir);grid<-build_replication_grid(sets,f5$display,exclude_anchor=TRUE)

# Five discovery-|log2FC| strata are used. Comparator genes are sampled with replacement
# within the same/nearest occupied stratum so every Shared gene has an effect-size-matched control.
breaks<-unique(quantile(sets$discovery_abs_log2FC,probs=seq(0,1,.2),na.rm=TRUE));if(length(breaks)<3)breaks<-pretty(range(sets$discovery_abs_log2FC),n=5)
sets<-sets%>%dplyr::mutate(bin=cut(discovery_abs_log2FC,breaks=breaks,include.lowest=TRUE,labels=FALSE))
grid<-grid%>%dplyr::select(-discovery_abs_log2FC)%>%dplyr::left_join(sets%>%dplyr::select(Gene,Set,discovery_abs_log2FC,bin),by=c("Gene","Set"))
shared<-grid%>%dplyr::filter(Set=="Shared")
match_sample<-function(pool,target_bins){avail<-sort(unique(pool$bin[!is.na(pool$bin)]));idx<-vapply(target_bins,function(b){bb<-avail[which.min(abs(avail-b))];sample(which(pool$bin==bb),1)},integer(1));pool[idx,,drop=FALSE]}
one_test<-function(ds,comp){s<-shared%>%dplyr::filter(Dataset==ds,measured);pool<-grid%>%dplyr::filter(Dataset==ds,Set==comp,measured);if(!nrow(s)||!nrow(pool))stop("No measured genes for ",ds," / ",comp);obs_s<-mean(s$same_direction,na.rm=TRUE);obs_c<-mean(pool$same_direction,na.rm=TRUE);null_c<-replicate(B,mean(match_sample(pool,s$bin)$same_direction,na.rm=TRUE));tibble::tibble(Dataset=ds,Comparator=comp,Shared_rate=obs_s,Comparator_rate=obs_c,Observed_difference=obs_s-obs_c,Null_mean_difference=mean(obs_s-null_c),P_one_sided=(1+sum((obs_s-null_c)<=0))/(B+1),Null=list(obs_s-null_c))}
res<-dplyr::bind_rows(lapply(unique(grid$Dataset),function(ds)dplyr::bind_rows(one_test(ds,"LV-only"),one_test(ds,"RV-only"))))
plot_res<-res%>%dplyr::mutate(label=paste0(Dataset,"\n",Comparator),Comparator=factor(Comparator,levels=c("LV-only","RV-only")))
pA<-ggplot(plot_res,aes(Observed_difference,factor(label,levels=rev(unique(label))),color=Comparator))+geom_vline(xintercept=0,linetype=2,color="grey50")+geom_point(size=2.6)+scale_color_manual(values=set_cols[c("LV-only","RV-only")])+
 labs(title="A. Observed replication advantage of Shared genes",x="Shared rate - comparator rate",y=NULL,color=NULL)+theme_manuscript(9)+theme(legend.position="top")
# Robust pooled null: average the dataset-specific matched differences at each draw.
pooled<-res%>%dplyr::select(Dataset,Comparator,Null)%>%dplyr::mutate(draw=lapply(Null,seq_along))%>%tidyr::unnest(c(Null,draw))%>%dplyr::group_by(Comparator,draw)%>%dplyr::summarise(Difference=mean(Null),.groups="drop")
obs_pool<-res%>%dplyr::group_by(Comparator)%>%dplyr::summarise(Observed=mean(Observed_difference),.groups="drop")
pB<-ggplot(pooled,aes(Difference,fill=Comparator))+geom_density(alpha=.42,linewidth=.35)+geom_vline(data=obs_pool,aes(xintercept=Observed,color=Comparator),linewidth=.8)+facet_wrap(~Comparator,ncol=1,scales="free_y")+scale_fill_manual(values=set_cols[c("LV-only","RV-only")])+scale_color_manual(values=set_cols[c("LV-only","RV-only")])+
 labs(title="B. Matched reference distributions",subtitle=paste0(B," random resampling iterations; vertical line = observed mean difference"),x="Mean matched difference across datasets",y="Density")+theme_manuscript(10)+theme(legend.position="none")
fig<-pA|pB;ggsave(file.path(out_dir,"FigureS2.pdf"),fig,width=12,height=8);ggsave(file.path(out_dir,"FigureS2.tiff"),fig,width=12,height=8,dpi=600,compression="lzw")
ggsave(file.path(out_dir,"FigS2A_observed_advantage.pdf"),pA,width=6.2,height=7.0)
ggsave(file.path(out_dir,"FigS2B_matched_null.pdf"),pB,width=5.8,height=7.0)
null_export<-res%>%dplyr::select(Dataset,Comparator,Null)%>%tidyr::unnest_longer(Null,values_to="Matched_difference")
openxlsx::write.xlsx(list(Test_summary=res%>%dplyr::select(-Null),Matched_null=null_export,Discovery_sets=sets,Meta=data.frame(N_resampling_iterations=B,Matching="Five quantile strata of discovery absolute log2FC; nearest occupied stratum; sampling with replacement")),file.path(out_dir,"FigureS2_source_data.xlsx"),overwrite=TRUE)
message("Figure S2 finished: ",out_dir)
