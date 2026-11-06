# Consensus stability conditional on a frozen discovery Shared candidate universe.
consensus_stability <- function(project_dir) {
 f5<-read_fig5_objects(project_dir);sets<-read_fig1_sets(project_dir)
 shared<-sets[sets$Set=='Shared',];genes<-unique(shared$Gene);ds<-unique(f5$display$Dataset)
 stopifnot(length(ds)==8L)
 g<-tidyr::expand_grid(Gene=genes,Dataset=ds)
 g<-dplyr::left_join(g,f5$display,by=c('Gene','Dataset'))
 g$Evaluable<-is.finite(g$log2FC)&is.finite(g$padj)&g$padj>=0&g$padj<=1
 g$Reference_direction<-shared$discovery_direction[match(g$Gene,shared$Gene)]
 anchor<-ds[grepl('GSE266139',ds)];stopifnot(length(anchor)==1)
 primary<-unique(f5$consensus$Gene)
 select<-function(lfc,required,omit=NULL,force_anchor=FALSE){
  z<-g[!g$Dataset%in%omit,];z$up<-z$Evaluable&z$padj<.05&z$log2FC>=lfc;z$down<-z$Evaluable&z$padj<.05&z$log2FC<=-lfc
  x<-z%>%dplyr::group_by(Gene)%>%dplyr::summarise(N_evaluable=sum(Evaluable),N_up=sum(up),N_down=sum(down),.groups='drop')%>%dplyr::mutate(N_support=pmax(N_up,N_down),Selected_direction=ifelse(N_up>=N_down,'Up','Down'))
  if(force_anchor){a<-g[g$Dataset==anchor,];a$pass<-a$Evaluable&a$padj<.05&abs(a$log2FC)>=lfc;x<-x[x$Gene%in%a$Gene[a$pass],];x<-x[x$Selected_direction==a$Reference_direction[match(x$Gene,a$Gene)],]}
  x[x$N_support>=required,,drop=FALSE]
 }
 base<-select(1,6,force_anchor=TRUE);stopifnot(setequal(base$Gene,primary))
 lodo<-do.call(rbind,lapply(ds,function(omit){x<-select(1,6,omit);x$Left_out<-omit;x}))
 summary<-do.call(rbind,lapply(ds,function(omit){selected<-lodo$Gene[lodo$Left_out==omit];lost<-setdiff(primary,selected);added<-setdiff(selected,primary);data.frame(Left_out=omit,n_selected=length(selected),n_overlap=length(intersect(primary,selected)),n_union=length(union(primary,selected)),Jaccard=length(intersect(primary,selected))/length(union(primary,selected)),N_lost=length(lost),N_added=length(added),Lost_genes=paste(sort(lost),collapse=';'),Added_genes=paste(sort(added),collapse=';'),Shared_membership_fixed=TRUE,Anchor_vote_removed=omit==anchor)}))
 full<-select(1,0);ret<-data.frame(Gene=genes,In_final=genes%in%primary,n_LODO=vapply(genes,function(x)sum(lodo$Gene==x),integer(1)))
 ret$Full_support<-full$N_support[match(ret$Gene,full$Gene)]
 ret$Expected_by_threshold<-ifelse(ret$Full_support>=7,8L,ifelse(ret$Full_support==6,2L,0L))
 stopifnot(all(ret$n_LODO[ret$In_final]==ret$Expected_by_threshold[ret$In_final]))
 # With 8 baseline votes and a 6-of-7 omission threshold, exactly 6/8 support
 # yields retention only when one of the two nonsupporting datasets is omitted.
 # At least 7/8 support survives every single omission. Lock these expected
 # counts for the current frozen 35-gene panel so output changes fail visibly.
 stopifnot(
   sum(ret$In_final & ret$Full_support==6L)==21L,
   sum(ret$In_final & ret$Full_support>=7L)==14L,
   sum(ret$In_final & ret$n_LODO==2L)==21L,
   sum(ret$In_final & ret$n_LODO==8L)==14L
 )
 ret<-ret[order(!ret$In_final,-ret$n_LODO,ret$Gene),]
 thresholds<-list();calls<-list();idx<-0
 for(lfc in c(.585,1,1.5))for(req in 5:7){
  idx<-idx+1;x<-select(lfc,req,force_anchor=TRUE);calls[[idx]]<-data.frame(Gene=genes,External_LFC=lfc,Support_required=req,Selected=genes%in%x$Gene)
  thresholds[[idx]]<-data.frame(LFC=lfc,Support_required=req,Candidate_n=length(genes),Selected=nrow(x),Retained_primary=sum(primary%in%x$Gene),Lost_primary=paste(sort(setdiff(primary,x$Gene)),collapse=';'),Added=paste(sort(setdiff(x$Gene,primary)),collapse=';'),Jaccard=length(intersect(primary,x$Gene))/length(union(primary,x$Gene)))
 }
 meta<-data.frame(Item=c('Candidate universe','Primary rule','LODO','Anchor omission','Interpretation','Threshold grid','Meta-analysis'),Value=c('Frozen 108 discovery Shared genes; membership is not recalculated after omission. Discovery information remains part of this conditional analysis.','Full panel requires discovery RV support and >=6/8 significant concordant votes.','Omit each of the eight datasets once and require >=6/7 same-direction votes. The 108-gene Shared universe remains fixed; no omitted dataset is forced as an anchor.','When GSE266139 RV is omitted, only its support vote is removed. Shared eligibility remains fixed; this is not discovery-independent validation.','Among the final 35 genes, 21 have exactly 6/8 support and survive exactly 2/8 omissions by construction; the other 14 have >=7/8 support and survive all eight. This is threshold arithmetic, not independent evidence of biological stability.','LFC 0.585/1/1.5 by support 5/6/7 out of 8; fixed Shared universe; threshold also applied to discovery RV anchor. No new discovery candidates admitted.','External rat studies only; discovery and human inputs excluded; one GSE240923 batch per fit. Direct fitted log2FC SE required; never reconstructed from adjusted P.'))
 list(grid=g,summary=summary,lodo=lodo,retention=ret,thresholds=do.call(rbind,thresholds),threshold_calls=do.call(rbind,calls),methods=meta,primary=primary)
}
# Intercept-only DerSimonian-Laird random-effects model with conservative modified
# Hartung-Knapp t interval. Base R only; no silently missing package fallback.
meta_dl <- function(y,se){
 k<-length(y);stopifnot(k>=3,all(is.finite(y)),all(is.finite(se)&se>0))
 v<-se^2;w<-1/v;mu<-sum(w*y)/sum(w);Q<-sum(w*(y-mu)^2);C<-sum(w)-sum(w^2)/sum(w)
 tau2<-max(0,(Q-(k-1))/C);wr<-1/(v+tau2);b<-sum(wr*y)/sum(wr)
 q<-sum(wr*(y-b)^2)/(k-1);sb<-sqrt(max(1,q)/sum(wr));t<-qt(.975,k-1)
 c(K=k,Estimate=b,SE=sb,Lower=b-t*sb,Upper=b+t*sb,P=2*pt(-abs(b/sb),k-1),Tau2=tau2,Q=Q,Q_P=pchisq(Q,k-1,lower.tail=FALSE),I2=if(Q>0)max(0,(Q-k+1)/Q)*100 else 0)
}
# Meaningful deterministic checks: homogeneous effects and effect-sign invariance.
stopifnot(abs(meta_dl(rep(1,4),rep(.2,4))['Estimate']-1)<1e-12,meta_dl(rep(1,4),rep(.2,4))['Tau2']==0)
stopifnot(abs(meta_dl(c(1,2,3),c(.2,.3,.4))['Estimate']+meta_dl(c(-1,-2,-3),c(.2,.3,.4))['Estimate'])<1e-12)
external_effect_synthesis <- function(project_dir,st){
 x<-read_fig5_objects(project_dir)$display
 if(!'lfcSE'%in%names(x)){
  p<-resolve_fig5_source(project_dir);x<-openxlsx::read.xlsx(p,sheet='Display_df_8datasets')
 }
 if(!'lfcSE'%in%names(x))stop('Run updated Fig5.R to export fitted lfcSE; no P-value imputation allowed')
 scenarios<-list(Primary_batch1=c('GSE133402','GSE186989','GSE242014','GSE240923_batch1'),Alternate_batch2=c('GSE133402','GSE186989','GSE242014','GSE240923_batch2'),Exclude_estimated_Hypoxia=c('GSE186989','GSE242014','GSE240923_batch1'))
 results<-list();detail<-list();k<-0L
 for(scenario in names(scenarios)){
  patterns<-scenarios[[scenario]];ds<-unique(x$Dataset);chosen<-ds[vapply(ds,function(d)any(vapply(patterns,function(p)grepl(p,d,fixed=TRUE),logical(1))),logical(1))]
  stopifnot(length(chosen)==length(patterns),!any(grepl('Human|GSE266139',chosen)))
  for(gene in unique(st$grid$Gene)){
   z<-x[x$Gene==gene&x$Dataset%in%chosen,];good<-is.finite(z$log2FC)&is.finite(z$lfcSE)&z$lfcSE>0
   ref<-unique(st$grid$Reference_direction[st$grid$Gene==gene]);sign<-if(ref=='Up')1 else -1
   k<-k+1;fit<-if(sum(good)>=3)meta_dl(sign*z$log2FC[good],z$lfcSE[good]) else setNames(rep(NA_real_,10),c('K','Estimate','SE','Lower','Upper','P','Tau2','Q','Q_P','I2'))
   fit['K']<-sum(good)
   results[[k]]<-data.frame(Gene=gene,In_primary=gene%in%st$primary,Scenario=scenario,as.list(fit),N_expected=length(chosen),Estimable=sum(good)>=3,Studied_datasets=paste(chosen,collapse=';'),Reference_direction=ref)
   if(nrow(z))detail[[k]]<-data.frame(Gene=gene,Scenario=scenario,Dataset=z$Dataset,log2FC=z$log2FC,lfcSE=z$lfcSE,Estimable=good,Direction_aligned=sign*z$log2FC)
  }
 }
 z<-do.call(rbind,results);z$FDR_within_scenario<-NA_real_
 for(sc in names(scenarios)){i<-z$Scenario==sc & z$Estimable;z$FDR_within_scenario[i]<-p.adjust(z$P[i],method='BH')}
 list(summary=z,inputs=do.call(rbind,detail))
}

stopifnot(abs(meta_dl(c(0,2,4),rep(1,3))['Tau2']-3)<1e-12,abs(meta_dl(c(0,2,4),rep(1,3))['Estimate']-2)<1e-12)
