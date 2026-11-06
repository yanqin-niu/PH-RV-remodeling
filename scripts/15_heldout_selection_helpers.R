# Fixed-discovery, leave-external-source-out reselection and held-out assessment.
heldout_panel_assessment <- function(project_dir) {
 f5<-read_fig5_objects(project_dir);sets<-read_fig1_sets(project_dir);sets<-sets[sets$Set=='Shared',]
 g<-build_replication_grid(sets,f5$native,exclude_anchor=FALSE)
 ds<-unique(g$Dataset);anchor<-ds[grepl('GSE266139',ds)];external<-setdiff(ds,anchor)
 stopifnot(length(ds)==8,length(anchor)==1,length(external)==7)
 accession<-function(x)sub('.*(GSE[0-9]+).*','\\1',x)
 publication<-function(x){a<-accession(x);ifelse(a%in%c('GSE240921','GSE240923','GSE242014'),'Boucherat_39196250',a)}
 units<-data.frame(Dataset=ds,Accession=accession(ds),Publication=publication(ds),Role=ifelse(ds==anchor,'Fixed discovery','External derivation cohort'),Input_scale=ifelse(grepl('GSE198618',ds),'Normalized-count approximation',ifelse(grepl('GSE133402',ds),'Estimated-count approximation','Deposited read counts')))
 candidates<-unique(sets$Gene);primary<-unique(f5$consensus$Gene)
 select_training<-function(training){
  stopifnot(anchor%in%training$Dataset)
  ntrain<-length(unique(training$Dataset));req<-ceiling(.75*ntrain)
  z<-training%>%dplyr::group_by(Gene)%>%dplyr::summarise(N_train_evaluable=sum(measured),N_train_support=sum(same_direction),.groups='drop')
  a<-training[training$Dataset==anchor,]
  z$Anchor_pass<-a$same_direction[match(z$Gene,a$Gene)]
  z$Train_units<-ntrain;z$Required_support<-req
  z$Selected<-z$Anchor_pass & z$N_train_evaluable>=req & z$N_train_support>=req
  z
 }
 baseline<-select_training(g);stopifnot(setequal(baseline$Gene[baseline$Selected],primary))
 definitions<-list(Publication_no_normalized=split(external[!grepl('GSE198618',external)],publication(external[!grepl('GSE198618',external)])),Publication_block=split(external,publication(external)),Accession_block=split(external,accession(external)),Cohort_secondary=as.list(setNames(external,external)))
 folds<-list();calls<-list();evaluations<-list();metrics<-list();i<-0L
 eval_metrics<-function(z,expected){
  z<-if(expected=='All')z else z[z$discovery_direction==expected,]
  n<-nrow(z);ea<-is.finite(z$log2FC);ta<-z$measured
  align<-z$log2FC*ifelse(z$discovery_direction=='Up',1,-1)
  same<-ea & align>0;sig_fdr<-ta & z$padj<.05 & same;strict<-z$same_direction
  data.frame(Expected_direction=expected,N_selected=n,N_effect_available=sum(ea),N_test_available=sum(ta),N_unavailable=n-sum(ta),Coverage_effect=if(n)sum(ea)/n else NA_real_,Coverage_test=if(n)sum(ta)/n else NA_real_,N_direction_concordant=sum(same),Direction_rate=if(sum(ea))sum(same)/sum(ea) else NA_real_,N_FDR_concordant=sum(sig_fdr),FDR_rate=if(sum(ta))sum(sig_fdr)/sum(ta) else NA_real_,N_strict_concordant=sum(strict),Strict_rate=if(sum(ta))sum(strict)/sum(ta) else NA_real_,Median_aligned_log2FC=if(sum(ea))median(align[ea]) else NA_real_)
 }
 for(scheme in names(definitions))for(unit in names(definitions[[scheme]])){
  i<-i+1L;held<-definitions[[scheme]][[unit]];active_ds<-if(scheme=='Publication_no_normalized')ds[!grepl('GSE198618',ds)] else ds;train_ds<-setdiff(active_ds,held)
  train<-g[g$Dataset%in%train_ds,];stopifnot(!any(train$Dataset%in%held))
  chosen<-select_training(train);selected<-chosen$Gene[chosen$Selected]
  # Leakage check: perturb every held-out outcome; the training selection must not change.
  perturbed<-g;ix<-perturbed$Dataset%in%held;perturbed$same_direction[ix]<-!perturbed$same_direction[ix];perturbed$log2FC[ix]<-1e6
  check<-select_training(perturbed[perturbed$Dataset%in%train_ds,]);stopifnot(identical(chosen,check))
  fid<-paste(scheme,unit,sep=' / ')
  folds[[i]]<-data.frame(Fold=fid,Scheme=scheme,Heldout_unit=unit,Heldout_datasets=paste(held,collapse=';'),Training_datasets=paste(train_ds,collapse=';'),N_training=length(train_ds),Required_support=ceiling(.75*length(train_ds)),N_selected=length(selected),Normalized_cohort_excluded=scheme=='Publication_no_normalized',Overlap_final35=sum(selected%in%primary),Candidates_fixed=length(candidates),Heldout_used_for_selection=FALSE)
  calls[[i]]<-cbind(Fold=fid,Scheme=scheme,chosen)
  ev<-g[g$Dataset%in%held & g$Gene%in%selected,];ev$Fold<-fid;ev$Scheme<-scheme;ev$In_full_panel<-ev$Gene%in%primary;evaluations[[i]]<-ev
  mm<-do.call(rbind,lapply(held,function(d){z<-ev[ev$Dataset==d,];rbind(cbind(Dataset=d,eval_metrics(z,'All')),cbind(Dataset=d,eval_metrics(z,'Up')),cbind(Dataset=d,eval_metrics(z,'Down')))}))
  metrics[[i]]<-cbind(Fold=fid,Scheme=scheme,mm)
 }
 folds<-do.call(rbind,folds);calls<-do.call(rbind,calls);ev<-do.call(rbind,evaluations);mm<-do.call(rbind,metrics)
 # Equal weighting of held-out cohorts within a source block, then equal block weights.
 per_fold<-mm[mm$Expected_direction=='All',]%>%dplyr::group_by(Fold,Scheme)%>%dplyr::summarise(N_heldout_cohorts=dplyr::n(),Direction_rate=mean(Direction_rate,na.rm=TRUE),Strict_rate=mean(Strict_rate,na.rm=TRUE),FDR_rate=mean(FDR_rate,na.rm=TRUE),Test_coverage=mean(Coverage_test,na.rm=TRUE),.groups='drop')
 overall<-per_fold%>%dplyr::group_by(Scheme)%>%dplyr::summarise(N_folds=dplyr::n(),Mean_direction=mean(Direction_rate),Mean_strict=mean(Strict_rate),Mean_FDR=mean(FDR_rate),Mean_test_coverage=mean(Test_coverage),.groups='drop')
 methods<-data.frame(Item=c('Goal','Primary splitting','Secondary splitting','Training','Evaluation','Missingness','No inner tuning','Leakage audit','Interpretation','No classification claim'),Value=c('Assess a discovery-conditioned selection procedure, not unbiased validation of the final fixed35-gene list.','Leave one source publication out among external cohorts. GSE240921,GSE240923 both batches,GSE242014 form one conservative source block; other external publications separate. Discovery always retained.','Accession-block omission keeps GSE240923 batches together. Cohort-wise omission is secondary and can retain a related batch/source in training; not primary independence evidence.','Frozen108 discovery Shared candidates; discovery RV support mandatory; require ceil(0.75*N_training) concordant significant votes including discovery. All training candidates evaluated; not restricted to final35.','Freeze selected genes and direction before viewing held-out effects. Report effect/test coverage, direction concordance, FDR-only concordance, and strict FDR<0.05 plus abs(log2FC)>=1 concordance. Up/down strata separate.','Unavailable tests excluded from evaluation denominators, not counted as negative; selected-gene and evaluable counts both reported. Training fixed-denominator support rule remains.','No data-adaptive hyperparameter search or inner tuning. This is outer held-out reselection with fixed rules, not a trained patient-level predictive model or full nested model tuning.','Selected-gene calls unchanged when held-out outcomes are perturbed; training/test dataset IDs explicitly checked disjoint. Baseline full-data selection reproduces35 genes.','Across-source rates descriptive; genes not independent replicates, few heterogeneous study blocks. No binomial precision claims. Discovery eligibility and source count-model limitations remain.','No patient-level expression/count model fitted across cohorts; no AUC, clinical diagnostic performance or unbiased fixed-panel validation is claimed.'))
 methods<-rbind(methods,data.frame(Item='Normalized-input sensitivity',Value='Publication_no_normalized excludes GSE198618 from both training and testing before folds are formed; ratios and ceil75% training thresholds are recalculated. Three external source blocks remain.'))
 list(folds=folds,training_calls=calls,heldout_genes=ev,metrics=mm,fold_summary=per_fold,scheme_summary=overall,unit_definitions=units,methods=methods)
}
