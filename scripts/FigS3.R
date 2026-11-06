# Supplementary analysis R2.1: fixed discovery-RV covariates, matched RV-only comparison,
# and continuous LV directional support. Does not alter Fig4 or panel selection.
suppressPackageStartupMessages({library(DESeq2);library(openxlsx);library(dplyr);library(ggplot2);library(patchwork);library(lme4)})
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
root<-normalizePath(file.path(sd,'..'));out<-file.path(root,'outputs','FigS3_outputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(sd,'12_count_filter_helpers.R'))
set.seed(20261102);B<-2000L
raw<-read.xlsx(file.path(root,'data/raw/GSE266139_12967_2025_6792_MOESM2_ESM.xlsx'))
cols<-grep('^(LP|LM|RP|RM)_\\d+$',names(raw),value=TRUE)
m<-as.matrix(raw[,cols]);storage.mode(m)<-'numeric';keep<-rowSums(is.na(m))==0
m<-round(m[keep,,drop=FALSE]);ann<-raw[keep,c('Gene.stable.ID','Gene.name')];rownames(m)<-ann$Gene.stable.ID
keep<-rowSums(m)>=10;m<-prefilter_counts(m[keep,,drop=FALSE],"GSE266139 all20")
fit_one<-function(a,b){s<-grepl(paste0('^(',a,'|',b,')_'),colnames(m));cd<-data.frame(condition=droplevels(factor(sub('_.*','',colnames(m)[s]),levels=c('LP','RP','LM','RM'))),row.names=colnames(m)[s]);dds<-DESeq(DESeqDataSetFromMatrix(m[,s],cd,~condition),quiet=TRUE);r<-as.data.frame(results(dds,contrast=c('condition',a,b)));r$Geneid<-rownames(r);r$dispersion<-dispersions(dds);r}
rv<-fit_one('RM','RP');lv<-fit_one('LM','LP')
old<-read.xlsx(file.path(root,'outputs/Fig1_outputs/Figure1_source_data.xlsx'),sheet='RMvsRP')
ix<-match(old$Geneid,rv$Geneid);stopifnot(all(abs(old$log2FoldChange-rv$log2FoldChange[ix])<1e-5,na.rm=TRUE))
v<-merge(rv,lv,by='Geneid',suffixes=c('_rv','_lv'));v$Gene<-toupper(ann$Gene.name[match(v$Geneid,ann$Gene.stable.ID)])
sets<-read.xlsx(file.path(root,'outputs/Fig1_outputs/Figure1_source_data.xlsx'),sheet='LV_RV_merge')
v$Set<-as.character(sets$set[match(v$Geneid,sets$Geneid)])
v<-v[!is.na(v$Gene)&nzchar(v$Gene)&v$Set%in%c('Shared','RV-only')&is.finite(v$lfcSE_rv)&v$lfcSE_rv>0&is.finite(v$lfcSE_lv)&v$lfcSE_lv>0,]
# Keep a deterministic representative stable ID for symbol-level external analyses.
v<-v[order(v$Gene,v$padj_rv,-abs(v$log2FoldChange_rv),v$Geneid),];dup<-v[duplicated(v$Gene)|duplicated(v$Gene,fromLast=TRUE),];v<-v[!duplicated(v$Gene),]
v$abs_RV<-abs(v$log2FoldChange_rv);v$log_SE<-log(v$lfcSE_rv);v$log_mean<-log2(v$baseMean_rv+1);v$RV_direction<-sign(v$log2FoldChange_rv)
v$LV_support<-v$RV_direction*v$log2FoldChange_lv/v$lfcSE_lv
# Frozen curated orthology reference, reciprocal-unique human queries only.
o<-read.delim(file.path(root,'data/reference/ORTHOLOGS_RAT.txt'),comment.char='#',check.names=FALSE,quote='');pairs<-unique(data.frame(rat=toupper(o[[1]]),human=toupper(o[[5]])));pairs<-pairs[nzchar(pairs$rat)&nzchar(pairs$human),];nr<-table(pairs$rat);nh<-table(pairs$human);good<-pairs[nr[pairs$rat]==1&nh[pairs$human]==1,]
x<-read.xlsx(file.path(root,'outputs/Fig5_outputs/Figure5_source_data.xlsx'),sheet='Native_results');ds<-setdiff(unique(x$Dataset),'GSE266139_RV (Rat_MCT)')
grid<-do.call(rbind,lapply(ds,function(d){h<-grepl('Human',d);q<-if(h)good$human[match(v$Gene,good$rat)] else v$Gene;z<-x[x$Dataset==d,];if(anyDuplicated(z$Gene))z<-z[!duplicated(z$Gene),];r<-z[match(q,z$Gene),];data.frame(Gene=v$Gene,Dataset=d,Set=v$Set,measured=is.finite(r$log2FC)&is.finite(r$padj),replicate=is.finite(r$log2FC)&is.finite(r$padj)&r$padj<.05&abs(r$log2FC)>=1&sign(r$log2FC)==v$RV_direction)}))
v$detectability<-vapply(v$Gene,function(g)mean(grid$measured[grid$Gene==g]),numeric(1))
vars<-c('abs_RV','log_SE','log_mean');z<-scale(v[,vars]);rownames(z)<-v$Gene
s<-which(v$Set=='Shared');c<-which(v$Set=='RV-only');matches<-list()
# Fixed within-variable caliper: 0.5 pooled SD, exact RV direction.
for(i in s){cand<-c[v$RV_direction[c]==v$RV_direction[i]];delta<-sweep(z[cand,,drop=FALSE],2,z[i,],'-');cand<-cand[apply(abs(delta)<=.5,1,all)];if(length(cand)){dist<-rowSums(sweep(z[cand,,drop=FALSE],2,z[i,],'-')^2);j<-cand[which.min(dist)];matches[[length(matches)+1L]]<-data.frame(Shared=v$Gene[i],Comparator=v$Gene[j],Distance=sqrt(min(dist)))}}
if(!length(matches))stop('No common-support matches under prespecified caliper')
matched<-do.call(rbind,matches);weights<-table(matched$Comparator)
v$weight<-ifelse(v$Set=='Shared',as.numeric(v$Gene%in%matched$Shared),as.numeric(weights[v$Gene]));v$weight[is.na(v$weight)]<-0
balance<-do.call(rbind,lapply(c(vars,'detectability'),function(k){a<-v[v$Set=='Shared',k];b<-v[v$Set=='RV-only',k];den<-sqrt((var(a)+var(b))/2);data.frame(Covariate=k,Before=(mean(a)-mean(b))/den,After=(weighted.mean(v[v$Set=='Shared',k],v$weight[v$Set=='Shared'])-weighted.mean(v[v$Set=='RV-only',k],v$weight[v$Set=='RV-only']))/den)}))
grid<-merge(grid,v[,c('Gene','weight','abs_RV','log_SE','log_mean','LV_support','detectability')],by='Gene')
# Poisson multiplier cluster bootstrap: one weight per unique gene, reused across
# all datasets and every reuse of a matched comparator; fixed matching design.
mult<-matrix(rpois(B*nrow(v),1),nrow=nrow(v));rownames(mult)<-v$Gene
rates<-function(d){g<-grid[grid$Dataset==d,];available<-g$Gene[g$measured];mp<-matched[matched$Shared%in%available & matched$Comparator%in%available,,drop=FALSE];ww<-table(c(mp$Shared,mp$Comparator));r<-g[g$Gene%in%names(ww),];r$weight<-as.numeric(ww[r$Gene]);a<-r$Set=='Shared';b<-!a;if(!any(a)||!any(b))return(NULL);wa<-r$weight[a];wb<-r$weight[b];da<-weighted.mean(r$replicate[a],wa);db<-weighted.mean(r$replicate[b],wb);ma<-mult[r$Gene[a],,drop=FALSE]*wa;mb<-mult[r$Gene[b],,drop=FALSE]*wb;draw<-colSums(ma*r$replicate[a])/colSums(ma)-colSums(mb*r$replicate[b])/colSums(mb);list(row=data.frame(Dataset=d,N_Shared=sum(a),N_complete_pairs=nrow(mp),N_distinct_comparator=sum(b),Shared_rate=da,Comparator_rate=db,Difference=da-db,Lower=quantile(draw,.025,na.rm=TRUE),Upper=quantile(draw,.975,na.rm=TRUE)),draw=draw)}
r<-lapply(ds,rates);r<-Filter(Negate(is.null),r);summary<-do.call(rbind,lapply(r,`[[`,'row'));draw<-rowMeans(do.call(cbind,lapply(r,`[[`,'draw')),na.rm=TRUE);pooled<-data.frame(Dataset='Equal-dataset mean',Difference=mean(summary$Difference),Lower=quantile(draw,.025,na.rm=TRUE),Upper=quantile(draw,.975,na.rm=TRUE))
md<-grid[grid$measured,];for(k in c(vars,'LV_support','detectability'))md[[paste0('z_',k)]]<-as.numeric(scale(md[[k]]))
# Continuous association among RV-selected genes, not a causal effect.
fit<-glmer(replicate~z_LV_support+z_abs_RV+z_log_SE+z_log_mean+(1|Dataset)+(1|Gene),data=md,family=binomial(),control=glmerControl(optimizer='bobyqa'))
cf<-fixef(fit);se<-sqrt(diag(vcov(fit)));or<-data.frame(Term=names(cf),OR=exp(cf),Lower=exp(cf-1.96*se),Upper=exp(cf+1.96*se))
fit2<-update(fit,.~.+z_detectability);cf2<-fixef(fit2);se2<-sqrt(diag(vcov(fit2)));or2<-data.frame(Term=names(cf2),OR=exp(cf2),Lower=exp(cf2-1.96*se2),Upper=exp(cf2+1.96*se2))
# Descriptive binned continuous association; model estimate is displayed separately.
md$Bin<-cut(md$LV_support,unique(quantile(md$LV_support,seq(0,1,.2))),include.lowest=TRUE)
continuous<-md%>%group_by(Bin)%>%summarise(LV_support=median(LV_support),Rate=mean(replicate),N=n(),.groups='drop')
# Manuscript palette: Shared red, RV-only green, effect-size orange,
# Human purple / Rat green. Colors encode diagnostics/species/summary,
# rather than suggesting separate estimates for a single contrast.
source(file.path(sd,"00_common_helpers.R"))
balance_plot<-rbind(data.frame(Covariate=balance$Covariate,SMD=abs(balance$Before),Stage="Before matching"),
                    data.frame(Covariate=balance$Covariate,SMD=abs(balance$After),Stage="After matching"))
balance_plot$Covariate<-factor(balance_plot$Covariate,levels=c('abs_RV','log_mean','log_SE','detectability'),
 labels=c('RV |log2FC|','RV expression abundance','RV standard error','External test availability'))
balance_plot$Stage<-factor(balance_plot$Stage,levels=c('Before matching','After matching'))
pA<-ggplot(balance_plot,aes(SMD,Covariate,color=Stage,shape=Stage))+
 geom_vline(xintercept=.1,linetype=2,color='grey55',linewidth=.4)+
 geom_point(size=2.8,stroke=.8)+
 scale_color_manual(values=c('Before matching'='#E68613','After matching'='#1F78B4'))+
 scale_shape_manual(values=c('Before matching'=1,'After matching'=16))+
 labs(title='A. Covariate balance',x='Absolute standardized difference',y=NULL,
 subtitle='Test availability audited, not matched',color=NULL,shape=NULL)+
 theme_manuscript(10)+theme(legend.position='top',plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
plot_summary<-summary;plot_summary$Species<-ifelse(grepl('Human',plot_summary$Dataset),'Human','Rat')
pB<-ggplot(plot_summary,aes(Difference,reorder(Dataset,Difference),color=Species))+
 geom_vline(xintercept=0,linetype=2,color='grey55',linewidth=.4)+
 geom_segment(aes(x=Lower,xend=Upper,yend=reorder(Dataset,Difference)),linewidth=.85,lineend='round')+
 geom_point(size=2.8)+scale_color_manual(values=species_cols)+
 labs(title='B. Matched replication difference',x='Shared minus RV-only',y=NULL,color=NULL)+
 theme_manuscript(10)+theme(legend.position='top',plot.title=element_text(hjust=0))
pC<-ggplot(pooled,aes(Difference,Dataset))+
 geom_vline(xintercept=0,linetype=2,color='grey55',linewidth=.4)+
 geom_segment(aes(x=Lower,xend=Upper,yend=Dataset),color=set_cols['Shared'],linewidth=1.15,lineend='round')+
 geom_point(shape=21,size=4.2,fill=set_cols['Shared'],color='white',stroke=.8)+
 labs(title='C. Equal-dataset average',x='Difference (gene-cluster bootstrap 95% CI)',y=NULL)+
 theme_manuscript(10)+theme(plot.title=element_text(hjust=0))
pD<-ggplot(continuous,aes(LV_support,Rate))+
 geom_line(color='#E68613',linewidth=.85)+
 geom_point(aes(size=N),shape=21,fill='#E68613',color='white',stroke=.65)+
 scale_size_continuous(range=c(2.8,5),name='Gene–dataset\nobservations')+
 labs(title='D. Continuous LV directional support',subtitle='Descriptive bins among RV-selected genes',
 x='sign(RV effect) × LV Wald statistic',y='Observed replication rate')+
 theme_manuscript(10)+theme(plot.title=element_text(hjust=0),plot.subtitle=element_text(size=9))
figS3<-(pA|pB)/(pC|pD)
ggsave(file.path(out,'FigureS3.pdf'),figS3,width=13,height=8)
ggsave(file.path(out,'FigureS3.png'),figS3,width=13,height=8,dpi=300)
for(panel in c('A','B','C','D')) ggsave(file.path(out,paste0('FigS3',panel,'.pdf')),get(paste0('p',panel)),width=6.5,height=4)
meta<-data.frame(Item=c('Matching','Bootstrap','Scope','Continuous model','Model diagnostics','Caution'),Value=c('1:1 nearest neighbor with replacement; exact RV direction; 0.5 pooled SD caliper for all3 covariates; no tuning to external outcomes',paste(B,'Poisson gene-cluster multiplier draws; fixed matches; equal-dataset pooled average'),'Shared vs RV-only; common-support subset, not whole Shared set; curated human orthology; finite external estimate and padj required','sign(RV log2FC)*LV log2FC/LV lfcSE; adjusted for RV effect, SE, baseMean; detectability added as sensitivity',paste('Singular primary:',isSingular(fit),'Singular sensitivity:',isSingular(fit2)), 'Complete matched pairs required in each external dataset. Association only; matching uncertainty and generalization beyond observed datasets not quantified'))
source(file.path(sd,'17_model_diagnostics_helpers.R'))
diag_primary<-audit_glmm(fit,'Continuous LV support');diag_availability<-audit_glmm(fit2,'Availability sensitivity')
write.xlsx(list(Model_diagnostics=rbind(diag_primary$Summary,diag_availability$Summary),Random_effects=rbind(diag_primary$Random_effects,diag_availability$Random_effects),Calibration=rbind(diag_primary$Calibration,diag_availability$Calibration),Discovery_RV=rv,Discovery_LV=lv,Candidate_genes=v,Matched_sample_sizes=data.frame(N_Shared_available=sum(v$Set=='Shared'),N_RV_only_available=sum(v$Set=='RV-only'),N_Shared_matched=nrow(matched),N_distinct_comparators=length(weights),Comparator_ESS=sum(weights)^2/sum(weights^2)),Duplicate_symbols=dup,Matched_pairs=matched,Balance=balance,Dataset_results=summary,Pooled=pooled,Continuous_OR=or,Detectability_OR=or2,Continuous_bins=continuous,Replication_grid=grid,Meta=meta,sessionInfo=data.frame(Text=capture.output(sessionInfo()))),file.path(out,'FigureS3_source_data.xlsx'),overwrite=TRUE)
writeLines(c(capture.output(summary),capture.output(pooled),capture.output(balance),capture.output(or),capture.output(meta)),file.path(out,'RESULTS.txt'))
message('Figure S3 finished: ',out)
