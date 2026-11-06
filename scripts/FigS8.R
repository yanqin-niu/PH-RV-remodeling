# Separate-batch R3.3 sensitivity; GSE240921 remains ONE consensus contribution.
suppressPackageStartupMessages({library(DESeq2);library(AnnotationDbi);library(org.Hs.eg.db);library(openxlsx);library(ggplot2);library(patchwork)})
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
project<-normalizePath(file.path(sd,'..'))
out<-file.path(project,'outputs/FigS8_outputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(project,'scripts/12_count_filter_helpers.R'));source(file.path(project,'scripts/11_orthology_helpers.R'))
meta<-read.csv(file.path(project,'data/metadata/GSE240921_sample_metadata.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE)
x<-read.xlsx(file.path(project,'data/raw/GSE240921_processed-data-human.xlsx'),'count matrix',check.names=FALSE)
mat<-as.matrix(x[,-1]);storage.mode(mat)<-'numeric';rownames(mat)<-x[[1]]
stopifnot(identical(colnames(mat),meta$Sample),all(is.finite(mat)),all(mat>=0),max(abs(mat-round(mat)))<1e-8)
ids<-sub('\\..*$','',rownames(mat));sy<-mapIds(org.Hs.eg.db,keys=ids,keytype='ENSEMBL',column='SYMBOL',multiVals='first');sy<-unname(sy[ids])
keep<-!is.na(sy)&nzchar(sy);mat<-mat[keep,,drop=FALSE];rownames(mat)<-toupper(trimws(sy[keep]));if(anyDuplicated(rownames(mat)))mat<-rowsum(mat,group=rownames(mat));mat<-round(mat)
meta$Group<-relevel(factor(meta$Group),ref='Ctrl');meta$Sex<-factor(meta$Sex,levels=c('F','M'));meta$Batch<-factor(meta$Batch,levels=c('PE25','SE15'));rownames(meta)<-meta$Sample
primary<-read.xlsx(file.path(project,'outputs/Fig5_outputs/Figure5_source_data.xlsx'),'Consensus_genes')
panel<-data.frame(Gene=primary$Gene,Expected=ifelse(primary$anchor_dir=='Up',1,-1));panel$Human_query<-rat_query(panel$Gene,project)$Query
stopifnot(nrow(panel)==35,!anyNA(panel$Human_query))
joint<-read.csv(file.path(project,'outputs/Fig5_outputs/GSE240921_batch_sex_DESeq2.csv'));names(joint)[names(joint)=='log2FoldChange']<-'log2FC';joint$Model<-'Joint40_Batch_Sex';joint$Significant<-is.finite(joint$log2FC)&is.finite(joint$padj)&joint$padj<.05&abs(joint$log2FC)>=1
testable<-function(z)is.finite(z$log2FC)&is.finite(z$padj)
strict<-function(z)testable(z)&z$padj<.05&abs(z$log2FC)>=1
as_result<-function(dds,contrast,model){z<-as.data.frame(results(dds,contrast=contrast));names(z)[names(z)=='log2FoldChange']<-'log2FC';z$Gene<-rownames(z);z$Model<-model;z$Significant<-strict(z);z}
panel_result<-function(z,model){p<-merge(panel,z,by.x='Human_query',by.y='Gene',all.x=TRUE,sort=FALSE);p$Model<-model;p$Covered<-is.finite(p$log2FC);p$Testable<-testable(p);p$Direction_match<-p$Covered&sign(p$log2FC)==p$Expected;p$Strict_concordant<-strict(p)&p$Direction_match;p$FDR_concordant<-p$Testable&p$padj<.05&p$Direction_match;p$Strict_discordant<-strict(p)&!p$Direction_match;p}
res<-list();panels<-list();summary_rows<-list();design_rows<-list();pca_rows<-list();pca_variance<-list();factors<-list()
for(batch in c('SE15','PE25')) {
 i<-meta$Batch==batch;md<-droplevels(meta[i,,drop=FALSE]);m<-prefilter_counts(mat[,i,drop=FALSE],paste('GSE240921',batch))
 for(adjust in c(FALSE,TRUE)) {
   nm<-paste0(batch,if(adjust)'_Sex_Group' else '_Group');formula<-if(adjust)~Sex+Group else ~Group
   design<-model.matrix(formula,md);stopifnot(qr(design)$rank==ncol(design))
   design_rows[[nm]]<-data.frame(Model=nm,Formula=paste(deparse(formula),collapse=''),Samples=nrow(md),Count_filtered_genes=nrow(m),Parameters=ncol(design),Rank=qr(design)$rank,Condition_number=kappa(design))
   message('Fitting ',nm)
   dds<-DESeq(DESeqDataSetFromMatrix(m,md,formula),quiet=TRUE)
   z<-as_result(dds,c('Group','dRV','Ctrl'),nm);res[[nm]]<-z;panels[[nm]]<-panel_result(z,nm)
   summary_rows[[nm]]<-data.frame(Model=nm,Genes=nrow(z),Testable=sum(testable(z)),FDR05=sum(testable(z)&z$padj<.05),Strict=sum(strict(z)),Strict_up=sum(strict(z)&z$log2FC>0),Strict_down=sum(strict(z)&z$log2FC<0))
   write.csv(z,file.path(out,paste0(nm,'_DESeq2.csv')),row.names=FALSE,na='')
   factors[[nm]]<-data.frame(Model=nm,Sample=colnames(dds),Size_factor=sizeFactors(dds))
   if(!adjust){vs<-assay(varianceStabilizingTransformation(dds,blind=TRUE));top<-head(order(apply(vs,1,var),decreasing=TRUE),2000);pc<-prcomp(t(vs[top,,drop=FALSE]),center=TRUE,scale.=FALSE);variance<-pc$sdev^2/sum(pc$sdev^2);pca_rows[[batch]]<-data.frame(md,PC1=pc$x[,1],PC2=pc$x[,2],PCA_panel=sprintf('%s (PC1 %.1f%%, PC2 %.1f%%)',batch,100*variance[1],100*variance[2]));pca_variance[[batch]]<-data.frame(Batch=batch,PC=seq_along(variance),Fraction=variance)}
   rm(dds);gc()
 }
}
res$Joint40_Batch_Sex<-joint;panels$Joint40_Batch_Sex<-panel_result(joint,'Joint40_Batch_Sex')
summary_rows$Joint40_Batch_Sex<-data.frame(Model='Joint40_Batch_Sex',Genes=nrow(joint),Testable=sum(testable(joint)),FDR05=sum(testable(joint)&joint$padj<.05),Strict=sum(strict(joint)),Strict_up=sum(strict(joint)&joint$log2FC>0),Strict_down=sum(strict(joint)&joint$log2FC<0))
summary_rows<-do.call(rbind,summary_rows)
panel_summary<-do.call(rbind,lapply(names(panels),function(nm){p<-panels[[nm]];data.frame(Model=nm,Fixed_panel=35,Covered=sum(p$Covered),Testable=sum(p$Testable),Direction_match=sum(p$Direction_match),FDR_concordant=sum(p$FDR_concordant),Strict_concordant=sum(p$Strict_concordant),Strict_discordant=sum(p$Strict_discordant))}))
comparison_pairs<-list(c('SE15_Group','PE25_Group'),c('SE15_Sex_Group','PE25_Sex_Group'),c('SE15_Sex_Group','Joint40_Batch_Sex'),c('PE25_Sex_Group','Joint40_Batch_Sex'))
comparisons<-list();cross<-list()
for(pair in comparison_pairs) {
 nm<-paste(pair,collapse='_vs_');aa<-res[[pair[1]]];bb<-res[[pair[2]]];z<-merge(aa,bb,by='Gene',suffixes=c('_a','_b'));ii<-is.finite(z$log2FC_a)&is.finite(z$log2FC_b);ev<-ii&is.finite(z$padj_a)&is.finite(z$padj_b)
 comparisons[[nm]]<-data.frame(A=pair[1],B=pair[2],Common_effects=sum(ii),Spearman=cor(z$log2FC_a[ii],z$log2FC_b[ii],method='spearman'),Direction_agreement=sum(sign(z$log2FC_a[ii])==sign(z$log2FC_b[ii]))/sum(ii),Both_strict=sum(z$Significant_a&z$Significant_b),Both_strict_same_direction=sum(z$Significant_a&z$Significant_b&sign(z$log2FC_a)==sign(z$log2FC_b)))
 pa<-panels[[pair[1]]];pb<-panels[[pair[2]]];p<-merge(pa,pb,by='Gene',suffixes=c('_a','_b'));ee<-p$Covered_a&p$Covered_b
 comparisons[[nm]]$Panel_common_effects<-sum(ee)
 comparisons[[nm]]$Panel_direction_agree<-sum(ee&sign(p$log2FC_a)==sign(p$log2FC_b),na.rm=TRUE)
 comparisons[[nm]]$Panel_both_expected_direction<-sum(p$Direction_match_a&p$Direction_match_b)
 comparisons[[nm]]$Panel_both_FDR_concordant<-sum(p$FDR_concordant_a&p$FDR_concordant_b)
 comparisons[[nm]]$Panel_both_strict_concordant<-sum(p$Strict_concordant_a&p$Strict_concordant_b)
 if(pair[1]=='SE15_Sex_Group'&&pair[2]=='PE25_Sex_Group')cross[[nm]]<-z
}
comparisons<-do.call(rbind,comparisons)
# Complementary disease-specific interaction: difference between batch effects.
# All40 raw counts, group Ctrl reference; test (dRV-Ctrl)SE - (dRV-Ctrl)PE.
all_counts<-prefilter_counts(mat,'GSE240921 all40 interaction')
formula_int<-~Batch+Sex+Group+Batch:Group;mm<-model.matrix(formula_int,meta);stopifnot(qr(mm)$rank==ncol(mm))
message('Fitting all40 batch-by-group interaction')
dds<-DESeq(DESeqDataSetFromMatrix(all_counts,meta,formula_int),quiet=TRUE)
terms<-resultsNames(dds);term<-grep('BatchSE15[.]GroupdRV',terms,value=TRUE);stopifnot(length(term)==1)
interaction<-as.data.frame(results(dds,name=term));interaction$Gene<-rownames(interaction);interaction$Contrast<-'(dRV-Ctrl)SE15 minus (dRV-Ctrl)PE25, sex adjusted'
int_panel<-merge(panel,interaction,by.x='Human_query',by.y='Gene',all.x=TRUE)
int_summary<-data.frame(Contrast=interaction$Contrast[1],Genes=nrow(interaction),Testable=sum(is.finite(interaction$padj)),FDR05=sum(is.finite(interaction$padj)&interaction$padj<.05),Fixed_panel=35,Panel_covered=sum(is.finite(int_panel$log2FoldChange)),Panel_testable=sum(is.finite(int_panel$padj)),Panel_FDR05=sum(is.finite(int_panel$padj)&int_panel$padj<.05))
design_rows$Interaction<-data.frame(Model='Interaction',Formula=paste(deparse(formula_int),collapse=''),Samples=nrow(meta),Count_filtered_genes=nrow(all_counts),Parameters=ncol(mm),Rank=qr(mm)$rank,Condition_number=kappa(mm))
write.csv(interaction,file.path(out,'Batch_by_Group_interaction.csv'),row.names=FALSE,na='')
rm(dds);gc()
panel_data<-do.call(rbind,panels);pc<-do.call(rbind,pca_rows)
model_levels<-names(res);labels<-c(SE15_Group='SE15 group',SE15_Sex_Group='SE15 sex + group',PE25_Group='PE25 group',PE25_Sex_Group='PE25 sex + group',Joint40_Batch_Sex='Joint40 batch + sex + group')
summary_rows$Model<-factor(summary_rows$Model,levels=model_levels);panel_data$Model<-factor(panel_data$Model,levels=model_levels);panel_data$Gene<-factor(panel_data$Gene,levels=rev(primary$Gene));panel_data$Aligned_log2FC<-panel_data$Expected*panel_data$log2FC
pA<-ggplot(pc,aes(PC1,PC2,color=Group,shape=Sex))+geom_point(size=2.6)+facet_wrap(~PCA_panel,scales='free')+scale_color_manual(values=c(Ctrl='#4D4D4D',cRV='#E68613',dRV='#C83436'))+labs(title='A. Within-batch VST PCA',subtitle='Axes fitted separately; no sample exclusions')+theme_bw(base_size=9)
zz<-cross[[1]];zz$Panel<-zz$Gene%in%panel$Human_query
pB<-ggplot(zz,aes(log2FC_a,log2FC_b))+geom_point(size=.45,alpha=.3,color='grey65')+geom_hline(yintercept=0,linewidth=.25)+geom_vline(xintercept=0,linewidth=.25)+geom_abline(slope=1,intercept=0,linetype=2)+geom_point(data=subset(zz,Panel),color='#C83436',size=1.5)+labs(title='B. Sex-adjusted effects in separate batches',x='SE15 dRV versus Ctrl log2FC',y='PE25 dRV versus Ctrl log2FC')+theme_bw(base_size=9)
pC<-ggplot(summary_rows,aes(Model,Strict,fill=Model))+geom_col(width=.65)+geom_text(aes(label=Strict),vjust=-.3)+scale_fill_manual(values=c('#8CAAC5','#1F78B4','#8DC591','#33A02C','#C83436'))+scale_x_discrete(labels=labels)+labs(title='C. Significant disease effects',x=NULL,y='FDR < 0.05 and |log2FC| >= 1')+theme_bw(base_size=9)+theme(legend.position='none',axis.text.x=element_text(angle=25,hjust=1))
pD<-ggplot(panel_data,aes(Model,Gene,fill=Aligned_log2FC))+geom_tile(color='white',linewidth=.2)+geom_text(aes(label=ifelse(Strict_concordant,'*','')),size=2.6)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,name='Aligned\nlog2FC')+scale_x_discrete(labels=labels)+labs(title='D. Fixed original 35 genes',subtitle='* Significant concordant support; grey = unavailable',x=NULL,y=NULL)+theme_bw(base_size=8)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=25,hjust=1))
figure<-(pA|pB)/(pC|pD)+plot_layout(heights=c(1,1.8))
ggsave(file.path(out,'FigureS8.pdf'),figure,width=14,height=12)
ggsave(file.path(out,'FigureS8.png'),figure,width=14,height=12,dpi=220)
methods<-data.frame(Item=c('Separate fits','Filtering','Comparison','Interaction','Panel','Scope','Limitations'),Definition=c('All samples within each technical stratum retained; dRV-versus-Ctrl with cRV retained; group-only and sex+group models in each batch','Same annotation/symbol summation; rowSums>=10 separately within each batch, so fitted gene/test availability can differ','All-common finite effects and fixed35 are reported separately; strict support uses FDR<0.05, abs(log2FC)>=1 and expected direction','Exploratory Wald test of batch-by-dRV interaction with Ctrl and PE25 references, adjusting common sex effect; BH over transcriptome, not35; non-significance is not equivalence','Original35 genes and discovery signs fixed; reciprocal-unique curated human queries; missing tests unavailable','GSE240921 continues to contribute ONE joint adjusted result; separate batches never inserted as extra votes; prior main files unchanged','Small batches, imbalanced sex, different collection/library backgrounds and missing clinical covariates; differences in DEG counts do not prove different biology or poor quality'))
wb<-list(Samples=meta,Sample_sizes=as.data.frame(with(meta,table(Batch,Group,Sex))),Model_design=do.call(rbind,design_rows),Model_summary=summary_rows,Fixed35_summary=panel_summary,Fixed35_genes=panel_data,Comparisons=comparisons,Interaction_summary=int_summary,Interaction_panel=int_panel,Interaction_allgenes=interaction,PCA=pc,PCA_variance=do.call(rbind,pca_variance),Size_factors=do.call(rbind,factors),Methods=methods,Software=data.frame(Text=capture.output(sessionInfo())))
for(nm in names(res))wb[[paste0(nm,'_results')]]<-res[[nm]]
write.xlsx(wb,file.path(out,'FigureS8_source_data.xlsx'),overwrite=TRUE)
for(nm in c('Model_summary','Fixed35_summary','Fixed35_genes','Comparisons','Interaction_summary','Interaction_panel'))write.csv(wb[[nm]],file.path(out,paste0(nm,'.csv')),row.names=FALSE,na='')
writeLines(c(capture.output(summary_rows),capture.output(panel_summary),capture.output(comparisons),capture.output(int_summary)),file.path(out,'RESULTS.txt'))
message('Separate-batch analysis completed: ',out)

writeLines('Figure S8. Technical-stratum sensitivity of GSE240921. (A) VST PCA fitted separately by stratum; no sample exclusions. (B) Sex-adjusted dRV-versus-Ctrl log2FC estimates in SE15 and PE25. Red points mark the fixed35-panel human orthologues. (C) Strict DEG counts for group-only, sex-adjusted separate fits and the joint Batch+Sex+Group fit. (D) Discovery-direction-aligned effects for the fixed35, with stars denoting FDR<0.05, |log2FC|>=1 and expected direction; grey denotes unavailable effects. The primary reference is the same fixed panel in all columns. The interaction analysis in source data tests (dRV-Ctrl)SE15 minus (dRV-Ctrl)PE25 using all40 samples with sex adjustment and transcriptome-wide FDR. Nonsignificance is not equivalence. Technical strata may also differ in donor/collection backgrounds and are not treated as independent study-level votes.',file.path(out,'FigureS8_legend.txt'))
