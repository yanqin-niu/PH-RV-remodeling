# Table S3: conditional LODO, threshold sensitivity and external-rat effect synthesis.
a<-commandArgs(FALSE);hit<-grep('^--file=',a);sd<-if(length(hit))dirname(sub('^--file=','',a[hit[1]])) else Sys.getenv('PH_RV_SCRIPT_DIR','')
if(!nzchar(sd))stop('Set PH_RV_SCRIPT_DIR when sourcing TableS3.R outside Rscript.')
source(file.path(sd,'00_common_helpers.R'));source(file.path(sd,'13_consensus_stability_helpers.R'))
root<-Sys.getenv('PH_RV_PROJECT_ROOT','');if(!nzchar(root))root<-find_project_dir(sd);out<-file.path(root,'outputs','TableS3_outputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
st<-consensus_stability(root);es<-external_effect_synthesis(root,st)
limits<-data.frame(Item=c('Effect model','Study weighting','Scope','Few studies','Selection'),Value=c('DerSimonian-Laird random effects; conservative modified Hartung-Knapp t95%CI with q bounded below by1; conventional Cochran-Q I2 and tau2. Positive estimate means discovery-consistent response.','One GSE240923 batch per fit; batch1 main, batch2 alternate. No pooled duplicate vote. Remaining animal-reuse uncertainty persists.','External rat hypoxia,SuHx,PAB,MCT only; discovery and humans excluded. Estimated-count hypoxia omitted in additional sensitivity. Fitted SEs used without requiring significance.','At least3 evaluable study units; heterogeneity and t intervals imprecise with3–4 studies. DL tau2 can be biased; no assertion of homogeneous biology.','All108 fixed discovery Shared candidates evaluated;35 primary genes flagged. Selection remains discovery-conditioned; BH adjustment is over estimable candidates within each scenario (89 of108 in this run), not only the35 primary genes. This does not validate the panel as a clinical assay.'))
write.xlsx(list(LODO_summary=st$summary,Gene_retention=st$retention,LODO_selected=st$lodo,Threshold_summary=st$thresholds,Threshold_gene_calls=st$threshold_calls,Availability=st$grid,Effect_synthesis=es$summary,Effect_inputs=es$inputs,Methods=st$methods,Effect_limits=limits,Session=data.frame(Text=capture.output(sessionInfo()))),file.path(out,'TableS3_consensus_stability.xlsx'),overwrite=TRUE)
write.csv(st$summary,file.path(out,'LODO_summary.csv'),row.names=FALSE);write.csv(st$thresholds,file.path(out,'Threshold_summary.csv'),row.names=FALSE);write.csv(es$summary,file.path(out,'Effect_synthesis.csv'),row.names=FALSE)
writeLines(c(
  'LODO rule and interpretation:',
  'The 108 discovery Shared genes were held fixed as the candidate universe. Each dataset was omitted once, retaining genes with at least 6/7 same-direction votes. Omitting GSE266139 removes its vote only; Shared membership is not recalculated.',
  'Of the final 35 genes, 21 had exactly 6/8 full-cohort support and were retained in exactly 2/8 omission runs by construction. The other 14 had at least 7/8 support and were retained in all eight runs. This threshold-conditioned retention is not independent evidence of biological stability.',
  '',
  'Per-omission Jaccard and retained-gene counts:',
  capture.output(st$summary[,1:7]),
  '',
  'Threshold sensitivity:',
  capture.output(st$thresholds[,1:6]),
  '',
  'Primary-gene LODO retention-frequency table:',
  capture.output(table(st$retention$n_LODO[st$retention$In_final]))
),file.path(out,'RESULTS.txt'))
message('Table S3 completed: ',out)
