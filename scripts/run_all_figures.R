# Canonical run order; isolated R processes prevent accidental object reuse.
a<-commandArgs(FALSE);hit<-grep('^--file=',a)
if(!length(hit))stop("Run with Rscript scripts/run_all_figures.R")
sd<-dirname(normalizePath(sub('^--file=','',a[hit[1]])))
root<-normalizePath(file.path(sd,'..'));logs<-file.path(root,'outputs','run_logs');dir.create(logs,recursive=TRUE,showWarnings=FALSE)
exe<-file.path(R.home('bin'),'Rscript.exe');if(!file.exists(exe))exe<-file.path(R.home('bin'),'Rscript')
figures<-c('Fig1','Fig2','Fig3','Fig5','Fig4','Fig6',paste0('FigS',1:8))
for(f in figures){message('Running ',f);log<-file.path(logs,paste0(f,'.log'));status<-system2(exe,shQuote(file.path(sd,paste0(f,'.R'))),stdout=log,stderr=log);if(status!=0)stop(f,' failed; inspect ',log)}
for(f in c(paste0('TableS',1:4),'software_versions')) {
 message('Running ',f);log<-file.path(logs,paste0(f,'.log'))
 status<-system2(exe,shQuote(file.path(sd,paste0(f,'.R'))),stdout=log,stderr=log)
 if(status!=0)stop(f,' failed; inspect ',log)
}
message('All 14 figures, four tables and software-version records completed')
