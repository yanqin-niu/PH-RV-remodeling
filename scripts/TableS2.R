# Table S2: sample-identifier audit and study/model weighting sensitivity.
a<-commandArgs(FALSE);sd<-dirname(normalizePath(sub('^--file=','',a[grep('^--file=',a)][1])))
p<-Sys.getenv('PHRV_PYTHON',unset=Sys.which('python'))
if(!nzchar(p))p<-Sys.which('python3')
if(!nzchar(p))stop('Set PHRV_PYTHON to a Python executable with openpyxl')
status<-system2(p,shQuote(file.path(sd,'audit_orthology_study.py')))
if(status!=0)stop('Orthology audit failed')
status<-system2(p,shQuote(file.path(sd,'TableS2_weighting.py')))
if(status!=0)stop('R2.3 weighting failed')
