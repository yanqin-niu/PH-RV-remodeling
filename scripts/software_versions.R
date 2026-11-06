# Export the actual environment; this does not install or update packages.
a<-commandArgs(FALSE);hit<-grep('^--file=',a)
sd<-dirname(normalizePath(sub('^--file=','',a[hit[1]])))
root<-normalizePath(file.path(sd,'..'));out<-file.path(root,'outputs','software_versions')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-list.files(sd,pattern='\\.R$',full.names=TRUE)
text<-paste(unlist(lapply(files,readLines,warn=FALSE)),collapse='\n')
ns<-unique(regmatches(text,gregexpr('[A-Za-z][A-Za-z0-9.]*::',text))[[1]])
libs<-unique(regmatches(text,gregexpr('library\\([A-Za-z][A-Za-z0-9.]*',text))[[1]])
packages<-sort(unique(c(sub('::$','',ns),sub('library\\(','',libs))))
versions<-data.frame(Package=packages,Version=vapply(packages,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else 'NOT INSTALLED',character(1)))
write.csv(versions,file.path(out,'R_package_versions.csv'),row.names=FALSE)
writeLines(c(R.version.string,capture.output(sessionInfo())),file.path(out,'R_sessionInfo.txt'))
python<-Sys.getenv('PHRV_PYTHON',unset=Sys.which('python'))
if(nzchar(python)) {
 py<-tryCatch(system2(python,c('-c',shQuote('import sys,openpyxl;print(sys.version);print("openpyxl "+openpyxl.__version__)')),stdout=TRUE,stderr=TRUE),error=function(e)conditionMessage(e))
 writeLines(py,file.path(out,'Python_versions.txt'))
}
message('Software environment exported to ',out)
