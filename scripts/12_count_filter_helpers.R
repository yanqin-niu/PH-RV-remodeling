# Common prefit rule. Applied to cohort expression matrices before contrast subsets.
prefilter_counts <- function(x,label='count matrix') {
 x<-as.matrix(x)
 if(any(!is.finite(x))||any(x<0))stop('Invalid count-like input: ',label)
 before<-nrow(x);x<-x[rowSums(x)>=10,,drop=FALSE]
 if(!nrow(x))stop('No features after rowSums>=10: ',label)
 if(any(rowSums(x)<10))stop('Prefilter failed: ',label)
 message('PREFILTER ',label,': ',before,' -> ',nrow(x),' features; ',ncol(x),' cohort columns; rowSums>=10')
 x
}
