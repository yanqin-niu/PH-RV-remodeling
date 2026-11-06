# Robust sample-level scoring inference. No additional packages required.
score_fit <- function(formula,data) {
 f<-lm(formula,data=data);X<-model.matrix(f);stopifnot(qr(X)$rank==ncol(X),df.residual(f)>0)
 h<-hatvalues(f);stopifnot(all(h<1-1e-8));u<-residuals(f)/(1-h)
 bread<-solve(crossprod(X));V<-bread%*%crossprod(X,X*u^2)%*%bread
 dimnames(V)<-list(names(coef(f)),names(coef(f)));list(fit=f,V=V)
}
score_contrast <- function(f,c) {
 b<-coef(f$fit);stopifnot(length(c)==length(b));est<-sum(c*b);se<-sqrt(as.numeric(t(c)%*%f$V%*%c));df<-df.residual(f$fit);ci<-qt(.975,df)*se
 data.frame(Estimate=est,SE=se,Lower=est-ci,Upper=est+ci,DF=df,P=if(se>0)2*pt(-abs(est/se),df) else NA_real_)
}
score_term <- function(f,term){c<-as.numeric(names(coef(f$fit))==term);stopifnot(sum(c)==1);score_contrast(f,c)}
score_wald <- function(f,terms){i<-match(terms,names(coef(f$fit)));stopifnot(!anyNA(i));b<-coef(f$fit)[i];q<-length(i);stat<-as.numeric(t(b)%*%solve(f$V[i,i,drop=FALSE],b))/q;data.frame(F=stat,DF1=q,DF2=df.residual(f$fit),P=pf(stat,q,df.residual(f$fit),lower.tail=FALSE))}
score_sig <- function(p)ifelse(is.na(p),'NA',ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns'))))
score_group_inference <- function(x) {
 contrasts<-list();groups<-list()
 for(st in levels(x$Set)){
  z<-x[x$Set==st,];z$Group<-droplevels(z$Group);has_sex<-'Sex'%in%names(z)
  f<-score_fit(if(has_sex)Score~Group+Sex else Score~Group,z)
  terms<-grep('^Group',names(coef(f$fit)),value=TRUE)
  groups[[st]]<-cbind(Set=st,score_wald(f,terms),Model=if(has_sex)'Score~Group+Sex' else 'Score~Group')
  cc<-do.call(rbind,lapply(levels(z$Group)[-1],function(group){cbind(Set=st,group1=levels(z$Group)[1],group2=group,score_term(f,paste0('Group',group)))}))
  contrasts[[st]]<-cc
 }
 cc<-do.call(rbind,contrasts);cc$p.adj<-p.adjust(cc$P,method='holm');cc$p.signif<-score_sig(cc$p.adj)
 list(contrasts=cc,groups=do.call(rbind,groups))
}
score_statistics <- function(b1,b2,gene_audit) {
 trend<-list();paired<-list();sex<-list();sex_stage<-list();sex_interaction<-list();counts<-list();idx<-0L
 for(batch in c('Batch1','Batch2')){
  dat<-if(batch=='Batch1')b1 else b2;has_sex<-batch=='Batch2'
  dat$Stage<-as.integer(dat$Group)-1L;stopifnot(!anyNA(dat$Stage),all(table(dat$Sample)==3))
  counts[[batch]]<-as.data.frame(table(Group=dat$Group[dat$Set=='Shared'],Sex=if(has_sex)dat$Sex[dat$Set=='Shared'] else rep('M',sum(dat$Set=='Shared'))))
  for(scale in c('Raw','Standardized')){
   for(component in c('Score','Unsigned_score','Up_component','Down_component')){
    z<-dat;z$Value<-z[[component]]
    if(scale=='Standardized')for(st in levels(z$Set)){i<-z$Set==st;z$Value[i]<-as.numeric(base::scale(z$Value[i]))}
    for(st in levels(z$Set)){
     q<-z[z$Set==st & is.finite(z$Value),];if(!nrow(q))next
     f<-score_fit(if(has_sex)Value~Stage+Sex else Value~Stage,q)
     idx<-idx+1L;trend[[idx]]<-cbind(Batch=batch,Set=st,Component=component,Scale=scale,N=nrow(q),score_term(f,'Stage'))
    }
    if(component=='Score'){
     a<-z[z$Set=='Shared',];
     for(comp in c('LV-only','RV-only')){
      b<-z[z$Set==comp,];stopifnot(setequal(a$Sample,b$Sample));a$Delta<-a$Value-b$Value[match(a$Sample,b$Sample)]
      fit<-score_fit(if(has_sex)Delta~Stage+Sex else Delta~Stage,a)
      paired[[length(paired)+1L]]<-cbind(Batch=batch,Comparator=comp,Scale=scale,N=nrow(a),score_term(fit,'Stage'))
      # Difference-score slope must equal the difference of component slopes.
      fa<-score_fit(if(has_sex)Value~Stage+Sex else Value~Stage,a);fb<-score_fit(if(has_sex)Value~Stage+Sex else Value~Stage,b)
      stopifnot(abs(coef(fit$fit)['Stage']-(coef(fa$fit)['Stage']-coef(fb$fit)['Stage']))<1e-10)
     }
     if(has_sex)for(st in levels(z$Set)){
      q<-z[z$Set==st,];f<-score_fit(Value~Group+Sex,q)
      sex[[length(sex)+1L]]<-cbind(Set=st,Scale=scale,Contrast='M minus F, adjusted for categorical stage',N=nrow(q),score_term(f,'SexM'))
      fi<-score_fit(Value~Stage*Sex,q);sex_interaction[[length(sex_interaction)+1L]]<-cbind(Set=st,Scale=scale,Term='Ordinal stage slope difference M minus F',score_term(fi,'Stage:SexM'))
      fs<-score_fit(Value~Group*Sex,q)
      for(group in levels(q$Group)){
       nd<-data.frame(Group=factor(rep(group,2),levels=levels(q$Group)),Sex=factor(c('F','M'),levels=c('F','M')))
       mm<-model.matrix(delete.response(terms(fs$fit)),nd)
       sex_stage[[length(sex_stage)+1L]]<-cbind(Set=st,Scale=scale,Group=group,Contrast='M minus F within stage',score_contrast(fs,mm[2,]-mm[1,]))
      }
     }
    }
   }
  }
 }
 trend<-do.call(rbind,trend);trend$FDR<-NA_real_;for(sc in unique(trend$Scale))for(co in unique(trend$Component)){i<-trend$Scale==sc&trend$Component==co;trend$FDR[i]<-p.adjust(trend$P[i],'BH')}
 paired<-do.call(rbind,paired);paired$Holm_P<-NA_real_;for(sc in unique(paired$Scale)){i<-paired$Scale==sc;paired$Holm_P[i]<-p.adjust(paired$P[i],'holm')}
 sex<-do.call(rbind,sex);si<-do.call(rbind,sex_interaction);ss<-do.call(rbind,sex_stage)
 sex$Holm_P<-NA_real_;si$Holm_P<-NA_real_;ss$Holm_P<-NA_real_
 for(sc in unique(sex$Scale)){
  i<-sex$Scale==sc;sex$Holm_P[i]<-p.adjust(sex$P[i],'holm')
  i<-si$Scale==sc;si$Holm_P[i]<-p.adjust(si$P[i],'holm')
  i<-ss$Scale==sc;ss$Holm_P[i]<-p.adjust(ss$P[i],'holm')
 }
 coverage<-gene_audit%>%dplyr::group_by(Batch,Set)%>%dplyr::summarise(N_reference=dplyr::n(),N_present=sum(Measured),N_used=sum(Used),N_up=sum(Used&Direction>0),N_down=sum(Used&Direction<0),.groups='drop')
 meta<-data.frame(Item=c('Unit','Stage','Primary trend','Comparison','Standardization','Uncertainty','Sex','Multiplicity','Interpretation'),Value=c('Animal/sample is the observational unit; genes are not independent score replicates; batches analyzed separately.','Batch1 Ctrl,cRV,early_dRV,late_dRV encoded0,1,2,3;Batch2 Ctrl,cRV,dRV encoded0,1,2. Equal ordinal steps are not equal time intervals.','OLS score~Stage;Batch2 additionally adjusts Sex. Aligned/unsigned/up/down components exported.','Shared-minus-comparator score difference within the same sample regressed on Stage (plus Sex in batch2); preserves pairing across signatures.','Raw score primary stage associations; score divided by its within-batch,set SD for standardized slope comparisons; no across-batch pooling.','HC3 heteroskedasticity-robust covariance and residual-df t95%CI; ordinal trend assumes a linear association. Cross-sectional data do not establish longitudinal tracking.','Batch2 M-F adjusted for categorical Group, ordinal Stage*Sex slope interaction, and categorical-stage-specific M-F effects;Batch1 all male cannot estimate sex. Nonsignificance is not equivalence.','BH over6 batch-by-set trends within each component/scale;Holm over4 paired comparisons per scale;Holm over3 sex main/interaction and9 within-stage effects per scale;displayed control contrasts Holm within batch.','95% intervals are pointwise. Standardized slope difference is not a measure of predictive accuracy or clinical superiority.'))
 list(trends=trend,paired=paired,sex=sex,sex_by_stage=ss,sex_interaction=si,coverage=coverage,counts=counts,methods=meta)
}
