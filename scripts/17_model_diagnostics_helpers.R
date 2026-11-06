# Transparent basic GLMM diagnostics; not a claim of definitive model adequacy.
audit_glmm <- function(fit,label) {
  stopifnot(inherits(fit,'merMod'))
  grad<-fit@optinfo$derivs$gradient;hess<-fit@optinfo$derivs$Hessian
  scaled<-tryCatch(solve(chol(hess),grad),error=function(e) NA_real_)
  eigenvalues<-tryCatch(eigen(hess,symmetric=TRUE,only.values=TRUE)$values,error=function(e) NA_real_)
  messages<-fit@optinfo$conv$lme4$messages
  pearson<-residuals(fit,type='pearson');y<-model.response(model.frame(fit));prob<-fitted(fit)
  summary<-data.frame(Model=label,Formula=paste(deparse(formula(fit)),collapse=' '),
    N=nobs(fit),Events=sum(y),Non_events=sum(1-y),
    N_genes=nlevels(fit@flist$Gene),N_datasets=nlevels(fit@flist$Dataset),
    Optimizer=paste(fit@optinfo$optimizer,collapse=';'),
    Optimizer_code=paste(fit@optinfo$conv$opt,collapse=';'),
    Convergence_message=if(length(messages))paste(messages,collapse=';') else 'None',
    Singular=lme4::isSingular(fit,tol=1e-4),
    Max_abs_gradient=if(length(grad))max(abs(grad)) else NA_real_,
    Max_abs_scaled_gradient=if(all(is.finite(scaled)))max(abs(scaled)) else NA_real_,
    Hessian_min_eigenvalue=min(eigenvalues),
    Pearson_X2=sum(pearson^2),Residual_df=df.residual(fit),
    Pearson_X2_per_df=sum(pearson^2)/df.residual(fit),
    Note='Conditional Pearson ratio and in-sample calibration are descriptive for Bernoulli GLMMs; they do not exclude misspecification or residual dependence.')
  variances<-as.data.frame(lme4::VarCorr(fit));variances$Model<-label
  # Rank deciles prevent quantile ties from producing invalid cut points.
  bins<-pmin(10L,ceiling(rank(prob,ties.method='first')/length(prob)*10))
  calibration<-do.call(rbind,lapply(sort(unique(bins)),function(b) {
    i<-bins==b;data.frame(Model=label,Bin=b,N=sum(i),Predicted=mean(prob[i]),Observed=mean(y[i]))
  }))
  list(Summary=summary,Random_effects=variances,Calibration=calibration)
}
