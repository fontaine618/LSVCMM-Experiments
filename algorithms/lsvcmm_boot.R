instance = LSVCMM::generate_synthetic_data(n_timepoints=10)

lsvcmm_boot_wrapper = function(
    data, job, instance,
    selection="ebich",
    cross_sectional=F,
    independent=F,
    kernel.name="gaussian",
    kernel.scale=NULL,
    kernel.rescale_boundary=T,
    penalty.alpha=1.,
    penalty.adaptive=1.,
    penalty.lambda=NULL,
    penalty.name="adaptive_sparse_group_lasso",
    ar1.estimated=F,
    ar1.correlation=1.
){
  t0 = proc.time()
  df = instance$data
  if(cross_sectional){
    k_args = list(name="epa", scale=min(diff(sort(instance$estimated_time)))/2, n_scale=1L, rescale_boundary=kernel.rescale_boundary)
  }else{
    k_args = list(name=kernel.name, scale=kernel.scale, n_scale=1L, rescale_boundary=kernel.rescale_boundary)
  }
  if(independent){
    wc_args = list(name="independent")
  }else{
    wc_args = list(name="compound_symmetry", estimate=T, ratio=1.)
  }
  if(ar1.correlation < 1.){
    wc_args = list(name="autoregressive", estimate=T, correlation=ar1.correlation)
  }
  fit = LSVCMM::lsvcmm(
    data=df,
    response=instance$colnames$response,
    subject=instance$colnames$subject,
    time=instance$colnames$index,
    vc_covariates=instance$colnames$vc_covariates,
    nvc_covariates=instance$colnames$nvc_covariates,
    offset=instance$colnames$offset,
    add_intercept=T,
    estimated_time=instance$estimated_time,
    working_covariance=wc_args,
    kernel=k_args,
    penalty=list(name=penalty.name,
                 adaptive=penalty.adaptive, alpha=penalty.alpha, penalize_intercept=T,
                 lambda=penalty.lambda, nlambda=ifelse(is.null(penalty.lambda), 100L, 1L)
    ),
    return_models=F
  )

  i = which.min(fit$results[[selection]])
  h = fit$results$kernel.scale[i]
  l = fit$results$penalty.lambda[i]
  k_args$scale = h

  boot = LSVCMM::lsvcmm.boot(
    data=df,
    response=instance$colnames$response,
    subject=instance$colnames$subject,
    time=instance$colnames$index,
    vc_covariates=instance$colnames$vc_covariates,
    nvc_covariates=instance$colnames$nvc_covariates,
    offset=instance$colnames$offset,
    add_intercept=T,
    estimated_time=instance$estimated_time,
    working_covariance=wc_args,
    kernel=k_args,
    penalty=list(name=penalty.name,
                 adaptive=penalty.adaptive, alpha=penalty.alpha, penalize_intercept=T,
                 lambda=l, nlambda=ifelse(is.null(penalty.lambda), 100L, 1L)
    ),
    n_samples=1000
  )

  band = LSVCMM:::confidence_band(boot, var=2)

  pvals = data.frame(time=instance$estimated_time, pval=band$pval_percentile)
  pvals$pval.adj = p.adjust(pvals$pval, method="hommel")
  pvals$decision = band$excludes_zero

  classification = data.frame(
    time=instance$estimated_time,
    classification=ifelse(
      instance$true_values$b1 == 0,
      ifelse(pvals$decision, "FP", "TN"), # Negative: Detection, No Detection
      ifelse(pvals$decision, "TP", "FN") # Positive: Detection, No Detection
    )
  )

  list(
    decision=pvals,
    classification=classification,
    time=proc.time()-t0
  )



}
