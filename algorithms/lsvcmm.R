lsvcmm_wrapper = function(
  data, job, instance,
  selection="ebich",
  cross_sectional=F,
  independent=F,
  kernel.name="gaussian",
  penalty.alpha=1.,
  penalty.adaptive=1.,
  penalty.lambda=NULL,
  kernel.scale=NULL
){
  t0 = proc.time()
  df = instance$data
  if(cross_sectional){
    k_args = list(name="epa", scale=min(diff(sort(instance$estimated_time)))/2, n_scale=1L)
  }else{
    k_args = list(name=kernel.name, scale=kernel.scale, n_scale=1L)
  }
  if(independent){
    wc_args = list(name="independent")
  }else{
    wc_args = list(name="compound_symmetry", estimate=T, ratio=1.)
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
    kernel=k_args,
    working_covariance=wc_args,
    penalty=list(adaptive=penalty.adaptive, alpha=penalty.alpha, penalize_intercept=T,
                 lambda=penalty.lambda, nlambda=ifelse(is.null(penalty.lambda), 100L, 1L)),
    return_models=F
  )

  i = which.min(fit$results[[selection]])
  B = fit$vc_path[,,i]
  res = fit$results[i, ]

  estimate = data.frame(
    time=instance$estimated_time,
    intercept=B[1, ],
    group_difference=B[2, ]
  )

  estimation_error = data.frame(
    time=instance$estimated_time,
    intercept=B[1, ] - instance$true_values$b0,
    group_difference=B[2, ] - instance$true_values$b1
  )

  decision = data.frame(
    time=instance$estimated_time,
    intercept=abs(B[1, ]) > 0,
    group_difference=abs(B[2, ]) > 0
  )

  classification_error = data.frame(
    time=instance$estimated_time,
    intercept=ifelse(
      instance$true_values$b0 == 0,
      ifelse(decision$intercept, "FP", "TN"), # Negative: Detection, No Detection
      ifelse(decision$intercept, "TP", "FN") # Positive: Detection, No Detection
    ),
    group_difference=ifelse(
      instance$true_values$b1 == 0,
      ifelse(decision$group_difference, "FP", "TN"), # Negative: Detection, No Detection
      ifelse(decision$group_difference, "TP", "FN") # Positive: Detection, No Detection
    )
  )

  list(
    estimate=estimate,
    estimation_error=estimation_error,
    decision=decision,
    classification_error=classification_error,
    results=res,
    fit=fit,
    time=proc.time()-t0
  )
}
