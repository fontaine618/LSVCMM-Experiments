lsvcmm = function(
    data, job, instance
){
  df = instance$data
  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    penalty=list(adaptive=1.)
  )
  # EBIC(0.5) Selection
  i = which.min(fit$results$ebich+fit$results$bich)
  return(fit$vc_path[2,,i])
}

lsvcmm_independent = function(
    data, job, instance
){
  df = instance$data
  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    penalty=list(adaptive=1.),
    working_covariance=list(estimate=F, ratio=0.)
  )
  # EBIC(0.5) Selection
  i = which.min(fit$results$ebich+fit$results$bich)
  return(fit$vc_path[2,,i])
}

lsvcmm_independent_cs = function(
    data, job, instance
){
  df = instance$data
  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    penalty=list(adaptive=1.),
    working_covariance=list(estimate=F, ratio=0.),
    kernel=list(scale=0.01)
  )
  # EBIC(0.5) Selection
  i = which.min(fit$results$ebich+fit$results$bich)
  return(fit$vc_path[2,,i])
}

cross_sectional = function(
    data, job, instance
){
  df = instance$data
  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    penalty=list(lambda=0.),
    working_covariance=list(estimate=F, ratio=0.),
    kernel=list(scale=0.01)
  )
  # EBIC(0.5) Selection
  i = which.min(fit$results$ebich+fit$results$bich)
  return(fit$vc_path[2,,i])
}


spfda = function(data, job, instance){
  require(dplyr)

  dat = instance$data_wide_imputed
  X = dat %>% select(group) %>% as.matrix()
  X = cbind(1, X)
  Y = dat %>% select(starts_with("t")) %>% as.matrix()
  t = instance$times

  W = spfda::spfda_weight(
    X=X,
    Y=Y,
    bandwidth=2,
    part=list(seq(length(t)))
  )

  alphas = seq(0.1, 0.9, by = 0.05)
  lambdas = 10^seq(log10(0.1), log10(10), length.out = 20)
  Ks = seq(5, 11)
  all_params = expand.grid(lambdas, alphas, Ks)
  names(all_params) = c("lambda", "alpha", "K")

  BICs = lapply(seq_len(nrow(all_params)), function(ii){
    param = all_params[ii,]
    lambda = param$lambda
    alpha = param$alpha
    K = param$K
    res = spfda::spfda(
      Y=Y,
      X=X,
      time=instance$times,
      lambda=lambda,
      alpha=alpha,
      nsp=K,
      ord=3,
      CI=F,
      W=W
    )
    BIC(res)
  })

  BICfd = cbind(all_params, BIC=unlist(BICs))
  parms = BICfd[which.min(BICfd$BIC),  ]

  out = spfda::spfda(
    Y=Y,
    X=X,
    time=instance$times,
    lambda=parms$lambda,
    alpha=parms$alpha,
    nsp=parms$K,
    ord=3,
    CI=F,
    W=W
  )

  return(out$get_coef()[2,])
}
