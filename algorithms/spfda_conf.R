spfda_conf_wrapper = function(
    data, job, instance,
    independent=F,
    oracle=F, K=NULL
){

  t0 = proc.time()
  dat = instance$data_wide_imputed
  X = dat %>% select(instance$colnames_wide$vc_covariates) %>% as.matrix()
  X = cbind(1, X)
  Y = dat %>% select(starts_with("t")) %>% as.matrix()
  t = instance$estimated_time
  nt = length(t)

  if(!independent){
    W = spfda::spfda_weight(
      X=X,
      Y=Y,
      bandwidth=2,
      part=list(seq(length(t)))
    )
  }else{
    W = NULL
  }

  alphas = c(0.5)
  lambdas = 10^seq(log10(10), log10(200), length.out = 100)
  if(is.null(K)) Ks = c(floor(nt/2)) else Ks = c(K)
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
      time=instance$estimated_time,
      lambda=lambda,
      alpha=alpha,
      nsp=K,
      ord=3,
      CI=F,
      W=W
    )
    df = sum(abs(res$get_coef()[2, ])>1e-10)
    c(df, BIC(res))
  })

  BICfd = cbind(
    all_params,
    BIC=unlist(BICs)[seq(2, length(unlist(BICs)), 2)],
    df=unlist(BICs)[seq(1, length(unlist(BICs)), 2)]
  )
  parms = BICfd[which.min(BICfd$BIC),  ]

  out = spfda::spfda(
    Y=Y,
    X=X,
    time=instance$estimated_time,
    lambda=parms$lambda,
    alpha=parms$alpha,
    nsp=parms$K,
    ord=3,
    CI=T,
    W=W
  )

  est = out$get_coef()
  se = out$get_se()

  band = data.frame(time=t, estimate=est[2,], se=se[2,])
  band$Lower = band$estimate - 2*band$se
  band$Upper = band$estimate + 2*band$se
  band$excludes_zero = band$Lower > 0 | band$Upper < 0

  pvals = data.frame(time=instance$estimated_time, decision=band$excludes_zero)

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
    fit=out,
    time=proc.time()-t0
  )

}
