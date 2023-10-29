spfda_wrapper = function(
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
  lambdas = 10^seq(log10(1.), log10(100), length.out = 100)
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
    BIC(res)
  })

  BICfd = cbind(all_params, BIC=unlist(BICs))
  parms = BICfd[which.min(BICfd$BIC),  ]
  # ggplot() +
  #   geom_line(
  #     data=BICfd,
  #     mapping=aes(x=lambda, color=as.factor(alpha), group=paste0(K, alpha), y=BIC, linetype=as.factor(K))
  #     ) +
  #   scale_x_log10()

  out = spfda::spfda(
    Y=Y,
    X=X,
    time=instance$estimated_time,
    lambda=parms$lambda,
    alpha=parms$alpha,
    nsp=parms$K,
    ord=3,
    CI=F,
    W=W
  )

  est = out$get_coef()

  estimate = data.frame(
    time=instance$estimated_time,
    intercept=est[1, ],
    group_difference=est[2, ]
  )

  estimation_error = data.frame(
    time=instance$estimated_time,
    intercept=est[1, ] - instance$true_values$b0,
    group_difference=est[2, ] - instance$true_values$b1
  )

  decision = data.frame(
    time=instance$estimated_time,
    intercept=abs(est[1, ]) > 0,
    group_difference=abs(est[2, ]) > 0
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
    results= BICfd[which.min(BICfd$BIC),  ],
    fit=BICfd,
    time=proc.time()-t0
  )

}
