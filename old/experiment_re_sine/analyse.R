library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment_re_sine/registry/"

registry = loadRegistry(
  file.dir=DIR,
  writeable=T
)
# ------------------------------------------------------------------------------



# ==============================================================================
# Gather results
estimate = function(result) result$estimate
estimates = reduceResultsList(fun = estimate) %>% bind_rows(.id="job.id")
estimates %<>% mutate(job.id = as.numeric(job.id))

estimation_error = function(result) result$estimation_error
estimation_errors = reduceResultsList(fun = estimation_error) %>% bind_rows(.id="job.id")
estimation_errors %<>% mutate(job.id = as.numeric(job.id))

classification = function(result) result$classification_error
classifications = reduceResultsList(fun = classification) %>% bind_rows(.id="job.id")
classifications %<>% mutate(job.id = as.numeric(job.id))

parameters = getJobPars() %>% unwrap()
estimates %<>% left_join(parameters, by="job.id")
estimation_errors %<>% left_join(parameters, by="job.id")
classifications %<>% left_join(parameters, by="job.id")

write.csv(estimates, file="./experiment_re_sine/estimates.csv", row.names=F)
write.csv(estimation_errors, file="./experiment_re_sine/estimation_errors.csv", row.names=F)
write.csv(classifications, file="./experiment_re_sine/classifications.csv", row.names=F)
write.csv(parameters, file="./experiment_re_sine/parameters.csv", row.names=F)


setting = function(result){
  i = which.min(result$fit$results[["ebich"]])
  res = result$fit$results[i, ]
  return(res)
}
settings = reduceResultsList(fun = setting) %>% bind_rows(.id="job.id")
settings %<>% mutate(job.id = as.numeric(job.id))
settings %<>% left_join(parameters, by="job.id")
write.csv(settings, file="./experiment_re_sine/settings.csv", row.names=F)
# ------------------------------------------------------------------------------





# ==============================================================================
# multiple ICs
b1true = LSVCMM::generate_synthetic_data(n_timepoints=31, grpdiff_function="sine")$true_values$b1
prepare_results = function(result){
  fit = result$fit
  # get quantities
  llk = fit$results$llk
  n = 100*31*0.75 # wrong, but I can't easily get the number of observations
  # df_max = 22
  df = fit$results$df
  k = fit$results$kernel.scale
  dfs = fit$results$df_kernel
  dfs_lognhat = fit$results$df_logn_kernel
  df_lognhat = fit$results$df_logn
  df_logn = df * log(n)
  dfs_logn = dfs * log(n)
  df_max = fit$results$df_max
  B1 = fit$vc_path[2,,]
  t0 = fit$scaled_time
  # compute ICs
  AIC = -2*llk + 2*df
  AICs = -2*llk + 2*dfs
  BIC = -2*llk + df_logn
  BICs = -2*llk + dfs_logn
  BICnhat = -2*llk + df_lognhat
  BICsnhat = -2*llk + dfs_lognhat
  EBIC = -2*llk + df_logn + 1. * df * log(df_max)
  EBICs = -2*llk + dfs_logn + 1. * dfs * log(df_max)
  EBICnhat = -2*llk + df_lognhat + 1. * df * log(df_max)
  EBICsnhat = -2*llk + dfs_lognhat + 1. * dfs * log(df_max)
  # get index
  iAIC = which.min(AIC)
  iAICs = which.min(AICs)
  iBIC = which.min(BIC)
  iBICs = which.min(BICs)
  iBICnhat = which.min(BICnhat)
  iBICsnhat = which.min(BICsnhat)
  iEBIC = which.min(EBIC)
  iEBICs = which.min(EBICs)
  iEBICnhat = which.min(EBICnhat)
  iEBICsnhat = which.min(EBICsnhat)
  # get estimate
  bAIC = B1[, iAIC]
  bAICs = B1[, iAICs]
  bBIC = B1[, iBIC]
  bBICs = B1[, iBICs]
  bBICnhat = B1[, iBICnhat]
  bBICsnhat = B1[, iBICsnhat]
  bEBIC = B1[, iEBIC]
  bEBICs = B1[, iEBICs]
  bEBICnhat = B1[, iEBICnhat]
  bEBICsnhat = B1[, iEBICsnhat]
  # error
  eAIC = abs(bAIC - b1true)
  eAICs = abs(bAICs - b1true)
  eBIC = abs(bBIC - b1true)
  eBICs = abs(bBICs - b1true)
  eBICnhat = abs(bBICnhat - b1true)
  eBICsnhat = abs(bBICsnhat - b1true)
  eEBIC = abs(bEBIC - b1true)
  eEBICs = abs(bEBICs - b1true)
  eEBICnhat = abs(bEBICnhat - b1true)
  eEBICsnhat = abs(bEBICsnhat - b1true)
  # decision
  dAIC = ifelse(abs(bAIC)<1e-10, 0, 1)
  dAICs = ifelse(abs(bAICs)<1e-10, 0, 1)
  dBIC = ifelse(abs(bBIC)<1e-10, 0, 1)
  dBICs = ifelse(abs(bBICs)<1e-10, 0, 1)
  dBICnhat = ifelse(abs(bBICnhat)<1e-10, 0, 1)
  dBICsnhat = ifelse(abs(bBICsnhat)<1e-10, 0, 1)
  dEBIC = ifelse(abs(bEBIC)<1e-10, 0, 1)
  dEBICs = ifelse(abs(bEBICs)<1e-10, 0, 1)
  dEBICnhat = ifelse(abs(bEBICnhat)<1e-10, 0, 1)
  dEBICsnhat = ifelse(abs(bEBICsnhat)<1e-10, 0, 1)
  # store
  dfAIC = data.frame(
    ICname="AIC",
    ICtype="AIC",
    df="raw",
    n="raw",
    estimate=bAIC,
    error=eAIC,
    decision=dAIC,
    time=t0
  )
  dfAICs = data.frame(
    ICname="AICs",
    ICtype="AIC",
    df="smooth",
    n="raw",
    estimate=bAICs,
    error=eAICs,
    decision=dAICs,
    time=t0
  )
  dfBIC = data.frame(
    ICname="BIC",
    ICtype="BIC",
    df="raw",
    n="raw",
    estimate=bBIC,
    error=eBIC,
    decision=dBIC,
    time=t0
  )
  dfBICs = data.frame(
    ICname="BICs",
    ICtype="BIC",
    df="smooth",
    n="raw",
    estimate=bBICs,
    error=eBICs,
    decision=dBICs,
    time=t0
  )
  dfBICnhat = data.frame(
    ICname="BICnhat",
    ICtype="BIC",
    df="raw",
    n="nhat",
    estimate=bBICnhat,
    error=eBICnhat,
    decision=dBICnhat,
    time=t0
  )
  dfBICsnhat = data.frame(
    ICname="BICsnhat",
    ICtype="BIC",
    df="smooth",
    n="nhat",
    estimate=bBICsnhat,
    error=eBICsnhat,
    decision=dBICsnhat,
    time=t0
  )
  dfEBIC = data.frame(
    ICname="EBIC",
    ICtype="EBIC",
    df="raw",
    n="raw",
    estimate=bEBIC,
    error=eEBIC,
    decision=dEBIC,
    time=t0
  )
  dfEBICs = data.frame(
    ICname="EBICs",
    ICtype="EBIC",
    df="smooth",
    n="raw",
    estimate=bEBICs,
    error=eEBICs,
    decision=dEBICs,
    time=t0
  )
  dfEBICnhat = data.frame(
    ICname="EBICnhat",
    ICtype="EBIC",
    df="raw",
    n="nhat",
    estimate=bEBICnhat,
    error=eEBICnhat,
    decision=dEBICnhat,
    time=t0
  )
  dfEBICsnhat = data.frame(
    ICname="EBICsnhat",
    ICtype="EBIC",
    df="smooth",
    n="nhat",
    estimate=bEBICsnhat,
    error=eEBICsnhat,
    decision=dEBICsnhat,
    time=t0
  )
  # combine
  df = rbind(dfAIC, dfAICs, dfBIC, dfBICs, dfBICnhat, dfBICsnhat, dfEBIC, dfEBICs, dfEBICnhat, dfEBICsnhat)
  # add classification
  df %<>% mutate(
    classification=ifelse(
      abs(time-0.5)>0.25,
      ifelse(decision, "FP", "TN"), # Negative: Detection, No Detection
      ifelse(decision, "TP", "FN") # Positive: Detection, No Detection
    )
  )
  return(df)
}

# only get LSVCMM jobs
job.ids = findExperiments(algo.pattern="LSVCMM") %>% select(job.id) %>% unique()
icresults = reduceResultsList(ids = job.ids, fun = prepare_results) %>% bind_rows(.id="job.id")
icresults %<>% mutate(job.id = as.numeric(job.id))
write.csv(icresults, file="./experiment_re_sine/icresults.csv", row.names=F)
# ------------------------------------------------------------------------------
