library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)


# ==============================================================================
# Setup batchtools registry
name = "experiment_pvalues"
DIR = paste0("./", name, "/")
DIR_REGISTRY = paste0("./", name, "/registry/")
if(dir.exists(DIR_REGISTRY)) unlink(DIR, recursive=T)
if(!dir.exists(DIR)) dir.create(DIR, recursive=T)
registry = makeExperimentRegistry(
  file.dir=DIR_REGISTRY,
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM")
)
registry$cluster.functions = makeClusterFunctionsMulticore(ncpus=15)
# ------------------------------------------------------------------------------



# ==============================================================================
# Setup problem
synthetic = function(
  data, job,
  scenario=c("A", "B", "C", "D"),
  n_subjects=100,
  n_timepoints=11,
  prop_observed=0.7,
  observation_variance=1.,
  random_effect_variance_ratio=1.,
  random_effect_ar1_correlation=0.9,
  effect_size=1,
  missingness="uniform",
  seed=1
){
  set.seed(seed)
  # generate errors
  corr = random_effect_ar1_correlation^(seq(0, n_timepoints-1) / (n_timepoints-1))
  corrmat = stats::toeplitz(corr)
  thetamat = mvtnorm::rmvnorm(n_subjects, sigma=corrmat) * sqrt(random_effect_variance_ratio) * sqrt(observation_variance)
  t0 = seq(0, n_timepoints-1) / (n_timepoints-1)
  timemat = matrix(t0, n_subjects, n_timepoints, byrow=T)
  # generate covariates and mean
  if(scenario=="A"){
    # single VC with binary variable
    # short group difference
    Xmat = matrix(rbinom(n_subjects, 1, 0.5), n_subjects, 1)
    f1raw = function(t) 1/(1+exp((0.8-t)*20))
    f1 = function(t) ifelse(abs(f1raw(t)) < 0.05, 0, f1raw(t)) # make small values exact 0s
    Bmat = matrix(f1(t0), 1, n_timepoints)
    mu = Xmat %*% Bmat
  }
  if(scenario=="B"){
    # single VC with binary variable
    # long group difference
    Xmat = matrix(rbinom(n_subjects, 1, 0.5), n_subjects, 1)
    f1raw = function(t) 1-1/(1+exp((0.8-t)*20))
    f1 = function(t) ifelse(abs(f1raw(t)) < 0.05, 0, f1raw(t)) # make small values exact 0s
    Bmat = matrix(f1(t0), 1, n_timepoints)
    mu = Xmat %*% Bmat
  }
  if(scenario=="C"){
    # four VC with binary variable, one DA
    # short group difference
    Xmat = matrix(rbinom(n_subjects*4, 1, 0.5), n_subjects, 4)
    f1raw = function(t) 1/(1+exp((0.8-t)*20))
    f1 = function(t) ifelse(abs(f1raw(t)) < 0.05, 0, f1raw(t)) # make small values exact 0s
    Bmat = matrix(c(f1(t0), 0*t0, 0*t0, 0*t0), 4, n_timepoints, byrow=T)
    mu = Xmat %*% Bmat
  }
  if(scenario=="D"){
    # four VC with binary variable, all DA
    # long group difference
    Xmat = matrix(rbinom(n_subjects*4, 1, 0.5), n_subjects, 4)
    f1raw = function(t) 1-1/(1+exp((0.8-t)*20))
    f1 = function(t) ifelse(abs(f1raw(t)) < 0.05, 0, f1raw(t)) # make small values exact 0s
    Bmat = matrix(c(f1(t0), -f1(1-t0), f1(1-t0), -f1(t0)), 4, n_timepoints, byrow=T)
    mu = Xmat %*% Bmat
  }
  # generate responses
  errormat = matrix(stats::rnorm(n_timepoints*n_subjects), n_subjects, n_timepoints) * sqrt(observation_variance)
  ymat = mu * effect_size + thetamat + errormat
  smat = matrix(seq(n_subjects), n_subjects, n_timepoints)
  # missing data
  if(missingness=="uniform") omat = matrix(stats::runif(n_timepoints*n_subjects) < prop_observed, n_subjects, n_timepoints)
  if(missingness=="contiguous"){
    n_missing = stats::rbinom(n_subjects, n_timepoints, 1-prop_observed)
    starting = sapply(n_missing, function(x) sample(1:(n_timepoints-x+1), 1))
    ending = starting + n_missing - 1
    omat = matrix(T, n_subjects, n_timepoints)
    for (i in 1:n_subjects) if(n_missing[i]>0) omat[i, starting[i]:ending[i]] = F
  }
  if(missingness=="sqrt"){
    obs_t = sample.int( n_timepoints, round(n_timepoints*sqrt(1-prop_observed)),F)
    obs_s = sample.int( n_subjects, round(n_subjects*sqrt(1-prop_observed)),F)
    omat = matrix(T, n_subjects, n_timepoints)
    omat[obs_s, obs_t] = F
  }
  if(missingness=="fixed_uniform"){
    n_obs = round(n_timepoints*prop_observed)
    omat = matrix(F, n_subjects, n_timepoints)
    for (i in 1:n_subjects) omat[i, sample.int(n_timepoints, n_obs)] = T
  }
  # build data frame
  data_full = data.frame(
    response=as.vector(ymat),
    time=as.vector(timemat),
    subject_id=as.vector(smat),
    observed=as.vector(omat)
  )
  for(j in seq(ncol(Xmat))){
    data_full[[paste0("X", j)]] = as.vector(Xmat[,j])
  }
  data_long = data_full %>% dplyr::filter(observed) %>% dplyr::select(-observed)

  # return
  true_values = list(
    B=Bmat*effect_size,
    da=abs(Bmat*effect_size)>1e-6,
    scenario=scenario,
    effect_size=effect_size,
    random_effect_variance_ratio=random_effect_variance_ratio,
    random_effect_ar1_correlation=random_effect_ar1_correlation,
    observation_variance=observation_variance,
    prop_observed=prop_observed,
    missingness=missingness,
    seed=seed,
    n_subjects=n_subjects,
    n_timepoints=n_timepoints
  )
  instance = list(
    data_long=data_long,
    time=t0,
    covariates=paste0("X", seq(ncol(Xmat))),
    true_values=true_values
  )
  return(instance)
}

addProblem(
  name="synthetic",
  fun=synthetic,
  data=NULL
)

# # for debugging
# instance = synthetic(
#   NULL, NULL,
#   "C",
#   n_subjects=100,
#   prop_observed=0.7,
#   observation_variance=1.,
#   random_effect_ar1_correlation=1.,
#   random_effect_variance_ratio=1.,
#   effect_size=1.,
#   n_timepoints=11,
#   missingness="sqrt",
#   seed=2
# )
# ggplot() +
#   theme_minimal() +
#   geom_line(data=instance$data_long, mapping=aes(x=time, y=response, group=subject_id, color=X1)) +
#   geom_smooth(data=instance$data_long, mapping=aes(x=time, y=response, group=X1, color=X1))
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup algorithms
lsvcmm_wrapper = function(
    data, job, instance
){
  # extract data
  df = instance$data_long
  t0 = instance$time
  covariates = instance$covariates
  # arguments
  k_args = list(name="gaussian", scale=NULL, n_scale=10L, rescale_boundary=T)
  wc_args = list(name="compound_symmetry", estimate=T, ratio=1.)
  p_args = list(name="adaptive_sparse_group_lasso",
                adaptive=0.5, alpha=1., penalize_intercept=T,
                lambda=NULL, nlambda=100L)
  # fit model
  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates=covariates,
    nvc_covariates=NULL,
    offset=NULL,
    add_intercept=T,
    estimated_time=instance$estimated_time,
    kernel=k_args,
    working_covariance=wc_args,
    penalty=p_args,
    return_models=F
  )
  # model selection
  i = which.min(fit$results$ebich)
  l = fit$results$penalty.lambda[i]
  h = fit$results$kernel.scale[i]
  B = fit$vc_path[,,i]
  # bootstrap
  p_args$lambda = l
  p_args$nlambda = 1L
  k_args$scale = h
  k_args$n_scale = 1L
  boot = LSVCMM::lsvcmm.boot(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates=covariates,
    nvc_covariates=NULL,
    offset=NULL,
    add_intercept=T,
    estimated_time=instance$estimated_time,
    kernel=k_args,
    working_covariance=wc_args,
    penalty=p_args,
    n_samples=1000
  )

  bands = LSVCMM:::confidence_band(boot, var=2:(1+length(covariates)))
  bands %<>% mutate(var=paste0("X", var-1))
  bands$true_values = as.vector(t(instance$true_values$B))
  bands %<>% mutate(in_band=(true_values > L) & (true_values < U))
  bands %<>% left_join(
    bands %>% group_by(var) %>% summarize(var_in_band=all(in_band)),
    by="var"
  )
  bands$all_in_band = all(bands$in_band)
  # omnibus p-values
  fisher = function(pval) 1 - pchisq(-2*sum(log(pval)), 2*length(pval), lower.tail=TRUE)
  omni_pvalues_pervar = bands %>% group_by(var) %>%
    summarize(
      normal_min=min(pval_normal),
      normal_fisher=fisher(pval_normal),
      percentile_min=min(pval_percentile),
      percentile_fisher=fisher(pval_percentile),
      percentile=mean(pval_joint)
    )
  omni_pvalues = bands %>%
    summarize(
      normal_min=min(pval_normal),
      normal_fisher=fisher(pval_normal),
      percentile_min=min(pval_percentile),
      percentile_fisher=fisher(pval_percentile),
      percentile=mean(pval_joint_all)
    )

  # return
  fit$bands = bands
  fit$omni_pvalues = omni_pvalues
  fit$omni_pvalues_pervar = omni_pvalues_pervar
  return(fit)
}

addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_wrapper
)
# ------------------------------------------------------------------------------



# ==============================================================================
# Experimental design
n_reps=1000
problems = list(
  `synthetic`=CJ(
    scenario=c("A", "B", "C", "D"),
    effect_size=c(0., 0.25, 0.5, 0.75, 1.),
    seed=seq(n_reps)
  )
)

algorithms = list(
  `LSVCMM`=data.table()
)

addExperiments(
  prob.designs=problems,
  algo.designs=algorithms,
  repls=1
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Run
summarizeExperiments()
submitJobs(resources=list(walltime=10000))
getStatus()
# ------------------------------------------------------------------------------






# ==============================================================================
# Load batchtools registry
library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

name = "experiment_pvalues"
DIR = paste0("./", name, "/")
DIR_REGISTRY = paste0("./", name, "/registry/")

registry = loadRegistry(
  file.dir=DIR_REGISTRY,
  writeable=T
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Results
estimate = function(result) result$omni_pvalues
estimates = reduceResultsList(fun = estimate) %>% bind_rows(.id="job.id")
estimates %<>% mutate(job.id = as.numeric(job.id))

parameters = getJobPars() %>% unwrap()
parameters %<>% mutate(job.id = as.numeric(job.id))

estimates %<>% left_join(parameters, by="job.id")
estimates %>%
  group_by(scenario, effect_size) %>%
  summarize(
    normal_min=mean(normal_min<0.05),
    normal_fisher=mean(normal_fisher<0.05),
    percentile_min=mean(percentile_min<0.05),
    percentile_fisher=mean(percentile_fisher<0.05),
    percentile=mean(percentile<0.05)
  )
# ------------------------------------------------------------------------------
