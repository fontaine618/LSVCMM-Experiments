source("/storage/work/spf5519/LSVCMM/renv/activate.R")

library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
setwd("/storage/work/spf5519/LSVCMM/LSVCMM-Experiments")

name = "experiment_sgl"
DIR = paste0("./", name, "/")
DIR_REGISTRY = paste0("./", name, "/registry/")
if(dir.exists(DIR_REGISTRY)) unlink(DIR_REGISTRY, recursive=T)
if(!dir.exists(DIR)) dir.create(DIR, recursive=T)
registry = makeExperimentRegistry(
  file.dir=DIR_REGISTRY,
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM")
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup problem
synthetic = function(data, job, ...){
  instance = LSVCMM::generate_synthetic_data_p(...)
  return(instance)
}

addProblem(
  name="synthetic",
  fun=synthetic,
  data=NULL
)

# for debugging
instance = synthetic(NULL, NULL, n_timepoints=51, n_features=5)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup algorithms
source("./algorithms/lsvcmm.R")
addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_wrapper_p
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Experimental design
n_reps=2
problems = list(
  `synthetic`=CJ(
    n_subjects=20,
    n_timepoints=51,
    n_features=5,
    feature_type="binary",
    observation_variance=1.,
    random_effect_variance_ratio=1.,
    random_effect_ar1_correlation=1.,
    effect_size=1.,
    effect_groupsparsity=c(0., 0.8), # all 5 or just 1
    effect_sparsity=c(0., 0.8), # narrow spike or large bump
    prop_observed=0.5,
    missingness="uniform",
    seed=seq(n_reps)
  )
)

algorithms = list(
  `LSVCMM`=data.table(
    cross_sectional=F,
    independent=F,
    penalty.adaptive=0.5,
    penalty.alpha=seq(0, 1, 0.1),
    kernel.scale=0.2
    )
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
getStatus()

resources = list(
  account="open",
  partition="open",
  memory="10g", # this is per cpu
  ncpus=1,
  walltime="3:00:00",
  chunks.as.arrayjobs=FALSE,
  job_name=name
)
njobs = findJobs() %>% nrow()
chunk_df = data.table(job.id=1:njobs, chunk=1:n_reps)
head(chunk_df)
submitJobs(chunk_df, resources)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup batchtools registry
registry = loadRegistry(
  file.dir=DIR_REGISTRY,
  writeable=T
)
# ------------------------------------------------------------------------------



# ==============================================================================
# Gather results
DIR_RESULTS = paste0("./", name, "/results/")
if(!dir.exists(DIR_RESULTS)) dir.create(DIR_RESULTS, recursive=T)

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
parameters %<>% mutate(job.id = as.numeric(job.id))

setting = function(result){
  i = which.min(result$fit$results[["ebich"]])
  res = result$fit$results[i, ]
  return(res)
}
settings = reduceResultsList(fun = setting) %>% bind_rows(.id="job.id")
settings %<>% mutate(job.id = as.numeric(job.id))

write.csv(parameters, file=paste0(DIR_RESULTS, "parameters.csv"), row.names=F)
write.csv(estimates, file=paste0(DIR_RESULTS, "estimates.csv"), row.names=F)
write.csv(estimation_errors, file=paste0(DIR_RESULTS, "estimation_errors.csv"), row.names=F)
write.csv(classifications, file=paste0(DIR_RESULTS, "classifications.csv"), row.names=F)
write.csv(settings, file=paste0(DIR_RESULTS, "settings.csv"), row.names=F)
# ------------------------------------------------------------------------------
