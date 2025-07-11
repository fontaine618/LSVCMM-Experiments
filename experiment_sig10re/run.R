library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)


# ==============================================================================
# Setup batchtools registry

setwd("/storage/work/spf5519/LSVCMM/LSVCMM-Experiments")
env_path = "/storage/work/spf5519/LSVCMM/renv/activate.R"

name = "experiment_sig10re"
DIR = paste0("./", name, "/")
DIR_REGISTRY = paste0("./", name, "/registry/")
if(dir.exists(DIR_REGISTRY)) unlink(DIR_REGISTRY, recursive=T)
if(!dir.exists(DIR)) dir.create(DIR, recursive=T)
registry = makeExperimentRegistry(
  file.dir=DIR_REGISTRY,
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM", "spfda", "splinectomeR", "SummarizedExperiment", "OmicsLonDA")
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup problem
synthetic = function(data, job, ...){
  instance = LSVCMM::generate_synthetic_data(...)
  return(instance)
}

addProblem(
  name="synthetic",
  fun=synthetic,
  data=NULL
)
#
# instance = synthetic(
#   NULL, NULL,
#   n_subjects=50,
#   prop_observed=0.5,
#   observation_variance=1.,
#   random_effect_ar1_correlation=1.,
#   random_effect_variance_ratio=1.,
#   effect_size=1.,
#   n_timepoints=10,
#   grpdiff_function="sine",
#   missingness="sqrt",
#   seed=1
# )
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup algorithms
source("./algorithms/lsvcmm_boot.R")
source("./algorithms/spfda_conf.R")
source("./algorithms/ssanova.R")
source("./algorithms/splines.R")
addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_boot_wrapper
)
addAlgorithm(
  name="SPFDA",
  fun=spfda_conf_wrapper
)
addAlgorithm(
  name="SSANOVA",
  fun=ssanova_wrapper
)
addAlgorithm(
  name="SPLINECTOMER",
  fun=splinectomer_wrapper
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Experimental design
n_reps=100
problems = list(
  `synthetic`=CJ(
    n_subjects=100,
    prop_observed=0.5,
    observation_variance=1.,
    random_effect_ar1_correlation=1.,
    random_effect_variance_ratio=seq(0, 2, 0.25),
    effect_size=1.,
    n_timepoints=10,
    grpdiff_function=c("sine"),
    missingness="sqrt",
    seed=seq(n_reps)
  )
)

algorithms = list(
  `LSVCMM`=data.table(cross_sectional=F, independent=F, penalty.adaptive=0.5, kernel.scale=0.2, selection="bich"),
  `SPFDA`=data.table(K=12),
  `SSANOVA`=data.table(),
  `SPLINECTOMER`=data.table()
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
  memory="7g", # this is per cpu
  ncpus=1,
  walltime="2:00:00",
  chunks.as.arrayjobs=FALSE,
  job_name=name,
  env=env_path
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
ids = findDone()

decision = function(result) result$decision
decisions = reduceResultsList(fun = decision) %>% bind_rows(.id="job.id")
decisions %<>% mutate(job.id = as.numeric(job.id))
decisions %<>% mutate(job.id=ids$job.id[job.id])

classification = function(result) result$classification
classifications = reduceResultsList(fun = classification) %>% bind_rows(.id="job.id")
classifications %<>% mutate(job.id = as.numeric(job.id))
classifications %<>% mutate(job.id=ids$job.id[job.id])

parameters = getJobPars() %>% unwrap()
parameters %<>% mutate(job.id = as.numeric(job.id))

write.csv(parameters, file=paste0(DIR_RESULTS, "parameters.csv"), row.names=F)
write.csv(classifications, file=paste0(DIR_RESULTS, "classifications.csv"), row.names=F)
write.csv(decisions, file=paste0(DIR_RESULTS, "decisions.csv"), row.names=F)
# ------------------------------------------------------------------------------
