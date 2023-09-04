library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment_missing/"
if(!dir.exists(DIR)) dir.create(DIR, recursive=T)
registry = makeExperimentRegistry(
  file.dir=paste0(DIR, "registry/"),
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM", "spfda")
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

# for debugging
instance = synthetic(NULL, NULL)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup algorithms
source("./algorithms/lsvcmm.R")
source("./algorithms/spfda.R")
addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="LSVCMM.Independent",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="LSVCMM.Cross-sectional",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="SPFDA",
  fun=spfda_wrapper
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Experimental design
n_reps=100
problems = list(
  `synthetic`=CJ(
    seed=seq(n_reps),
    observation_variance=1.,
    random_effect_ar1_correlation=1.,
    random_effect_variance_ratio=2.,
    prop_observed=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.)
  )
)

algorithms = list(
  `LSVCMM`=data.table(cross_sectional=F, independent=F, penalty.adaptive=1.),
  `LSVCMM.Independent`=data.table(cross_sectional=F, independent=T, penalty.adaptive=1.),
  `LSVCMM.Cross-sectional`=data.table(cross_sectional=T, independent=T, penalty.adaptive=1.),
  `SPFDA`=data.table()
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
  account="stats_dept1",
  partition="standard",
  memory="6g", # this is per cpu
  ncpus=1,
  walltime="60:00",
  chunks.as.arrayjobs=FALSE,
  job_name="LSVCMM_missing"
)

chunk_df = data.table(job.id=1:2800, chunk=1:100)
head(chunk_df)
submitJobs(chunk_df, resources)
# ------------------------------------------------------------------------------
