library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment_cov_fixed/"
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
addAlgorithm(
  name="Longitudinal.h_fixed.l_selected",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="Longitudinal.h_fixed.l_fixed",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="Independent.h_fixed.l_selected",
  fun=lsvcmm_wrapper
)
addAlgorithm(
  name="Independent.h_fixed.l_fixed",
  fun=lsvcmm_wrapper
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Experimental design
n_reps=100
problems = list(
  `synthetic`=CJ(
    seed=seq(n_reps),
    observation_variance=1.,
    random_effect_ar1_correlation=c(0., 0.1, 0.2, 0.4, 0.6, 0.8, 1.),
    random_effect_variance_ratio=2.,
    prop_observed=1.
  )
)

algorithms = list(
  `Longitudinal.h_fixed.l_selected`=data.table(cross_sectional=F, independent=F, penalty.adaptive=1., kernel.scale=0.2),
  `Longitudinal.h_fixed.l_fixed`=data.table(cross_sectional=F, independent=F, penalty.adaptive=1., kernel.scale=0.2, penalty.lambda=10.),
  `Independent.h_fixed.l_selected`=data.table(cross_sectional=F, independent=T, penalty.adaptive=1., kernel.scale=0.2),
  `Independent.h_fixed.l_fixed`=data.table(cross_sectional=F, independent=T, penalty.adaptive=1., kernel.scale=0.2, penalty.lambda=10.)
)

addExperiments(
  prob.designs=problems,
  algo.designs=algorithms,
  repls=1
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Run
registry$cluster.functions = makeClusterFunctionsSocket(ncpus = 10)
summarizeExperiments()
submitJobs(resources=list(walltime=1000))
getStatus()
# ------------------------------------------------------------------------------
