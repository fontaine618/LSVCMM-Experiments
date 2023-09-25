library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment_missing/registry/"

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


write.csv(estimates, file="./experiment_missing/estimates.csv", row.names=F)
write.csv(estimation_errors, file="./experiment_missing/estimation_errors.csv", row.names=F)
write.csv(classifications, file="./experiment_missing/classifications.csv", row.names=F)

setting = function(result){
  i = which.min(result$fit$results[["ebich"]])
  res = result$fit$results[i, ]
  return(res)
}
settings = reduceResultsList(fun = setting) %>% bind_rows(.id="job.id")
settings %<>% mutate(job.id = as.numeric(job.id))
settings %<>% left_join(parameters, by="job.id")
write.csv(settings, file="./experiment_missing/settings.csv", row.names=F)
# ------------------------------------------------------------------------------
