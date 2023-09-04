library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment_snr/registry/"

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

write.csv(estimates, file="./experiment_snr/estimates.csv", row.names=F)
write.csv(estimation_errors, file="./experiment_snr/estimation_errors.csv", row.names=F)
write.csv(classifications, file="./experiment_snr/classifications.csv", row.names=F)
# ------------------------------------------------------------------------------
