library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment1/batchtools/"

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
decision = estimates %>% mutate(
  group_difference=ifelse(abs(group_difference)>1e-6, 1, 0),
  intercept=ifelse(abs(intercept)>1e-4, 1, 0)
)


parameters = getJobPars() %>% unwrap()
estimates %<>% left_join(parameters, by="job.id")
estimation_errors %<>% left_join(parameters, by="job.id")

# ------------------------------------------------------------------------------



# ==============================================================================
# Plot results
ggplot() +
  theme_minimal() +
  geom_line(
    data=estimation_errors %>% group_by(time, algorithm) %>% summarise(group_difference=mean(group_difference^2)),
    aes(x=time, y=group_difference, color=algorithm, group=paste0(algorithm)),
    alpha=1.
  )

ggplot() +
  theme_minimal() +
  geom_boxplot(
      data=estimation_errors,
      mapping=aes(x=time,
                  y=(group_difference)^2,
                  color=algorithm,
                  fill=algorithm,
                  group=paste(algorithm, time)),
      outlier.alpha=0.05, alpha=0.8
    )
# ------------------------------------------------------------------------------

