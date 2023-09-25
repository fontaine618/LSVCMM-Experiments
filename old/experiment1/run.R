library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

# ==============================================================================
# Setup batchtools registry
DIR = "./experiment1/batchtools/"

registry = makeExperimentRegistry(
  file.dir=DIR,
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM", "spfda")
)
registry$cluster.functions = makeClusterFunctionsSocket(ncpus = 10)
export = list()
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
  name="LSVCMM.Non-adaptive",
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
    random_effect_ar1_correlation=1.,
    random_effect_variance_ratio=2.,
    prop_observed=0.5
  )
)

algorithms = list(
  `LSVCMM`=data.table(cross_sectional=F, independent=F, penalty.adaptive=1.),
  `LSVCMM.Independent`=data.table(cross_sectional=F, independent=T, penalty.adaptive=1.),
  `LSVCMM.Cross-sectional`=data.table(cross_sectional=T, independent=T, penalty.adaptive=1.),
  `LSVCMM.Non-adaptive`=data.table(cross_sectional=F, independent=F, penalty.adaptive=0.),
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
submitJobs(resources=list(walltime=1000))
getStatus()
# ------------------------------------------------------------------------------





#
#
#
#
# # ==============================================================================
# # Gather results
# instance = synthetic(NULL, NULL)
# t0 = instance$true_values$time
# tnames = paste0("T_", t0)
#
# estimate = function(result) result
# estimates = reduceResultsDataTable(fun = estimate) %>% unwrap()
# colnames(estimates) = c("job.id", paste0("T_", t0))
#
# parameters = getJobPars() %>% unwrap()
#
# estmat = estimates %>% apply(1, function(row) row) %>% t() %>% data.frame()
# estmat %<>% left_join(parameters %>% select(job.id, algorithm, seed), by="job.id")
# estmat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="estimate")
# estmat %<>% mutate(time=as.numeric(stringr::str_split_i(time, "_", 2)))
# estmat %<>% filter(algorithm!="Cross-sectional")
# estmat %<>% mutate(algorithm=ifelse(algorithm=="Lasso", "Cross-sectional", algorithm))
#
# errmat = estimates %>% apply(1, function(row) abs(row - c(0, instance$true_values$b1))) %>% t() %>% data.frame()
# errmat %<>% left_join(parameters %>% select(job.id, algorithm, seed), by="job.id")
# errmat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="estimate")
# errmat %<>% mutate(time=as.numeric(stringr::str_split_i(time, "_", 2)))
# errmat %<>% filter(algorithm!="Cross-sectional")
# errmat %<>% mutate(algorithm=ifelse(algorithm=="Lasso", "Cross-sectional", algorithm))
# # ------------------------------------------------------------------------------
#
#
#
# # ==============================================================================
# # Plot
# theme_set(theme_minimal())
#
# gest = ggplot() +
#   geom_line(
#     data=instance$true_values,
#     mapping=aes(x=time, y=b1),
#     linewidth=1, alpha=1
#   ) +
#   geom_boxplot(
#     data=estmat,
#     mapping=aes(x=time, y=estimate,
#                 color=algorithm,
#                 fill=algorithm,
#                 group=paste(algorithm, time)),
#     outlier.alpha=0.05, alpha=0.8
#   ) + ylab("Estimate") + xlab("Time") +
#   scale_x_continuous(breaks=0:5/5, labels=0:5/5) +
#   theme(
#     legend.position="none",
#     text=element_text(family="Helvetica"),
#     # panel.grid = element_line(color = "white"),
#     panel.grid = element_line(color = "black"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.x=element_blank(),
#     panel.grid.minor.y = element_blank()
#   ) + ylim(-0.5, 1.5) +
#   scale_fill_manual(
#     values=c("LSVCMM"="#FFCB05", "Independent"="#000000", "SPFDA"="#00274C", "Cross-sectional"="#aaaaaa"),
#     name="Algorithm",
#     aesthetics=c("fill", "color")
#   )
#
#
# gerr = ggplot() +
#   geom_boxplot(
#     data=errmat,
#     mapping=aes(x=time, y=estimate,
#                 color=algorithm,
#                 fill=algorithm,
#                 group=paste(algorithm, time)),
#     outlier.alpha=0.05, alpha=0.8
#   ) + ylab("Absolute\nerror") + xlab("Time") +
#   scale_x_continuous(breaks=0:5/5, labels=0:5/5) +
#   theme(
#     legend.position="none",
#     text=element_text(family="Helvetica"),
#     # panel.grid = element_line(color = "white"),
#     panel.grid = element_line(color = "black"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.x=element_blank(),
#     panel.grid.minor.y = element_blank()
#   ) + ylim(0., 1.1) +
#   scale_fill_manual(
#     values=c("LSVCMM"="#FFCB05", "Independent"="#000000", "SPFDA"="#00274C", "Cross-sectional"="#aaaaaa"),
#     name="Algorithm",
#     aesthetics=c("fill", "color")
#   )
#
# gprop = ggplot() +
#   geom_bar(
#     data=estmat %>% group_by(algorithm, time) %>% summarize(propnz=mean(estimate != 0)),
#     mapping=aes(x=time, y=propnz,
#                 color=algorithm,
#                 fill=algorithm,
#                 group=paste(algorithm, time)),
#     alpha=0.8, stat="identity", position="dodge"
#   ) + ylab("Proportion\nselected") + xlab("Time") +
#   scale_x_continuous(breaks=0:5/5, labels=0:5/5) +
#   theme(
#     legend.position="bottom",
#     text=element_text(family="Helvetica"),
#     # panel.grid = element_line(color = "white"),
#     panel.grid = element_line(color = "black"),
#     axis.title.x=element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank()
#   ) + labs(fill="Method", color="Method") +
#   scale_fill_manual(
#     values=c("LSVCMM"="#FFCB05", "Independent"="#000000", "SPFDA"="#00274C", "Cross-sectional"="#aaaaaa"),
#     name="Method",
#     aesthetics=c("fill", "color")
#   )
#
#
# g = cowplot::plot_grid(gest, gerr, gprop,
#                        ncol=1, nrow=3, align="vh", axis="lr", rel_heights=c(1, 0.8, 1.4))
# ggsave(paste0(DIR, "results.pdf"), g, width=6, height=3.5)
# # ------------------------------------------------------------------------------
