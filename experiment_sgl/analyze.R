library(ggplot2)
library(dplyr)
library(magrittr)

dir = "experiment_sgl"

experiments = list(
  "Dense, group dense"=c(effect_sparsity=0., effect_groupsparsity=0.),
  "Dense, group sparse"=c(effect_sparsity=0., effect_groupsparsity=0.8),
  "Sparse, group dense"=c(effect_sparsity=0.8, effect_groupsparsity=0.),
  "Sparse, group sparse"=c(effect_sparsity=0.8, effect_groupsparsity=0.8)
)

estimates = read.csv(paste0(dir, "/results/estimates.csv"))
estimation_errors = read.csv(paste0(dir, "/results/estimation_errors.csv"))
classifications = read.csv(paste0(dir, "/results/classifications.csv"))
parameters = read.csv(paste0(dir, "/results/parameters.csv"))
results = read.csv(paste0(dir, "/results/settings.csv"))

estimates %<>% left_join(parameters, by="job.id")
estimation_errors %<>% left_join(parameters, by="job.id")
classifications %<>% left_join(parameters, by="job.id")
results %<>% left_join(parameters %>% select(job.id, seed, effect_sparsity, effect_groupsparsity), by="job.id")

xvar = "penalty.alpha"
xdisplay = expression(alpha)

# y vars:
# - ebich: plot all curves
# - MAE of the estimates
# - classification error

gs = list()

for(exp in names(experiments)){
  effect_sparsity = experiments[[exp]][1]
  effect_groupsparsity = experiments[[exp]][2]

  # plot true values
  instance = LSVCMM::generate_synthetic_data_p(
    n_subjects=20,
    n_timepoints=51,
    n_features=5,
    feature_type="binary",
    observation_variance=1.,
    random_effect_variance_ratio=1.,
    random_effect_ar1_correlation=1.,
    effect_size=1.,
    effect_groupsparsity=effect_groupsparsity, # all 5 or just 1
    effect_sparsity=effect_sparsity, # narrow spike or large bump
    prop_observed=0.5,
    missingness="uniform",
    seed=1
  )
  B = instance$true_values[, -2] # drop intercept
  Bdf = as.data.frame(B)
  # wide to long X1 ... X5
  Bdf = Bdf %>% tidyr::pivot_longer(cols=2:6, names_to="feature", values_to="value")
  g = ggplot() +
    geom_line(data=Bdf, aes(x=time, y=value, color=feature)) +
    theme_minimal() +
    ylim(-1.1, 1.1) +
    labs(title=exp, x="Time", y="True effects") +
    theme(legend.position="none")
  gs[[length(gs)+1]] = g

  # plot histogram of selected alpha
  df = results %>% filter(effect_sparsity==!!effect_sparsity, effect_groupsparsity==!!effect_groupsparsity)
  # subset to rows with ebich within 1 of the best ebich per seed
  dfmin = df %>% group_by(seed) %>% filter(ebich < min(ebich) + 2) %>% ungroup()
  g = ggplot() +
    geom_histogram(data=dfmin, aes(x=penalty.alpha), fill="grey", color="black", breaks=seq(-0.05, 1.05, 0.1)) +
    labs(title=exp, x=xdisplay, y="Frequency lowest EBIC") +
    theme_minimal()
  gs[[length(gs)+1]] = g
}

cowplot::plot_grid(plotlist=gs, ncol=length(experiments), byrow=F)
