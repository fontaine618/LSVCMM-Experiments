library(ggplot2)
library(dplyr)
library(magrittr)

dir = "experiment_sgl"

experiments = list(
  "Dense, group dense"=c(effect_sparsity=-1.5, effect_groupsparsity=0.),
  "Sparse, group dense"=c(effect_sparsity=0.8, effect_groupsparsity=0.),
  "Dense, group sparse"=c(effect_sparsity=-1.5, effect_groupsparsity=0.8),
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

for(col in seq_along(experiments)){
  exp = names(experiments)[col]
  effect_sparsity = experiments[[exp]][1]
  effect_groupsparsity = experiments[[exp]][2]

  # plot true values
  instance = LSVCMM::generate_synthetic_data_p(
    n_subjects=200,
    n_timepoints=21,
    n_features=5,
    feature_type="binary",
    observation_variance=1.,
    random_effect_variance_ratio=1.,
    random_effect_ar1_correlation=1.,
    effect_size=1.,
    effect_groupsparsity=effect_groupsparsity, # all 5 or just 1
    effect_sparsity=effect_sparsity, # narrow spike or large bump
    prop_observed=0.7,
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
    labs(title=exp, x="Time", y=ifelse(col==1, "True VCs", "")) +
    theme(legend.position="none")
  gs[[length(gs)+1]] = g

  # plot estimated curves
  df = estimates %>% filter(effect_sparsity==!!effect_sparsity, effect_groupsparsity==!!effect_groupsparsity)
  # wide to long
  df = df %>% tidyr::pivot_longer(cols=3:7, names_to="feature", values_to="estimate")
  df = df %>% group_by(time, feature) %>% summarise(
    median=median(estimate),
    lower=quantile(estimate, 0.025),
    upper=quantile(estimate, 0.975)
    ) %>% ungroup()
  g = ggplot() +
    geom_ribbon(data=df, aes(x=time, ymin=lower, ymax=upper, fill=feature), alpha=0.1) +
    geom_line(data=df, aes(x=time, y=median, color=feature, group=feature)) +
    # geom_smooth(data=df, aes(x=time, y=estimate, color=feature, group=feature), method="gam") +
    theme_minimal() +
    coord_cartesian(ylim=c(-1.1, 1.1)) +
    labs(x="Time", y=ifelse(col==1, "Estimated VCs", "")) +
    theme(legend.position="none")
  gs[[length(gs)+1]] = g

  # plot histogram of selected alpha
  df = results %>% filter(effect_sparsity==!!effect_sparsity, effect_groupsparsity==!!effect_groupsparsity)
  # subset to rows with ebich within 1 of the best ebich per seed
  dfmin = df %>% group_by(seed) %>% filter(ebich < min(ebich) + 2) %>% ungroup()
  g = ggplot() +
    geom_histogram(data=dfmin, aes(x=penalty.alpha), fill="grey", color="black", breaks=seq(-0.05, 1.05, 0.1)) +
    labs(x=xdisplay, y=ifelse(col==1, "Freq. lowest EBIC", "")) +
    theme_minimal()
  gs[[length(gs)+1]] = g

  # plot MAE of the estimates
  df = estimation_errors %>% filter(effect_sparsity==!!effect_sparsity, effect_groupsparsity==!!effect_groupsparsity)
  df = df %>% group_by(seed, penalty.alpha) %>% summarise(mae=mean(abs(X1) + abs(X2) + abs(X3) + abs(X4) + abs(X5))/51) %>% ungroup()
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=mae), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "MAE", "")) +
    theme_minimal() +
    expand_limits(y=0)
  gs[[length(gs)+1]] = g

  # plot classification error
  df = classifications %>% filter(effect_sparsity==!!effect_sparsity, effect_groupsparsity==!!effect_groupsparsity)
  df = df %>% group_by(seed, penalty.alpha) %>%
    summarise(
      tn=sum(X1=="TN") + sum(X2=="TN") + sum(X3=="TN") + sum(X4=="TN") + sum(X5=="TN"),
      fp=sum(X1=="FP") + sum(X2=="FP") + sum(X3=="FP") + sum(X4=="FP") + sum(X5=="FP"),
      fn=sum(X1=="FN") + sum(X2=="FN") + sum(X3=="FN") + sum(X4=="FN") + sum(X5=="FN"),
      tp=sum(X1=="TP") + sum(X2=="TP") + sum(X3=="TP") + sum(X4=="TP") + sum(X5=="TP")
    ) %>%
    mutate(
      precision=tp/(tp+fp),
      recall=tp/(tp+fn),
      f1=2*precision*recall/(precision+recall),
      accuracy=(tp+tn)/(tp+tn+fp+fn),
      fprrate=fp/(fp+tn),
      fnrrate=fn/(fn+tp),
      fdrrate=fp/(fp+tp)
    )
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=accuracy), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "Accuracy", "")) +
    theme_minimal() +
    expand_limits(y=1)
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=recall), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "Recall", "")) +
    theme_minimal()
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=precision), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "Precision", "")) +
    theme_minimal()
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=tp), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "TP", "")) +
    theme_minimal() +
    expand_limits(y=1)
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=fp), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "FP", "")) +
    theme_minimal() +
    expand_limits(y=1)
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=tn), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "TN", "")) +
    theme_minimal() +
    expand_limits(y=1)
  gs[[length(gs)+1]] = g
  g = ggplot() +
    geom_smooth(data=df, aes(x=penalty.alpha, y=fn), method="loess") +
    labs(x=xdisplay, y=ifelse(col==1, "FN", "")) +
    theme_minimal() +
    expand_limits(y=1)
  gs[[length(gs)+1]] = g
}

g = cowplot::plot_grid(plotlist=gs, ncol=length(experiments), byrow=F, align="v")
ggsave("sim_sgl.pdf", g, width=10, height=16)
