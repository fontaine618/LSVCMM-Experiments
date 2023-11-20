library(tidyverse)
# setup
n_timepoints = 31
t0 = seq(0, 1, length.out=n_timepoints)
fn = function(t, center) 1 / (1 + exp(-100 * (t - center)))
center = 0.2
b = fn(t0, center)

shifts = seq(-100, 100)/50

# penalties
lasso = abs
scad = function(x, lambda=0.5, a=3.7) {
  x = abs(x)
  out = lambda*x*(x <= lambda) +
    (2*a*lambda*x - x^2 - lambda^2) / (2*(a-1))*(x <= a*lambda & x > lambda) +
    lambda^2*(a+1)/2*(x > a*lambda)
  return(out)
}
bridge = function(x, power = 0.5) abs(x) ^ power

lasso = data.frame(
  name="lasso",
  shift=shifts,
  penalty=outer(b, shifts, "+") %>% lasso %>% colSums()
)

scad0.1 = data.frame(
  name="scad(l=0.1)",
  shift=shifts,
  penalty=10*(outer(b, shifts, "+") %>% scad(., 0.1) %>% colSums())
)

scad0.5 = data.frame(
  name="scad(l=0.5)",
  shift=shifts,
  penalty=outer(b, shifts, "+") %>% scad(., 0.5) %>% colSums()
)

bridge0.5 = data.frame(
  name="bridge0.5",
  shift=shifts,
  penalty=outer(b, shifts, "+") %>% bridge(., 0.5) %>% colSums()
)

bridge0.9 = data.frame(
  name="bridge0.9",
  shift=shifts,
  penalty=outer(b, shifts, "+") %>% bridge(., 0.9) %>% colSums()
)

# combine
penalties = bind_rows(lasso, scad0.1, scad0.5, bridge0.5, bridge0.9)

# get shift and penalty at minimum
min_penalty = penalties %>%
  group_by(name) %>%
  filter(penalty == min(penalty)) %>%
  ungroup()
# plot
g = penalties %>%
  ggplot(aes(shift, penalty, color=name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(data=min_penalty, aes(shift, penalty, color=name), size=3) +
  ggtitle("Mode=1")
# ggtitle(paste0("Mode=", round(median(b), 2)))
ggsave("./penalty/penalty_mode1.pdf", g, width=6, height=4)

# get estimates at minimum
lasso_b = data.frame(
  name="lasso",
  t=t0,
  b=b + min_penalty$shift[min_penalty$name == "lasso"]
)
scad0.1_b = data.frame(
  name="scad(l=0.1)",
  t=t0,
  b=b + min_penalty$shift[min_penalty$name == "scad(l=0.1)"]
)
scad0.5_b = data.frame(
  name="scad(l=0.5)",
  t=t0,
  b=b + min_penalty$shift[min_penalty$name == "scad(l=0.5)"]
)
bridge0.5_b = data.frame(
  name="bridge0.5",
  t=t0,
  b=b + min_penalty$shift[min_penalty$name == "bridge0.5"]
)
bridge0.9_b = data.frame(
  name="bridge0.9",
  t=t0,
  b=b + min_penalty$shift[min_penalty$name == "bridge0.9"]
)
true = data.frame(
  name="true",
  t=t0,
  b=b
)

# combine
penalties_b = bind_rows(lasso_b, scad0.1_b, scad0.5_b, bridge0.5_b, bridge0.9_b)

# plot
glasso = lasso_b %>%
  ggplot(aes(t, b, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(data=true, aes(t, b), size=1, color="black") +
  ggtitle(paste0("Best Lasso = ", round(min_penalty$penalty[min_penalty$name == "lasso"], 2),
                 " (penalty at truth = ", round(lasso$penalty[lasso$shift == 0], 2), ")"))
gscad0.1 = scad0.1_b %>%
  ggplot(aes(t, b, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(data=true, aes(t, b), size=1, color="black") +
  ggtitle(paste0("Best SCAD (l=0.1) = ", round(min_penalty$penalty[min_penalty$name == "scad(l=0.1)"], 2),
                 " (penalty at truth = ", round(scad0.1$penalty[scad0.1$shift == 0], 2), ")"))
gscad0.5 = scad0.5_b %>%
  ggplot(aes(t, b, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(data=true, aes(t, b), size=1, color="black") +
  ggtitle(paste0("Best SCAD (l=0.5) = ", round(min_penalty$penalty[min_penalty$name == "scad(l=0.5)"], 2),
                 " (penalty at truth = ", round(scad0.5$penalty[scad0.5$shift == 0], 2), ")"))
gbridge0.5 = bridge0.5_b %>%
  ggplot(aes(t, b, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(data=true, aes(t, b), size=1, color="black") +
  ggtitle(paste0("Best Bridge (l=0.5) = ", round(min_penalty$penalty[min_penalty$name == "bridge0.5"], 2),
                 " (penalty at truth = ", round(bridge0.5$penalty[bridge0.5$shift == 0], 2), ")"))
gbridge0.9 = bridge0.9_b %>%
  ggplot(aes(t, b, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(data=true, aes(t, b), size=1, color="black") +
  ggtitle(paste0("Best Bridge (l=0.9) = ", round(min_penalty$penalty[min_penalty$name == "bridge0.9"], 2),
                 " (penalty at truth = ", round(bridge0.9$penalty[bridge0.9$shift == 0], 2), ")"))
# mermge
gridExtra::grid.arrange(glasso, gscad0.1, gscad0.5, gbridge0.5, gbridge0.9, ncol=2)




# LSVCMM experiments
library(tidyverse)
n_timepoints = 31
t0 = seq(0, 1, length.out=n_timepoints)
med0 = function(t) 1 / (1 + exp(-20 * (t - 0.8)))
med0.5 = function(t) 1 / (1 + exp(-20 * (t - 0.5)))
med1 = function(t) 1 / (1 + exp(-20 * (t - 0.2)))
b0 = med0(t0)
b0.2 = med0.2(t0)
b0.5 = med0.5(t0)
b1 = med1(t0)



re = 0.5

out = lapply(fns, function(fn){

  # generate data
  instance = LSVCMM::generate_synthetic_data(
    n_subjects=100,
    n_timepoints=31,
    prop_observed=1.0,
    observation_variance=1.,
    random_effect_variance_ratio=re,
    random_effect_ar1_correlation=1.,
    effect_size=1.,
    grpdiff_function=fn,
    seed=1
  )
  df = instance$data

  out_true = data.frame(
    name="truth",
    t=t0,
    b=instance$true_values$b1
  )

  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    kernel=list(scale=0.2, name="gaussian"),
    penalty=list(adaptive=1., n_lambda=100, penalize_intercept=T, lambda_factor=1e-4),
    control=list(max_iter=1000, max_rounds=20, rel_tol=1e-6, update_method="PGD", verbose=0),
    working_covariance=list(name="independent")
  )

  i = which.min(fit$results$ebich)
  b_ind = fit$vc_path[2,,i]
  out_ind = data.frame(
    name="independent",
    t=t0,
    b=b_ind
  )


  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    kernel=list(scale=0.2, name="gaussian"),
    penalty=list(adaptive=1., n_lambda=100, penalize_intercept=T, lambda_factor=1e-4),
    control=list(max_iter=1000, max_rounds=20, rel_tol=1e-6, update_method="PGD", verbose=0),
    working_covariance=list(estimate=T, ratio=1, name="compound_symmetry")
  )

  i = which.min(fit$results$ebich)
  b_re = fit$vc_path[2,,i]
  out_are = data.frame(
    name="re (adaptive)",
    t=t0,
    b=b_re
  )


  fit = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    kernel=list(scale=0.2, name="gaussian"),
    penalty=list(adaptive=0., n_lambda=100, penalize_intercept=T, lambda_factor=1e-4),
    control=list(max_iter=1000, max_rounds=20, rel_tol=1e-6, update_method="PGD", verbose=0),
    working_covariance=list(estimate=T, ratio=1, name="compound_symmetry")
  )

  i = which.min(fit$results$ebich)
  b_re = fit$vc_path[2,,i]
  out_re = data.frame(
    name="re",
    t=t0,
    b=b_re
  )


  fithttp://127.0.0.1:25001/graphics/plot_zoom_png?width=1200&height=744 = LSVCMM::lsvcmm(
    data=df,
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    kernel=list(scale=0.2, name="gaussian"),
    penalty=list(lambda=c(0.)),
    control=list(max_iter=1000, max_rounds=20, rel_tol=1e-6, update_method="PGD", verbose=0),
    working_covariance=list(estimate=T, ratio=1, name="compound_symmetry")
  )

  i = which.min(fit$results$ebich)
  b_re = fit$vc_path[2,,i]
  out_ure = data.frame(
    name="re (unpenalized)",
    t=t0,
    b=b_re
  )


  out = bind_rows(out_ind, out_re, out_are, out_ure, out_true)
}) %>% bind_rows(.id="median")

# plot results: one plot per value of median
# x axis: t
# y axis: b
# curves and color: name

gs = lapply(names(fns), function(fn) {
  df_fn = out %>% filter(median==fn)
  g = ggplot() +
    geom_line(data=df_fn, aes(t, b, color=name)) +
    theme_bw() +
    ggtitle(fn) +
    geom_hline(yintercept=0, color="black", linetype="dashed")
  return(g)
})

g = cowplot::plot_grid(plotlist=gs, ncol=1)
# save
ggsave(paste0("./penalty/comparison_re_", re, ".pdf"), g, width=6, height=8, units="in")


