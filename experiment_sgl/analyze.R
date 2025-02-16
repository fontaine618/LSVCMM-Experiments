mse = estimation_errors %>% select(-time) %>% group_by(job.id) %>%
  summarize(
    mse = mean(mean(X1^2) + mean(X2^2) + mean(X3^2) + mean(X4^2) + mean(X5^2))
  )
parms = parameters %>% select(job.id, penalty.alpha, effect_groupsparsity, effect_sparsity)
out = mse %>% left_join(parms)
