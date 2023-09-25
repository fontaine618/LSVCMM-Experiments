# investigating drop in accuracy for sigmoid experiment with prop observed

estimates = read.csv(paste0("experiment_missing_sine/estimates.csv"))
classifications = read.csv(paste0("experiment_missing_sine/classifications.csv"))
settings = read.csv(paste0("experiment_re_sine/settings.csv"))
parameters = read.csv(paste0("experiment_missing_sine/parameters.csv"))

df = classifications %>%
  filter(algorithm=="LSVCMM") %>%
  filter(observation_variance==0.5)

ggplot() +
  geom_histogram(
    data=df,
    mapping=aes(x=time, fill=group_difference),
    bins=31
  ) +
  theme_minimal()

ggplot() +
  geom_boxplot(
    data=settings,
    mapping=aes(x=as.factor(observation_variance), y=kernel.scale, fill=algorithm)
  ) +
  scale_y_log10() +
  theme_minimal()


ggplot() +
  geom_boxplot(
    data=settings,
    mapping=aes(x=as.factor(observation_variance), y=penalty.lambda, fill=algorithm)
  ) +
  scale_y_log10() +
  theme_minimal()
