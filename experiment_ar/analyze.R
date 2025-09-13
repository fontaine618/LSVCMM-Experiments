library(ggplot2)
library(dplyr)
library(magrittr)

dir = "experiment_ar"

experiments = seq(0, 1, length.out=5)

estimates = read.csv(paste0(dir, "/results/estimates.csv"))
estimation_errors = read.csv(paste0(dir, "/results/estimation_errors.csv"))
classifications = read.csv(paste0(dir, "/results/classifications.csv"))
parameters = read.csv(paste0(dir, "/results/parameters.csv"))
results = read.csv(paste0(dir, "/results/settings.csv"))
parameters %<>% select(job.id, seed, random_effect_ar1_correlation, algorithm, ar1.correlation)

estimates %<>% left_join(parameters, by="job.id")
estimation_errors %<>% left_join(parameters, by="job.id")
classifications %<>% left_join(parameters, by="job.id")
results %<>% left_join(parameters, by="job.id")

xvar = "ar1.correlation"
xdisplay = "AR(1) Correlation"

colors = c(
  "LSVCMM"="red",
  "ALasso"="lightblue",
  "LSVCM"="lightgreen",
  "SPFDA"="bisque3"
)

gs = list()

for(col in seq_along(experiments)){
  corr = experiments[col]

  # plot histogram of selected alpha
  df = results %>% filter(random_effect_ar1_correlation==!!corr, algorithm=="LSVCMM")
  # compute how often a value of ar1.correlation is within 2 of the ebich per seed
  dfmin = df %>% group_by(seed) %>% filter(ebich < min(ebich) + 2) %>% ungroup()


  g = ggplot() +
    geom_histogram(data=dfmin, aes(x=ar1.correlation), fill="grey", color="black", breaks=seq(-0.0625, 1.0625, 0.125)) +
    labs(x=xdisplay, y=ifelse(col==1, "Freq. lowest EBIC", "")) +
    theme_minimal() + labs(x=NULL) +
    ggtitle(paste0("Corr.=", corr)) +
    xlim(-0.0625, 1.3125)
  gs[[length(gs)+1]] = g

  # plot MAE of the estimates
  df = estimation_errors %>% filter(random_effect_ar1_correlation==!!corr)
  df %<>% mutate(
    ar1.correlation=ifelse(algorithm=="LSVCMM", ar1.correlation, ifelse(algorithm=="SPFDA", 1.125, 1.25)),
  )
  df = df %>% group_by(seed, ar1.correlation, algorithm) %>% summarise(mae=mean(abs(group_difference))) %>% ungroup()
  g = ggplot() +
    geom_boxplot(data=df, aes(x=ar1.correlation, y=mae, group=ar1.correlation, fill=algorithm), outliers=F) +
    labs(x=xdisplay, y=ifelse(col==1, "MAE", "")) +
    theme_minimal() +
    expand_limits(y=0) +
    # scale_y_log10() +
    theme(legend.position="none") +
    scale_x_continuous(breaks=c(0, 0.5, 1, 1.125, 1.25), labels=c("0", "0.5", "1", "SPFDA", "ALasso")) +
    scale_fill_manual(values=c("LSVCMM"="red", "SPFDA"="bisque3", "ALasso"="lightblue")) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + labs(x=NULL)
  gs[[length(gs)+1]] = g

  # plot classification error
  df = classifications %>% filter(random_effect_ar1_correlation==!!corr)
  df %<>% mutate(
    ar1.correlation=ifelse(algorithm=="LSVCMM", ar1.correlation, ifelse(algorithm=="SPFDA", 1.125, 1.25)),
  )
  df = df %>% group_by(seed, ar1.correlation, algorithm) %>%
    summarise(
      tn=sum(group_difference=="TN"),
      fp=sum(group_difference=="FP"),
      fn=sum(group_difference=="FN"),
      tp=sum(group_difference=="TP")
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
    geom_boxplot(data=df, aes(x=ar1.correlation, y=accuracy, group=ar1.correlation, fill=algorithm), outliers=F) +
    labs(x=xdisplay, y=ifelse(col==1, "Accuracy", "")) +
    theme_minimal() +
    expand_limits(y=1) +
    theme(legend.position="none") +
    scale_x_continuous(breaks=c(0, 0.5, 1, 1.125, 1.25), labels=c("0", "0.5", "1", "SPFDA", "ALasso")) +
    scale_fill_manual(values=c("LSVCMM"="red", "SPFDA"="bisque3", "ALasso"="lightblue")) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_boxplot(data=df, aes(x=ar1.correlation, y=recall, group=ar1.correlation, fill=algorithm), outlier.alpha=0.2) +
  #   labs(x=xdisplay, y=ifelse(col==1, "Recall", "")) +
  #   theme_minimal() +
  #   theme(legend.position="none") +
  #   scale_x_continuous(breaks=c(0, 0.5, 1, 1.125, 1.25), labels=c("0", "0.5", "1", "SPFDA", "ALasso")) +
  #   scale_fill_manual(values=c("LSVCMM"="red", "SPFDA"="bisque3", "ALasso"="lightblue"))
  # gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_boxplot(data=df, aes(x=ar1.correlation, y=precision, group=ar1.correlation, fill=algorithm), outlier.alpha=0.2) +
  #   labs(x=xdisplay, y=ifelse(col==1, "Precision", "")) +
  #   theme_minimal() +
  #   theme(legend.position="none") +
  #   scale_x_continuous(breaks=c(0, 0.5, 1, 1.125, 1.25), labels=c("0", "0.5", "1", "SPFDA", "ALasso")) +
  #   scale_fill_manual(values=c("LSVCMM"="red", "SPFDA"="bisque3", "ALasso"="lightblue"))
  # gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_smooth(data=df, aes(x=penalty.alpha, y=tp), method="loess") +
  #   labs(x=xdisplay, y=ifelse(col==1, "TP", "")) +
  #   theme_minimal() +
  #   expand_limits(y=1)
  # gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_smooth(data=df, aes(x=penalty.alpha, y=fp), method="loess") +
  #   labs(x=xdisplay, y=ifelse(col==1, "FP", "")) +
  #   theme_minimal() +
  #   expand_limits(y=1)
  # gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_smooth(data=df, aes(x=penalty.alpha, y=tn), method="loess") +
  #   labs(x=xdisplay, y=ifelse(col==1, "TN", "")) +
  #   theme_minimal() +
  #   expand_limits(y=1)
  # gs[[length(gs)+1]] = g
  # g = ggplot() +
  #   geom_smooth(data=df, aes(x=penalty.alpha, y=fn), method="loess") +
  #   labs(x=xdisplay, y=ifelse(col==1, "FN", "")) +
  #   theme_minimal() +
  #   expand_limits(y=1)
  # gs[[length(gs)+1]] = g
}

g = cowplot::plot_grid(plotlist=gs, ncol=length(experiments), byrow=F, align="v")
ggsave("experiment_ar/sim_ar.pdf", g, width=10, height=6)
