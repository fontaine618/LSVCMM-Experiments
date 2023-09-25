library(ggplot2)
library(dplyr)
library(magrittr)

experiments = list(
  snr=list(
    dir="snr_sine",
    name="Signal strength",
    xvar="observation_variance",
    xname="Total variance",
    col=1
  ),
  missing=list(
    dir="missing_sine",
    name="Partially observed trajectories",
    xvar="prop_observed",
    xname="Proportion observed",
    col=2
  ),
  re=list(
    dir="re_sine",
    name="Random effect importance",
    xvar="random_effect_variance_ratio",
    xname="Random effect variance ratio",
    col=3
  )
  # cov=list(
  #   dir="cov",
  #   name="Covariance misspecification",
  #   xvar="random_effect_ar1_correlation",
  #   xname="AR(1) correlation",
  #   col=4
  # )
)

gs = list()

for(exp in experiments){
  estimates = read.csv(paste0("experiment_", exp$dir, "/estimates.csv"))
  estimation_errors = read.csv(paste0("experiment_", exp$dir, "/estimation_errors.csv"))
  classifications = read.csv(paste0("experiment_", exp$dir, "/classifications.csv"))

  gby1 = c("algorithm", "seed", exp$xvar)
  gby2 = c("algorithm", exp$xvar)

  # MAE
  df = estimation_errors %>%
    # filter(time>=0.5) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae), sd=sd(mae), se=sd(mae)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ggtitle(exp$name) +
    ylim(0, 0.5)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae")]] = g

  # MAE Null
  df = estimation_errors %>%
    filter(abs(time-0.5)>0.25) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae), sd=sd(mae), se=sd(mae)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE (null)") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ylim(0, 0.3)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae_null")]] = g

  # MAE Non-null
  df = estimation_errors %>%
    filter(abs(time-0.5)<0.25) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae), sd=sd(mae), se=sd(mae)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE (non-null)") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ylim(0, 0.75)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae_nonnull")]] = g

  # Classification metrics

  dfall = classifications %>%
    group_by(across(all_of(gby1))) %>%
    summarise(
      tn=sum(group_difference=="TN"),
      fp=sum(group_difference=="FP"),
      fn=sum(group_difference=="FN"),
      tp=sum(group_difference=="TP")
    ) %>%
    mutate(
      ppv=(tp+fp),
      acc=(tp+tn)/(tp+tn+fp+fn),
      fdr=fp/pmax(tp+fp, 1),
      tpr=tp/pmax(tp+fn, 1),
    )

  # ppv
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(ppv), sd=sd(ppv), se=sd(ppv)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("N selected") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ylim(0, 31)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "ppv")]] = g

  # accuracy
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(acc), sd=sd(acc), se=sd(acc)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Accuracy") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ylim(0.5, 1)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "acc")]] = g

  # power
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(tpr), sd=sd(tpr), se=sd(tpr)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("TPR (Power)") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ylim(0, 1)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "tpr")]] = g

  # fdr
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(fdr), sd=sd(fdr), se=sd(fdr)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    geom_line(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("FDR") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    ) +
    ylim(0, 0.3)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "fdr")]] = g

}

gtmp = ggplot() +
  theme_minimal() +
  geom_line(
    data=df,
    mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, linetype=algorithm),
  ) +
  geom_point(
    data=df,
    mapping=aes(x=!!sym(exp$xvar), y=mean, color=algorithm, shape=algorithm),
  ) +
  geom_ribbon(
    data=df,
    mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=algorithm),
    alpha=0.2
  ) +
  xlab(exp$xname) + ylab("FDR") +
  labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
  theme(legend.direction="horizontal")
glegend = cowplot::get_legend(gtmp)
glegend = ggpubr::as_ggplot(glegend)

ne = length(experiments)
g = cowplot::plot_grid(
  plotlist=gs,
  ncol=ne, nrow=7,
  byrow=F,
  align="vh", axis="tblr"
)

gg = cowplot::plot_grid(glegend, g, glegend, ncol=1, nrow=3, rel_heights=c(3, 20, 3))

ggsave("./experiment_sine_plot.pdf", gg, width=ne*3+1, height=20)


