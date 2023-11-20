library(ggplot2)
library(dplyr)
library(magrittr)

experiments = list(
  missing=list(
    dir="missing",
    name="Partially observed trajectories",
    xvar="prop_observed",
    xname="Proportion observed",
    col=1
  ),
  cov=list(
    dir="n_timepoints",
    name="Nb. of timepoints",
    xvar="n_timepoints",
    xname="Nb. of timepoints",
    col=2
  ),
  cov=list(
    dir="re_ratio",
    name="RE size",
    xvar="random_effect_variance_ratio",
    xname="RE variance ratio",
    col=3
  ),
  cov=list(
    dir="n_subjects",
    name="Nb. subjects",
    xvar="n_subjects",
    xname="Nb. subjects",
    col=4
  ),
  cov=list(
    dir="ar1",
    name="Misspecified working cov.",
    xvar="random_effect_ar1_correlation",
    xname="RE AR(1) correlation",
    col=5
  ),
  cov=list(
    dir="snr100",
    name="Signal strength (100, sparse)",
    xvar="observation_variance",
    xname="Variance",
    col=6
  ),
  cov=list(
    dir="missing100large",
    name="Partially observed trajectories (100)",
    xvar="prop_observed",
    xname="Proportion observed",
    col=7
  ),
  cov=list(
    dir="ar1_100",
    name="Misspecified working cov. (100)",
    xvar="random_effect_ar1_correlation",
    xname="RE AR(1) correlation",
    col=8
  ),
  cov=list(
    dir="re_ratio100",
    name="RE size (100)",
    xvar="random_effect_variance_ratio",
    xname="RE variance ratio",
    col=9
  )
  # cov=list(
  #   dir="missing20",
  #   name="Partially observed trajectories (20)",
  #   xvar="prop_observed",
  #   xname="Proportion observed",
  #   col=8
  # )
  # cov=list(
  #   dir="smalln_propobserved",
  #   name="Partially observed trajectories (n=20)",
  #   xvar="prop_observed",
  #   xname="Proportion observed",
  #   col=9
  # ),
  # cov=list(
  #   dir="n_subjects100",
  #   name="Nb. subjects (100)",
  #   xvar="n_subjects",
  #   xname="Nb. subjects",
  #   col=10
  # )
)

fn = "sigmoid"
# fn = "sine"

gs = list()

for(exp in experiments){
  estimates = read.csv(paste0("experiment_", exp$dir, "/results/estimates.csv"))
  estimation_errors = read.csv(paste0("experiment_", exp$dir, "/results/estimation_errors.csv"))
  classifications = read.csv(paste0("experiment_", exp$dir, "/results/classifications.csv"))
  parameters = read.csv(paste0("experiment_", exp$dir, "/results/parameters.csv"))

  estimates %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
  estimation_errors %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
  classifications %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)

  gby1 = c("algorithm", "seed", exp$xvar)
  gby2 = c("algorithm", exp$xvar)

  # MAE
  df = estimation_errors %>%
    # filter(time>=0.5) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae, na.rm=T), sd=sd(mae, na.rm=T), se=sd(mae, na.rm=T)/sqrt(n())
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
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    ) +
    ggtitle(exp$name) +
    ylim(0, 0.5)
  # if(exp$col>1) g = g + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank(),
  #   axis.title.y=element_blank()
  # )
  gs[[paste0(exp$dir, "mae")]] = g


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
      ppv=(tp+fp)/(tp+tn+fp+fn),
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
    xlab(exp$xname) + ylab("Prop. selected") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    ) +
    ylim(0, 1)
  # if(exp$col>1) g = g + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank(),
  #   axis.title.y=element_blank()
  # )
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
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    ) +
    ylim(0.4, 1)
  # if(exp$col>1) g = g + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank(),
  #   axis.title.y=element_blank()
  # )
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
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    ) +
    ylim(0, 1)
  # if(exp$col>1) g = g + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank(),
  #   axis.title.y=element_blank()
  # )
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
    ylim(0, 0.5)
  # if(exp$col>1) g = g + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank(),
  #   axis.title.y=element_blank()
  # )
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

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(experiments), nrow=5,
  byrow=F,
  align="vh", axis="tblr"
)

gg = cowplot::plot_grid(glegend, g, glegend, ncol=1, nrow=3, rel_heights=c(1, 20, 1))

ggsave(paste0("./", fn, "31_results.pdf"), gg, width=length(experiments)*4+1, height=16)

