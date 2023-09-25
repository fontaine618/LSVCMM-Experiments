library(ggplot2)
library(dplyr)
library(magrittr)

experiments = list(
  # snr=list(
  #   dir="snr",
  #   name="Signal strength",
  #   xvar="observation_variance",
  #   xname="Total variance",
  #   col=1
  # ),
  missing=list(
    dir="missing",
    name="Partially observed trajectories",
    xvar="prop_observed",
    xname="Proportion observed",
    col=2
  ),
  # re=list(
  #   dir="re",
  #   name="Random effect importance",
  #   xvar="random_effect_variance_ratio",
  #   xname="Random effect variance ratio",
  #   col=3
  # ),
  # cov=list(
  #   dir="cov_fixed",
  #   name="Covariance misspecification",
  #   xvar="random_effect_ar1_correlation",
  #   xname="AR(1) correlation",
  #   col=1
  # )
)

gs = list()

for(exp in experiments){
  estimates = read.csv(paste0("experiment_", exp$dir, "/estimates.csv"))
  estimation_errors = read.csv(paste0("experiment_", exp$dir, "/estimation_errors.csv"))
  classifications = read.csv(paste0("experiment_", exp$dir, "/classifications.csv"))

  estimation_errors %<>% mutate(
    Dependency=ifelse(grepl("Independent", algorithm),"Independent","Longitudinal"),
    Lambda=ifelse(grepl("l_selected", algorithm),"Selected","Fixed"),
  )

  classifications %<>% mutate(
    Dependency=ifelse(grepl("Independent", algorithm),"Independent","Longitudinal"),
    Lambda=ifelse(grepl("l_selected", algorithm),"Selected","Fixed"),
  )

  gby1 = c("algorithm", "Dependency", "Lambda", "seed", exp$xvar)
  gby2 = c("algorithm", "Dependency", "Lambda", exp$xvar)

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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    ) +
    ggtitle(exp$name)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae")]] = g

  # MAE Null
  df = estimation_errors %>%
    filter(time<0.5) %>%
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE (null)") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    )
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae_null")]] = g

  # MAE Non-null
  df = estimation_errors %>%
    filter(time>=0.5) %>%
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Estimation MAE (non-null)") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    )
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("N selected") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    )
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("Accuracy") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    )
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("TPR (Power)") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
    )
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
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
    ) +
    geom_point(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
    ) +
    geom_ribbon(
      data=df,
      mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
      alpha=0.2
    ) +
    xlab(exp$xname) + ylab("FDR") +
    labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.x=element_blank(),
    )
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
    mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, linetype=Lambda),
  ) +
  geom_point(
    data=df,
    mapping=aes(x=!!sym(exp$xvar), y=mean, color=Dependency, shape=Lambda),
  ) +
  geom_ribbon(
    data=df,
    mapping=aes(x=!!sym(exp$xvar), ymin=mean-se, ymax=mean+se, fill=Dependency, group=algorithm),
    alpha=0.2
  ) +
  xlab(exp$xname) + ylab("FDR") +
  labs(color="Dependency", linetype="Lambda", shape="Lambda", fill="Dependency") +
  theme(legend.direction="vertical")
glegend = cowplot::get_legend(gtmp)
glegend = ggpubr::as_ggplot(glegend)

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(experiments), nrow=7,
  byrow=F,
  align="vh", axis="tblr"
)

gg = cowplot::plot_grid(glegend, g, glegend, ncol=1, nrow=3, rel_heights=c(1, 10, 1))

ggsave("./experiment_fixed_results.pdf", gg, width=6, height=20)

