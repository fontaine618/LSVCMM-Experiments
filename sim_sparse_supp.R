library(ggplot2)
library(dplyr)
library(magrittr)

experiments = list(
  snr=list(
    dir="snr100",
    name="(a) Signal strength",
    xvar="observation_variance",
    xname="Variance",
    col=1,
    transform="none",
    ref=1.
  ),
  sparse=list(
    dir="missing100",
    name="(b) Sparsity",
    xvar="prop_observed",
    xname="Proportion observed",
    col=2,
    transform="none",
    ref=0.1
  ),
  # cov=list(
  #   dir="ar1_100",
  #   name="(c) Cov. misspecification",
  #   xvar="random_effect_ar1_correlation",
  #   xname="RE AR(1) correlation",
  #   col=3,
  #   transform="none",
  #   ref=1.
  # ),
  re=list(
    dir="re_ratio100",
    name="(c) RE size",
    xvar="random_effect_variance_ratio",
    xname="RE variance ratio",
    col=3,
    transform="sqrt",
    ref=1.
  )
)

gs = list()
colors = c(
  "LSVCMM"="red",
  "ALasso"="lightblue",
  "LSVCM"="lightgreen",
  "SPFDA"="bisque3"
)

display_names = c(
  "LSVCMM"="LSVCMM",
  "LSVCMM.Cross-sectional"="ALasso",
  "LSVCMM.Independent"="LSVCM",
  "SPFDA"="SPFDA"
)

fn = "sigmoid"

for(exp in experiments){
  estimates = read.csv(paste0("experiment_", exp$dir, "/results/estimates.csv"))
  estimation_errors = read.csv(paste0("experiment_", exp$dir, "/results/estimation_errors.csv"))
  classifications = read.csv(paste0("experiment_", exp$dir, "/results/classifications.csv"))
  parameters = read.csv(paste0("experiment_", exp$dir, "/results/parameters.csv"))

  # patch names
  parameters$algorithm = display_names[parameters$algorithm]

  estimates %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
  estimation_errors %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
  classifications %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)

  gby1 = c("algorithm", "seed", exp$xvar)
  gby2 = c("algorithm", exp$xvar)

  # MAE Null
  df = estimation_errors %>%
    filter(time < 0.45) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae, na.rm=T), sd=sd(mae, na.rm=T), se=sd(mae, na.rm=T)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    # geom_vline(xintercept=exp$ref, linetype="dashed", color="grey") +
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
    xlab(exp$xname) + ylab("MAE (null)") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ggtitle(exp$name) +
    ylim(0, 0.2)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "mae_null")]] = g

  # MAE Non-Null
  df = estimation_errors %>%
    filter(time > 0.45) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference))) %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(mae, na.rm=T), sd=sd(mae, na.rm=T), se=sd(mae, na.rm=T)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    # geom_vline(xintercept=exp$ref, linetype="dashed", color="grey") +
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
    xlab(exp$xname) + ylab("MAE (non-null)") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0, 0.8)
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
      ppv=(tp+fp)/(tp+tn+fp+fn),
      acc=(tp+tn)/(tp+tn+fp+fn),
      fdr=fp/pmax(tp+fp, 1),
      tpr=tp/pmax(tp+fn, 1),
    )

  # prop selected
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(ppv), sd=sd(ppv), se=sd(ppv)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    # geom_vline(xintercept=exp$ref, linetype="dashed", color="grey") +
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
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0, 1)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "prop_selected")]] = g

  # fdr
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(fdr), sd=sd(fdr), se=sd(fdr)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    # geom_vline(xintercept=exp$ref, linetype="dashed", color="grey") +
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
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0, 0.5)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "fdr")]] = g

  # tpr
  df = dfall %>%
    group_by(across(all_of(gby2))) %>%
    summarise(
      mean=mean(tpr), sd=sd(tpr), se=sd(tpr)/sqrt(n())
    )
  g = ggplot() +
    theme_minimal() +
    # geom_vline(xintercept=exp$ref, linetype="dashed", color="grey") +
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
    xlab(exp$xname) + ylab("TPR") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0, 1)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "tpr")]] = g

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
  scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
  xlab(exp$xname) + ylab("FDR") +
  labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
  theme(legend.direction="horizontal")
glegend = cowplot::get_legend(gtmp)
glegend = ggpubr::as_ggplot(glegend)

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(experiments), nrow=5,
  byrow=F,
  align="none", axis="tblr",
  rel_widths=c(1, rep(0.9, length(experiments)-1)),
  rel_heights=c(1, 0.8, 0.8, 0.8, 1)
)


gg = cowplot::plot_grid(g, glegend, ncol=1, nrow=2, rel_heights=c(20, 1))

ggsave(paste0("./sim_sparse_supp.pdf"), gg, width=length(experiments)*3+1, height=15)

