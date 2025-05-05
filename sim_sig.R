library(ggplot2)
library(dplyr)
library(magrittr)

experiments = list(
  n10=list(
    dir="sig10",
    name="(a) Missing values (10 timepoints)",
    xvar="n_subjects",
    xname="Nb. subjects",
    col=1,
    transform="none",
    ref=100
  )
  # n100=list(
  #   dir="sig100",
  #   name="(b) Irregular sampling (100 timepoints)",
  #   xvar="n_subjects",
  #   xname="Nb. subjects",
  #   col=1,
  #   transform="none",
  #   ref=100
  # )
)

gs = list()
colors = c(
  "LSVCMM"="red",
  "ALasso"="lightblue",
  "LSVCM"="lightgreen",
  "SPFDA"="bisque3",
  "SS-ANOVA"="lightpink",
  "SplinectomeR"="khaki"
)
fn = "sine"

display_names = c(
  "LSVCMM"="LSVCMM",
  "LSVCMM.Cross-sectional"="ALasso",
  "LSVCMM.Independent"="LSVCM",
  "SPFDA"="SPFDA",
  "SSANOVA"="SS-ANOVA",
  "SPLINECTOMER"="SplinectomeR"
)



for(exp in experiments){
  decisions = read.csv(paste0("experiment_", exp$dir, "/results/decisions.csv"))
  classifications = read.csv(paste0("experiment_", exp$dir, "/results/classifications.csv"))
  parameters = read.csv(paste0("experiment_", exp$dir, "/results/parameters.csv"))

  # patch names
  parameters$algorithm = display_names[parameters$algorithm]

  decisions %<>% left_join(parameters, by="job.id")
  classifications %<>% left_join(parameters, by="job.id")

  gby1 = c("algorithm", "seed", exp$xvar)
  gby2 = c("algorithm", exp$xvar)


  # Classification metrics

  dfall = classifications %>%
    group_by(across(all_of(gby1))) %>%
    summarise(
      tn=sum(classification=="TN"),
      fp=sum(classification=="FP"),
      fn=sum(classification=="FN"),
      tp=sum(classification=="TP")
    ) %>%
    mutate(
      ppv=(tp+fp)/(tp+tn+fp+fn),
      acc=(tp+tn)/(tp+tn+fp+fn),
      fdr=fp/pmax(tp+fp, 1),
      tpr=tp/pmax(tp+fn, 1),
    )

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
    ggtitle(exp$name) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0., 0.3)
  if(exp$col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(exp$dir, "fdr")]] = g

  # power
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
    xlab(exp$xname) + ylab("Power") +
    labs(color="Algorithm", linetype="Algorithm", shape="Algorithm", fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      panel.border = element_rect(colour = "grey", fill=NA, linewidth=1),
    ) +
    scale_fill_manual(values=colors, aesthetics=c("fill", "color")) +
    ylim(0., 1)
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
  theme(
    legend.direction="horizontal",
    text=element_text(family="Helvetica")
  )
glegend = cowplot::get_legend(gtmp)
glegend = ggpubr::as_ggplot(glegend)

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(experiments), nrow=2,
  byrow=F,
  align="none", axis="tblr",
  rel_widths=c(1, rep(0.9, length(experiments)-1)),
  rel_heights=c(1, 1)
)


gg = cowplot::plot_grid(g, glegend, ncol=1, nrow=2, rel_heights=c(10, 1))

ggsave(paste0("./sim_sig.pdf"), gg, width=length(experiments)*2.5+1, height=3)
