library(ggplot2)
library(dplyr)
library(magrittr)

dir = "misspecification"
fn = "sine"

display_names = c(
  "LSVCMM"="LSVCMM",
  "LSVCMM.Cross-sectional"="ALasso",
  "LSVCMM.Independent"="LSVCM",
  "SPFDA"="SPFDA"
)

estimates = read.csv(paste0("experiment_", dir, "/results/estimates.csv"))
estimation_errors = read.csv(paste0("experiment_", dir, "/results/estimation_errors.csv"))
classifications = read.csv(paste0("experiment_", dir, "/results/classifications.csv"))
parameters = read.csv(paste0("experiment_", dir, "/results/parameters.csv"))

# patch names
parameters$algorithm = display_names[parameters$algorithm]

estimates %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
estimation_errors %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)
classifications %<>% left_join(parameters, by="job.id") %>% filter(grpdiff_function==fn)

# experiment setup
random_effect_ar1_correlation=c(0., 0.3, 0.7, 1.)  # axis
random_effect_variance_ratio=c(1., 2., 4.)  # columns
xvar = "random_effect_ar1_correlation"
xname = "AR(1) correlation"
gs = list()
colors = c(
  "LSVCMM"="red",
  "ALasso"="lightblue",
  "LSVCM"="lightgreen",
  "SPFDA"="bisque3"
)

for(col in seq_along(random_effect_variance_ratio)){
  re_ratio = random_effect_variance_ratio[col]
  name = paste0("RE Variance ratio = ", re_ratio)

  gby1 = c("algorithm", "seed", xvar)

  # MAE
  df = estimation_errors %>%
    filter(random_effect_variance_ratio==re_ratio) %>%
    group_by(across(all_of(gby1))) %>%
    summarise(mae=mean(abs(group_difference)))
  g = ggplot() +
    theme_minimal() +
    geom_boxplot(
      data=df,
      mapping=aes(x=as.factor(random_effect_ar1_correlation), y=mae, fill=algorithm,
                  group=paste(algorithm, random_effect_ar1_correlation)),
      width=0.5
    ) +
    xlab(xname) + ylab("MAE") +
    labs(fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors) +
    ggtitle(name) +
    ylim(0, 0.4)
  if(col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(re_ratio, "mae")]] = g


  # Classification metrics

  df = classifications %>%
    filter(random_effect_variance_ratio==re_ratio) %>%
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

  # accuracy
  g = ggplot() +
    theme_minimal() +
    geom_boxplot(
      data=df,
      mapping=aes(x=as.factor(random_effect_ar1_correlation), y=acc, fill=algorithm,
                  group=paste(algorithm, random_effect_ar1_correlation)),
      width=0.5
    ) +
    xlab(xname) + ylab("Accuracy") +
    labs(fill="Algorithm") +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    ) +
    scale_fill_manual(values=colors) +
    ylim(0.4, 1)
  if(col>1) g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  gs[[paste0(re_ratio, "acc")]] = g

}

gtmp = ggplot() +
  theme_minimal() +
  geom_boxplot(
    data=df,
    mapping=aes(x=as.factor(random_effect_ar1_correlation), y=acc, fill=algorithm,
                group=paste(algorithm, random_effect_ar1_correlation)),
    width=0.5
  ) +
  scale_fill_manual(values=colors) +
  labs(fill="Algorithm") +
  theme(legend.direction="horizontal")
glegend = cowplot::get_legend(gtmp)
glegend = ggpubr::as_ggplot(glegend)

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(random_effect_variance_ratio), nrow=2,
  byrow=F,
  align="none", axis="tblr",
  rel_widths=c(1, rep(0.9, length(random_effect_variance_ratio)-1)),
  rel_heights=c(1, 1)
)


gg = cowplot::plot_grid(g, glegend, ncol=1, nrow=2, rel_heights=c(10, 1))

ggsave(paste0("./sim_misspecification.pdf"), gg, width=length(experiments)*3+1, height=6)

