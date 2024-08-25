library(ggplot2)
library(tidyverse)
library(magrittr)

name = "experiment_pvalues"
DIR_RESULTS = paste0("./", name, "/results/")
DIR_FIGURES = paste0("./", name, "/")

# get omni p-values results
pvalues = read.csv(paste0(DIR_RESULTS, "omni_pvalues.csv"))
pvalues_long = pvalues %>%
  pivot_longer(
    cols = -c(scenario, effect_size),
    names_to = "method",
    values_to = "pvalue"
    ) %>% filter(method!="percentile")
pvalues_long %<>% mutate(
  # extract up to underscore
  pvalue_method=stringr::str_extract(method, ".*(?=_)"),
  # extract after underscore
  pvalue_combination=stringr::str_extract(method, "(?<=_).*")
)
n = 1000
pvalues_long %<>% mutate(
  pvale_se = sqrt(pvalue*(1-pvalue)/n)
)

# plot omnibus p-values
nrows = 2
ncols = 2
scenarios = list(
  A="A. Single VC, 4/11 DA",
  B="B. Single VC, 10/11 DA",
  C="C. Four VCs, One 4/11 DA",
  D="D. Four VCs, All 10/11 DA"
)
methods_display = c(
  percentile="Percentile",
  normal="Normal"
)
combination_display = c(
  "min"="Min",
  "fisher"="Fisher"
)
methods_colors = c(
  percentile="darkblue",
  normal="red"
)
combination_linestyle = c(
  "min"="solid",
  "fisher"="dashed"
)
xaxis_display = "Effect Size"
yaxis_display = "Proportion positive"

gs = list()
for (scenario in names(scenarios)) {
  pvalues_long_scenario = pvalues_long %>% filter(scenario==!!scenario)
  g = ggplot() +
    theme_minimal() +
    geom_hline(yintercept=0.05, linetype="dotted") +
    geom_line(
      data=pvalues_long_scenario,
      aes(x=effect_size, y=pvalue, color=pvalue_method, linetype=pvalue_combination, group=method),
      linewidth=1
      ) +
    # geom_ribbon(
    #   data=pvalues_long_scenario,
    #   aes(x=effect_size, ymin=pvalue-pvale_se, ymax=pvalue+pvale_se, fill=pvalue_method, group=method),
    #   alpha=0.2
    #   ) +
    labs(
      title=scenarios[[scenario]],
      x=xaxis_display,
      y=yaxis_display,
      color="Pointwise\np-value",
      # fill="p-value",
      linetype="Combination"
    ) + xlim(c(0, 1)) + ylim(c(0, 1)) +
    scale_color_manual(values=methods_colors, aesthetics=c( "color"), labels=methods_display) +
    scale_linetype_manual(values=combination_linestyle, labels=combination_display)
  # extract legend
  glegend = ggpubr::as_ggplot(cowplot::get_legend(g))
  g = g +
    theme(
      legend.position="none",
      text=element_text(family="Helvetica"),
      panel.border = element_rect(colour = "grey", fill=NA, size=1),
    )
  if(scenario %in% c("B", "D"))  g = g + theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )
  if(scenario %in% c("A", "B"))  g = g + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank()
  )
  gs[[scenario]] = g
}

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=ncols, nrow=nrows,
  byrow=T,
  align="none", axis="tblr",
  rel_widths=c(1, rep(0.9, ncols-1)),
  rel_heights=c(1, rep(1.1, nrows-1))
)
g
# add legend
gg = cowplot::plot_grid(g, glegend, ncol=2, nrow=1, rel_widths=c(8, 1))
# save
ggsave(paste0(DIR_FIGURES, "omni_pvalues.pdf"), gg, width=10, height=8)
