library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggpattern)
library(latex2exp)

# ==============================================================================
# Prepare data
DIR_FIGURES = paste0("./dmbt1/figures/")
source("./dmbt1/prepare_data.R") # only adds the pseq object, which is at the otu level
if(taxa_are_rows(pseq)) pseq = t(pseq)
prevalent_otus = microbiome::core(pseq, detection=0, prevalence=0.05) %>% phyloseq::taxa_names()
pseq %<>% microbiome::transform(transform="clr")
otus = pseq %>% phyloseq::taxa_names()
pseq %<>% phyloseq::subset_taxa(otus %in% prevalent_otus)
t0 = c(0, 4, 8, 12, 16, 22)
otus = pseq %>% phyloseq::taxa_names()
rm(prevalent_otus)
tax = phyloseq::tax_table(pseq)
clr = phyloseq::otu_table(pseq) %>% data.frame()
meta = phyloseq::sample_data(pseq)
data = list(
  clr=clr,
  tax=tax,
  meta=meta,
  t0=t0,
  otus=otus
)
# ------------------------------------------------------------------------------



# ==============================================================================
# Plot missing data

# plot 1: missing data per KO and SCC
df = bind_cols(data.frame(data$meta), data$clr[, 1])
colnames(df) = c("subject_id", "Week", "Type", "Gender", "Diagnosis", "SCC", "CLR")
df %<>% rename(Sex=Gender, Genotype=Type)
df_wide = df %>% spread(Week, CLR, sep="_")
n_observed = df_wide %>% select(-subject_id, -Sex, -Diagnosis) %>%
  group_by(SCC, Genotype) %>%
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(Week, "n", -SCC, -Genotype) %>%
  rename(Diagnosis=SCC) %>%
  mutate(Week=gsub("Week_", "", Week) %>% as.numeric())
n_total = df_wide %>% select(-subject_id, -Sex, -Diagnosis) %>%
  group_by(SCC, Genotype) %>%
  summarise_all(funs(n())) %>%
  gather(Week, "n", -SCC, -Genotype) %>%
  rename(Diagnosis=SCC) %>%
  mutate(Week=gsub("Week_", "", Week) %>% as.numeric())

gmeta = ggplot() +
  theme_minimal() +
  geom_bar_pattern(
    data=n_total,
    mapping=aes(x=Week, y=n, fill=Genotype, pattern=Diagnosis,
                alpha="Missing",
                group=paste0(Genotype, Diagnosis, Week)),
    stat="identity",
    position="dodge",
    color="black",
    width=3, linewidth=0.5,
    pattern_fill="black",
    pattern_angle=45,
    pattern_density=0.5,
    pattern_spacing=0.02,
    pattern_size=0,
    pattern_alpha=0.5
  ) +
  geom_bar_pattern(
    data=n_observed,
    mapping=aes(x=Week, y=n, fill=Genotype, pattern=Diagnosis,
                alpha="Observed",
                group=paste0(Genotype, Diagnosis, Week)),
    stat="identity",
    position="dodge",
    color="black",
    width=3, linewidth=0.5,
    pattern_fill="black",
    pattern_angle=45,
    pattern_density=0.5,
    pattern_spacing=0.02,
    pattern_size=0,
    pattern_alpha=0.5
  ) +
  ylab("Nb. samples") +
  xlab("Week") +
  scale_x_continuous(breaks=t0, minor_breaks=NULL) +
  theme(
    legend.position="right",
    text=element_text(family="Helvetica"),
    plot.title=element_text(size=12)
  ) +
  scale_pattern_manual(
    values=c("SCC"="stripe", "HP/CIS"="none"),
    labels=c("SCC"="SCC", "HP/CIS"="ED/CIS"),
    name="Diagnosis"
  ) +
  scale_fill_manual(
    c("KO", "WT"),
    values=c("#FFCB05", "#00274C"),
    aesthetics=c("color", "fill"),
    name="DMBT1",
    guide = guide_legend(override.aes = list(pattern = "none"))
  ) +
  scale_alpha_manual(
    c("Missing", "Observed"),
    values=c(0.2, 1),
    name="Samples",
    guide = guide_legend(override.aes = list(pattern = "none", color="black", fill="black"))
  ) +
  ggtitle("")

# plot 2: patterns
df_pattern = df_wide %>%
  mutate(across(starts_with("Week_"), function(x) 1*is.na(x))) %>%
  unite("pattern", starts_with("Week_"), remove=F, sep="_")
df_pattern %<>%
  group_by(pattern) %>%
  summarize(n=n()) %>%
  tidyr::separate_wider_delim(pattern, "_", names=paste0("Week_", t0)) %>%
  mutate(pattern_id=seq(n()))
df_pattern_n = df_pattern %>% select(pattern_id, n)
df_pattern %<>%
  tidyr::pivot_longer(paste0("Week_", t0), names_to="Week", values_to="value") %>%
  mutate(Week=gsub("Week_", "", Week) %>% as.numeric())

gpatterns = ggplot() +
  theme_minimal() +
  geom_tile(
    data=df_pattern,
    mapping=aes(x=as.factor(Week), y=as.factor(pattern_id), fill=value)
  ) +
  scale_fill_manual(
    values=c("0"="#000000", "1"="#bbbbbb"),
    labels=c("0"="Observed", "1"="Missing")
    ) +
  scale_y_discrete(
    breaks=df_pattern_n$pattern_id,
    labels=df_pattern_n$n
  ) +
  theme(
    legend.position="right",
    text=element_text(family="Helvetica"),
    plot.title=element_text(size=12)
  ) +
  labs(fill="", x="Week", y="Pattern occurences") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# merge
g = egg::ggarrange(gmeta, gpatterns, ncol=2, widths=c(2, 1))
ggsave(paste0(DIR_FIGURES, "missingness.pdf"), g, width=10, height=3)
# ------------------------------------------------------------------------------





# ==============================================================================
# Example
otu = "Otu0014"
df = bind_cols(data.frame(data$meta), data$clr[, otu])
colnames(df) = c("Mouse", "Week", "Type", "Gender", "Diagnosis", "SCC", "CLR")
df %<>% rename(Sex=Gender, Genotype=Type)
df %<>% mutate(KO=ifelse(Genotype=="KO", 1, 0))

fit = LSVCMM::lsvcmm(
  data=df,
  response="CLR",
  subject="Mouse",
  time="Week",
  vc_covariates=c("KO"),
  add_intercept=T,
  penalty=list(alpha=0.5, adaptive=0.5, penalize_intercept=F),
  kernel=list(scale=0.2),
  return_models=FALSE
)

best_per_h = fit$results %>%
  mutate(model_id=seq(1, nrow(.))) %>%
  select(ebich, kernel.scale, penalty.lambda, model_id, df) %>%
  group_by(kernel.scale) %>%
  slice(which.min(ebich)) %>%
  ungroup() %>%
  mutate(overall_best=ebich==min(ebich))

ggplot(
  mapping=aes(x=penalty.lambda, y=ebich, color=as.factor(kernel.scale), group=kernel.scale)
) +
  geom_line(data=fit$results) +
  geom_point(data=best_per_h, mapping=aes(shape=overall_best), size=2) +
  scale_x_log10() +
  theme_bw() +
  xlab(TeX("Regularization parameter ($\\lambda$)")) +
  ylab("EBIC") +
  scale_shape_manual(values=c(1, 19)) +
  guides(
    color=guide_legend(title=TeX("Kernel scale ($h$)")),
    shape=guide_legend(title="Overall min.")
  )

B = fit$vc_path[,,best_per_h$model_id]

# EBIC selection
i = which.min(fit$results$ebich)
l = fit$results$penalty.lambda[i]
h = fit$results$kernel.scale[i]

# Bootstrap
boot = LSVCMM::lsvcmm.boot(
  data=df,
  response="CLR",
  subject="Mouse",
  time="Week",
  vc_covariates=c("KO"),
  add_intercept=T,
  kernel=list(scale=h),
  penalty=list(alpha=0.5, adaptive=0.5, penalize_intercept=F, lambda=l),
  n_samples=1000
)
LSVCMM:::confidence_band(boot, var=2)
grp_mean = matrix(c(
  1, 0,
  1, 1
), 2, 2, T)
boot_grp_mean = LSVCMM:::transform_obj(boot, grp_mean)
grp_means = list(
  WT=LSVCMM:::confidence_band(boot_grp_mean, var=1),
  KO=LSVCMM:::confidence_band(boot_grp_mean, var=2)
) %>% bind_rows(.id="Genotype")

yrange = range(df$CLR)

gdata = ggplot() +
  theme_minimal() +
  geom_line(
    data=df,
    mapping=aes(x=Week, y=CLR, group=Mouse, color=Genotype),
    alpha=0.5
  ) +
  ylab("CLR-transformed abundance") +
  xlab("Week") +
  scale_x_continuous(breaks=t0, minor_breaks=NULL) +
  theme(
    legend.position="none",
    text=element_text(family="Helvetica"),
    plot.title=element_text(size=12)
  ) +
  scale_color_manual(
    c("KO", "WT"),
    values=c("#FFCB05", "#00274C"),
    name="Genotype"
  ) +
  ggtitle(paste0(otu, ": Observed trajectories")) +
  ylim(yrange[1], yrange[2])

est_mean = data.frame(
  Week=t0,
  WT=B[1,],
  KO=B[1,]+B[2,]
) %>% gather("Genotype", "est_mean", -Week)

gmean = ggplot() +
  theme_minimal() +
  geom_line(
    data=grp_means,
    mapping=aes(x=estimated_time, y=median, group=Genotype, color=Genotype),
    alpha=1, linewidth=1.5
  ) +
  geom_ribbon(
    data=grp_means,
    mapping=aes(x=estimated_time, ymin=L, ymax=U, group=Genotype, fill=Genotype),
    alpha=0.5
  ) +
  ylab("Estimated mean") +
  xlab("Week") +
  scale_x_continuous(breaks=t0, minor_breaks=NULL) +
  theme(
    legend.position="right",
    text=element_text(family="Helvetica"),
    plot.title=element_text(size=12)
  ) +
  scale_color_manual(
    c("KO", "WT"),
    values=c("#FFCB05", "#00274C"),
    name="DMBT1",
    aesthetics=c("color", "fill")
  ) +
  ggtitle(paste0(otu, ": Fitted values")) +
  ylim(yrange[1], yrange[2])
g = cowplot::plot_grid(gdata, gmean, ncol=2, nrow=1, rel_widths=c(4, 5), align="hv", axis="tbl")
ggsave(paste0(DIR_FIGURES, otu, ".pdf"), g, width=8, height=2.5)


# plot with boxplots to compare
ghist = ggplot() +
  theme_minimal() +
  geom_boxplot(
    data=df,
    mapping=aes(x=Week, y=CLR, group=paste(Genotype, Week), color=Genotype)
  )

# ------------------------------------------------------------------------------


