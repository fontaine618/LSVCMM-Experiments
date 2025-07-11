library(batchtools)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(shadowtext)
library(ggrepel)


# ==============================================================================
# Setup batchtools registry
DIR = paste0("./dmbt1/registry/")
DIR_FIGURES = paste0("./dmbt1/figures/")
registry = loadRegistry(file.dir=DIR)
t0 = c(0, 4, 8, 12, 16, 22)
tax_level = "OTU"
# ------------------------------------------------------------------------------



# ==============================================================================
# Prepare data
DIR_FIGURES = paste0("./dmbt1/figures/")
source("./dmbt1/prepare_data.R") # only adds the pseq object, which is at the otu level
if(taxa_are_rows(pseq)) pseq = t(pseq)
prevalent_otus = microbiome::core(pseq, detection=0, prevalence=0.05) %>% phyloseq::taxa_names()
pseq_raw = pseq
pseq %<>% microbiome::transform(transform="clr")
otus = pseq %>% phyloseq::taxa_names()
pseq %<>% phyloseq::subset_taxa(otus %in% prevalent_otus)
pseq_raw %<>% phyloseq::subset_taxa(otus %in% prevalent_otus)
t0 = c(0, 4, 8, 12, 16, 22)
otus = pseq %>% phyloseq::taxa_names()
rm(prevalent_otus)
tax = phyloseq::tax_table(pseq)
clr = phyloseq::otu_table(pseq) %>% data.frame()
counts = phyloseq::otu_table(pseq_raw) %>% data.frame()
rel_abundance = counts/rowSums(counts)
meta = phyloseq::sample_data(pseq)
data = list(
  clr=clr,
  abundance=counts,
  rel_abundance=rel_abundance,
  tax=tax,
  meta=meta,
  t0=t0,
  otus=otus
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Load results
estimate = function(result) result
estimates = reduceResultsList(fun = estimate) %>% bind_rows(.id="job.id")
estimates %<>% mutate(job.id = as.numeric(job.id)) %>% rename(otu2=otu)
parameters = getJobPars() %>% unwrap()
parameters %<>% mutate(job.id = as.numeric(job.id))
estimates %<>% left_join(parameters, by="job.id")
# patch SPFDA from 1.96 to bonferonni sim. band
cval = qnorm(1-0.025/6)
estimates %<>%
  mutate(
    lower=ifelse(algo_name=="SPFDA", estimate-cval*se, lower),
    upper=ifelse(algo_name=="SPFDA", estimate+cval*se, upper),
  )
estimates %<>%
  mutate(
    da=ifelse(algo_name=="SPFDA", 1*((lower>0) | (upper<0)), da)
  )
sig_lsvcmm = estimates %>% filter(algo_name=="LSVCMM", da>0) %>% pull(otu2) %>% unique()
# ------------------------------------------------------------------------------



# ==============================================================================
# Plot



cols = list(
  LSVCMM=list(display="LSVCMM", pos="left"),
  # LSVCM=list(display="LSVCM", pos="middle"),
  ALasso=list(display="ALasso", pos="middle"),
  # OLS=list(display="OLS", pos="middle"),
  SPFDA=list(display="SPFDA", pos="right")
)
estimates %<>% filter(algo_name %in% names(cols))
rows = list(
  KO=list(display="KO-WT", pos="top"),
  SCC=list(display="SCC-ED/CIS", pos="middle"),
  KO_SCC=list(display="Interaction", pos="bottom")
)


ncols = length(cols)
rel_widths=c(2, 1, 1, 2)
nrows = length(rows)
rel_heights = c(1, 1, 1)


gs = list()
for(variable in names(rows)){
  otu_da = estimates %>% filter(da>0, variable==!!variable) %>% pull(otu2) %>% unique()
  otu_da = sort(otu_da)
  ntaxas = length(otu_da)
  rows[[variable]]$ntaxas = ntaxas
  vdisplay = rows[[variable]]$display
  rpos = rows[[variable]]$pos
  for(method in names(cols)){
    display = cols[[method]]$display
    pos = cols[[method]]$pos
    dfc = estimates %>% filter(variable==!!variable, algo_name==!!method, otu2 %in% otu_da) %>% arrange(otu2, week)
    dfc %<>% mutate(estimate=pmax(-1.5, pmin(1.5, estimate)))
    g = ggplot() +
      theme_minimal() +
      geom_tile(
        data=dfc,
        mapping=aes(x=as.factor(week), y=otu2, fill=estimate)
      ) +
      xlab("Week") +
      ggtitle(display) +
      scale_fill_gradientn(
        colors=c("#00274C", "white", "#FFCB05"),
        values=c(0, 0.5, 1),
        breaks=c(-1, 0, 1),
        limits=c(-1.5, 1.5),
        aesthetics=c("fill"), name=vdisplay
      ) + theme(
        text=element_text(family="Helvetica"),
        plot.title=element_text(size=10),
        panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
      ) +
      geom_shadowtext(
        data=dfc %>% filter(da>0),
        mapping=aes(x=as.factor(week), y=otu2, label="*"),
        color="white", size=5, fontface="bold", hjust=0.5, vjust=0.75,
        bg.r=0.1
      ) + scale_y_discrete(limits=rev(otu_da)) +
      # add grid
      geom_hline(yintercept=seq(-0.5, ntaxas+0.5), color="grey80", linewidth=0.5)

    if(pos == "left"){
      g = g + ylab(tax_level)+ theme(axis.ticks.y=element_line())
    } else {
      g = g + theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank()
      )
    }
    if(rpos != "bottom"){
      g = g + xlab(NULL) + theme(
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()
      )
    }
    if(rpos != "top"){
      g = g + ggtitle("")
    }
    if(pos != "right")g = g + theme(legend.position="none")
    gs[[length(gs)+1]] = g
  }
}

heights = sapply(rows, function(x) x$ntaxas)
heights[1] = heights[1] + 2
heights[length(heights)] = heights[length(heights)] + 2

g = egg::ggarrange(plots=gs, ncol=ncols, nrow=nrows, heights=heights, widths=rep(5, ncols))

ggsave(
  paste0(DIR_FIGURES, "comparison_all.pdf"), g,
  width=6, height=8
  )
# ------------------------------------------------------------------------------





# ==============================================================================
# Plot sparsity and mean abundance
count_stats = data.frame(
  otu = rel_abundance %>% colnames(),
  mean = apply(rel_abundance, 2, mean),
  sparsity = apply(rel_abundance, 2, function(x) sum(x==0)/length(x)),
  prevalence = apply(rel_abundance, 2, function(x) sum(x>0)/length(x)),
  nz_mean = apply(rel_abundance, 2, function(x) mean(x[x>0]))
)

g = ggplot() +
  theme_minimal() +
  geom_point(data=count_stats, aes(x=1-sparsity, y=nz_mean)) +
  labs(x="Prevalence", y="Mean rel. abundance (non-zero)") +
  scale_y_log10(labels=scales::percent) +
  scale_x_continuous(labels=scales::percent, limits=c(0,1)) +
  geom_point(data=count_stats %>% filter(otu %in% sig_lsvcmm), aes(x=1-sparsity, y=nz_mean), color="darkred") +
  geom_label_repel(data=count_stats %>% filter(otu %in% sig_lsvcmm),
                  aes(x=1-sparsity, y=nz_mean, label=otu), size=2, color="darkred",
                  min.segment.length = 0, box.padding=0.5)
ggsave(
  paste0(DIR_FIGURES, "dmbt1_sparsity_mean.pdf"), g,
  width=6, height=4
)
# ------------------------------------------------------------------------------
