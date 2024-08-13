library(batchtools)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(shadowtext)


# ==============================================================================
# Setup batchtools registry
DIR = paste0("./dmbt1/registry_defense/")
DIR_FIGURES = paste0("./dmbt1/figures_defense/")
registry = loadRegistry(file.dir=DIR)
t0 = c(0, 4, 8, 12, 16, 22)
tax_level = "OTU"
# ------------------------------------------------------------------------------




# ==============================================================================
# Load results
estimate = function(result) result
estimates = reduceResultsList(fun = estimate) %>% bind_rows(.id="job.id")
estimates %<>% mutate(job.id = as.numeric(job.id)) %>% rename(otu2=otu)
parameters = getJobPars() %>% unwrap()
parameters %<>% mutate(job.id = as.numeric(job.id))
estimates %<>% left_join(parameters, by="job.id")
# # find missing jobs
# edf = estimates %>% select(job.id, otu2, algo_name) %>%
#   group_by(job.id, otu2, algo_name) %>%
#   summarise(n=n()) %>%
#   select(-n) %>% rename(otu=otu2, algorithm=algo_name, ejob.id=job.id)
# pdf = parameters %>% select(job.id, otu, algorithm) %>%
#   group_by(job.id, otu, algorithm) %>%
#   summarise(n=n()) %>%
#   select(-n) %>% rename(pjob.id=job.id)
# missing = full_join(edf, pdf, by=c("otu", "algorithm")) %>% filter(is.na(ejob.id))
# missing_ids = missing$pjob.id

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
# ------------------------------------------------------------------------------



# ==============================================================================
# Plot



cols = list(
  LSVCMM=list(display="LSVCMM", pos="left"),
  # LSVCM=list(display="LSVCM", pos="middle"),
  # ALasso=list(display="ALasso", pos="middle"),
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
