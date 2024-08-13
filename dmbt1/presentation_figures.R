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
for(otu in c(
  "Otu0004",
  "Otu0007",
  "Otu0018",
  "Otu0025",
  "Otu0061",
  "Otu0078",
  "Otu0133",
  "Otu0141",
  "Otu0337",
  "Otu0458",
  "Otu0538",
  "Otu0012",
  "Otu0038",
  "Otu0058",
  "Otu0083",
  "Otu0097",
  "Otu0139",
  "Otu0318",
  "Otu0025",
  "Otu0054",
  "Otu0083",
  "Otu0107",
  "Otu0218",
  "Otu0318",
  "Otu0354",
  "Otu0648"
)){
  cols = list(
    KO=list(display="KO-WT", pos="left"),
    SCC=list(display="SCC-ED/CIS", pos="middle"),
    KO_SCC=list(display="Interaction", pos="right")
  )
  rows = list(
    ALasso=list(display="ALasso", pos="middle"),
    # LSVCM=list(display="LSVCM", pos="middle"),
    LSVCMM=list(display="LSVCMM", pos="left"),
    SPFDA=list(display="SPFDA", pos="right")
  )
  ncols = length(cols)
  rel_widths=c(1.1, 1, 1)
  nrows = length(rows)
  rel_heights = c(1.1, 1, 1.1)

  yrange = estimates %>% filter(otu2==!!otu) %>% select(lower, upper) %>% unlist() %>% range()

  gs = list()
  for(row in seq_along(rows)){
    algo_name = names(rows)[row]
    adisplay = rows[[algo_name]]$display
    for(col in seq_along(cols)){
      variable = names(cols)[col]
      vdisplay = cols[[variable]]$display
      prow = rows[[algo_name]]$pos
      pcol = cols[[variable]]$pos
      dfc = estimates %>% filter(algo_name==!!algo_name, variable==!!variable, otu2==!!otu) %>% arrange(week)

      g = ggplot() +
        theme_minimal() +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        geom_ribbon(
          data=dfc,
          mapping=aes(x=week, ymin=lower, ymax=upper, alpha=0.2),
          alpha=0.5
        ) +
        geom_line(
          data=dfc,
          mapping=aes(x=week, y=estimate)
        ) +
        labs(color="Algorithm", fill="Algorithm") +
        theme(
          text=element_text(family="Helvetica"),
          plot.title=element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
        ) +
        scale_x_continuous(breaks=t0, minor_breaks=NULL) +
        scale_y_continuous(limits=yrange)

      if(row!=nrows){
        g = g + theme(
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
          )
      }else{
        g = g + labs(x="Week")
      }
      if(row==1){
        g = g + labs(title=vdisplay)
      }
      if(col!=1){
        g = g + theme(
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
          )
      }else{
        g = g + labs(y=adisplay)
      }
      gs[[length(gs)+1]] = g
    }
  }



  g = egg::ggarrange(
    plots=gs, ncol=ncols, nrows=nrows,
    align="v", axis="tb", byrow=T,
    widths=rel_widths, heights=rel_heights)


  ggsave(
    paste0(DIR_FIGURES, "comparison_", otu, ".png"), g,
    width=7, height=7
    )

}
# ------------------------------------------------------------------------------
