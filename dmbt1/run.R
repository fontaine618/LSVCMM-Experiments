library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)
library(phyloseq)


# ==============================================================================
# Prepare data
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
# Setup batchtools registry
DIR = paste0("./dmbt1/registry/")
if(dir.exists(DIR)) unlink(DIR, recursive=T)

registry = makeExperimentRegistry(
  file.dir=DIR,
  seed=1,
  packages=c("dplyr", "magrittr", "LSVCMM", "spfda", "tidyr")
)
registry$cluster.functions = makeClusterFunctionsMulticore(ncpus=8)
# registry = loadRegistry(DIR, writeable=T)
# ------------------------------------------------------------------------------




# ==============================================================================
# Setup problem
taxawise = function(data, job, otu){
  df = bind_cols(data.frame(data$meta), data$clr[, otu])
  colnames(df) = c("subject_id", "Week", "Type", "Gender", "Diagnosis", "SCC", "CLR")
  df %<>% mutate(
    KO=ifelse(Type=="KO", 1, 0),
    Sex=ifelse(Gender=="F", 1, 0),
    SCC=ifelse(Diagnosis=="SCC", 1, 0)
  ) %>% mutate(
    KO_SCC=KO*SCC
  )

  taxonomy = as.vector(data$tax[otu, ])
  while(is.na(tail(taxonomy, 1))) taxonomy = head(taxonomy, -1)
  best_taxonomy = tail(taxonomy, 1)
  instance = list(
    df=df,
    otu=otu,
    taxonomy=taxonomy,
    best_taxonomy=best_taxonomy
  )
  return(instance)
}

addProblem(
  name="taxawise",
  fun=taxawise,
  data=data
)
# ------------------------------------------------------------------------------





# ==============================================================================
# Setup algorithms
instance = taxawise(data, NULL, "Otu0354")

lsvcmm = function(
    data, job, instance,
    cross_sectional=F,
    independent=F
){
  method = "LSVCMM"
  if(independent){
    wc_args = list(name="independent")
    method = "LSVCM"
  }else{
    wc_args = list(name="compound_symmetry", estimate=T, ratio=1.)
  }
  if(cross_sectional){
    k_args = list(name="epa", scale=0.05, n_scale=1L)
    method = "ALasso"
  }else{
    k_args = list(name="gaussian", scale=0.2, n_scale=1L)
  }
  fit = LSVCMM::lsvcmm(
    data=instance$df %>% select(subject_id, Week, CLR, Sex, KO, SCC, KO_SCC),
    response="CLR",
    subject="subject_id",
    time="Week",
    vc_covariates=c("KO", "SCC", "KO_SCC", "Sex"),
    add_intercept=T,
    penalty=list(alpha=0.5, adaptive=0.5, penalize_intercept=F),
    kernel=k_args,
    working_covariance=wc_args,
    return_models=FALSE
  )
  i = which.min(fit$results$ebich)
  l = fit$results$penalty.lambda[i]
  h = fit$results$kernel.scale[i]
  boot = LSVCMM::lsvcmm.boot(
    data=instance$df %>% select(subject_id, Week, CLR, Sex, KO, SCC, KO_SCC),
    response="CLR",
    subject="subject_id",
    time="Week",
    vc_covariates=c("KO", "SCC", "KO_SCC", "Sex"),
    add_intercept=T,
    penalty=list(alpha=0.5, adaptive=0.5, penalize_intercept=F, lambda=l),
    kernel=k_args,
    working_covariance=wc_args,
    n_samples=1000
  )
  # grp_mean = matrix(c(
  #   1, 0, 0, 0, 0,
  #   1, 1, 0, 0, 0,
  #   1, 0, 1, 0, 0,
  #   1, 1, 1, 1, 0
  # ), 4, 5, T)
  # boot_grp_mean = LSVCMM:::transform_obj(boot, grp_mean)
  # grp_diff = matrix(c(
  #   0, 1, 0, 0, 0,
  #   0, 1, 0, 1, 0,
  #   0, 0, 1, 0, 0,
  #   0, 0, 1, 1, 0
  # ), 4, 5, T)
  # boot_grp_diff = LSVCMM:::transform_obj(boot, grp_diff)

  # out = list(
  #   estimates=list(
  #     KO=LSVCMM:::confidence_band(boot, var=2),
  #     SCC=LSVCMM:::confidence_band(boot, var=3),
  #     KO_SCC=LSVCMM:::confidence_band(boot, var=4)
  #   ),
  #   grp_mean=list(
  #     WT_HPCIS=LSVCMM:::confidence_band(boot_grp_mean, var=1),
  #     KO_HPCIS=LSVCMM:::confidence_band(boot_grp_mean, var=2),
  #     WT_SCC=LSVCMM:::confidence_band(boot_grp_mean, var=3),
  #     KOSCC=LSVCMM:::confidence_band(boot_grp_mean, var=4)
  #   ),
  #   grp_diff=list(
  #     HPCIS_KOmWT=LSVCMM:::confidence_band(boot_grp_diff, var=1),
  #     SCC_KOmWT=LSVCMM:::confidence_band(boot_grp_diff, var=2),
  #     WT_SCCmHPCIS=LSVCMM:::confidence_band(boot_grp_diff, var=3),
  #     KO_SCCmHPCIS=LSVCMM:::confidence_band(boot_grp_diff, var=4)
  #   )
  # )
  estimates = list(
    KO=LSVCMM:::confidence_band(boot, var=2),
    SCC=LSVCMM:::confidence_band(boot, var=3),
    KO_SCC=LSVCMM:::confidence_band(boot, var=4)
  ) %>% bind_rows(.id="Variable")
  estimates %<>% select(Variable, estimated_time, median, L, U, pval, excludes_zero) %>%
    rename(variable=Variable, week=estimated_time, estimate=median, lower=L,
           upper=U, pvalue=pval, da=excludes_zero)
  estimates %<>% mutate(algo_name=method, otu=instance$otu)

  return(estimates)
}

ols = function(
    data, job, instance
){
  # Long to Wide
  df = instance$df %>%
    tidyr::spread(Week, CLR, sep="_")
  # run OLS
  lms = list(
    Week_0 = lm(Week_0 ~ KO*SCC + Sex, data=df),
    Week_4 = lm(Week_4 ~ KO*SCC + Sex, data=df),
    Week_8 = lm(Week_8 ~ KO*SCC + Sex, data=df),
    Week_12 = lm(Week_12 ~ KO*SCC + Sex, data=df),
    Week_16 = lm(Week_16 ~ KO*SCC + Sex, data=df),
    Week_22 = lm(Week_22 ~ KO*SCC + Sex, data=df)
  )
  # Get coefficients
  ols_res = lapply(lms, function(x) {
    out = data.frame(summary(x)$coefficients)
    out = out %>% rownames_to_column("term")
  }) %>% bind_rows(.id="Week")
  colnames(ols_res) = c("Week", "term", "estimate", "std.error", "statistic", "p.value")
  # do bonferonni within term
  ols_res %<>% group_by(term) %>% mutate(p.value.adj=p.adjust(p.value, method="bonferroni"))
  ols_res %<>% arrange(term)
  # reformat
  ols_res %<>% mutate(Week=gsub("Week_", "", Week) %>% as.numeric())
  ols_res %<>% select(term, Week, estimate, p.value.adj) %>%
    rename(variable=term, week=Week, estimate=estimate, pvalue=p.value.adj) %>%
    mutate(da=ifelse(pvalue<0.05, 1, 0)) %>%
    filter(variable!="(Intercept)") %>%
    filter(variable!="Sex") %>%
    mutate(variable=ifelse(variable=="KO:SCC", "KO_SCC", variable))

  ols_res %<>% mutate(algo_name="OLS", otu=instance$otu)
  return(ols_res)
}


spfda = function(
    data, job, instance
){
  # Long to Wide
  df = instance$df %>% tidyr::spread(Week, CLR, sep="_")
  # Impute with fPCA
  Y = df %>% select(starts_with("Week")) %>% as.matrix()
  t0 = as.numeric(gsub("Week_", "", colnames(Y)))
  Yhat = refund::fpca.sc(Y=Y, argvals=t0, nbasis=5)$Yhat
  Yhat[!is.na(Y)] = Y[!is.na(Y)]
  # Covariate matrix
  X = df %>% select(KO, SCC, KO_SCC, Sex) %>% as.matrix()
  X = cbind(1, X)
  # Weight matrix
  nt = length(t0)
  W = spfda::spfda_weight(
    X=X,
    Y=Yhat,
    bandwidth=2,
    part=list(seq(nt))
  )
  # Grid of tuning parameters
  alphas = c(0.5)
  lambdas = 10^seq(log10(1), log10(100), length.out = 100)
  Ks = c(5, 6, 7, 8)
  all_params = expand.grid(lambdas, alphas, Ks)
  names(all_params) = c("lambda", "alpha", "K")
  # Fit models
  BICs = lapply(seq_len(nrow(all_params)), function(ii){
    param = all_params[ii,]
    lambda = param$lambda
    alpha = param$alpha
    K = param$K
    res = spfda::spfda(
      Y=Yhat,
      X=X,
      time=t0,
      lambda=lambda,
      alpha=alpha,
      nsp=K,
      ord=3,
      CI=F,
      W=W
    )
    BIC(res)
  })
  BICfd = cbind(
    all_params,
    BIC=unlist(BICs)
  )
  # Best BIC
  parms = BICfd[which.min(BICfd$BIC),  ]
  # get final model
  spfda = spfda::spfda(
    Y=Yhat,
    X=X,
    time=t0,
    lambda=parms$lambda,
    alpha=parms$alpha,
    nsp=parms$K,
    ord=3,
    CI=T,
    W=W
  )
  est = spfda$get_coef()
  se = spfda$get_se()
  # reformat
  out = list(
    KO=data.frame(week=t0, estimate=est["KO", ], se=se[2, ]),
    SCC=data.frame(week=t0, estimate=est["SCC", ], se=se[3, ]),
    KO_SCC=data.frame(week=t0, estimate=est["KO_SCC", ], se=se[4, ])
  ) %>% bind_rows(.id="variable")
  cval = qnorm(1-0.025/6)
  out %<>% mutate(lower=estimate-cval*se, upper=estimate+cval*se) %>%
    mutate(da=pmax((lower>0), (upper<0)))
  out %<>% mutate(algo_name="SPFDA", otu=instance$otu)
  return(out)
}


addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm
)
addAlgorithm(
  name="LSVCM",
  fun=lsvcmm
)
addAlgorithm(
  name="ALasso",
  fun=lsvcmm
)
addAlgorithm(
  name="SPFDA",
  fun=spfda
)
addAlgorithm(
  name="OLS",
  fun=ols
)
# ------------------------------------------------------------------------------





# ==============================================================================
# Experimental design
problems = list(
  taxawise=data.table(otu=otus)
)

algorithms = list(
  `LSVCMM`=data.table(cross_sectional=F, independent=F),
  `LSVCM`=data.table(cross_sectional=F, independent=T),
  `ALasso`=data.table(cross_sectional=T, independent=T),
  `SPFDA`=data.table(),
  `OLS`=data.table()
)

addExperiments(
  prob.designs=problems,
  algo.designs=algorithms,
  repls=1
)
# ------------------------------------------------------------------------------




# ==============================================================================
# Run
summarizeExperiments()
submitJobs(resources=list(walltime=10000))
getStatus()


# rerun some jobs
missing_ids = findNotDone()$job.id
submitJobs(resources=list(walltime=10000), ids=missing_ids)
# ------------------------------------------------------------------------------
