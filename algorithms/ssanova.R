ssanova_wrapper = function(
    data, job, instance,
    n.perm=1000
){
  time0 = proc.time()

  df = instance$data
  df$group = factor(df$group)
  df$subject_id = factor(df$subject_id)

  se_obj = SummarizedExperiment::SummarizedExperiment(
    assays=list(matrix(df$response, nrow=1)),
    colData=data.frame(
      Group=df$group,
      Subject=df$subject_id,
      Time=df$time
    )
  )

  # prepare intervals
  # if d is the interval length, then we set the intervals to
  # 0, d/2, d + d/2, 2d + d/2, ..., (n-1)d + d/2, 1
  # so that the middle point of the interval corresponds to the estimated time
  t0 = instance$estimated_time
  gap = mean(diff(t0))
  points = t0[-1] - gap/2
  points = c(t0[1], points, t0[length(t0)])

  # run
  fit = OmicsLonDA::omicslonda(
    se_object=se_obj,
    n.perm=n.perm,
    fit.method="ssgaussian",
    points=points,
    text="Group"
  )

  # extract results

  pvals = data.frame(time=t0, pval=fit$details$intervals.pvalue)
  pvals$pval.adj = p.adjust(pvals$pval, method="hommel")
  pvals$decision = pvals$pval.adj < 0.05

  classification = data.frame(
    time=instance$estimated_time,
    classification=ifelse(
      instance$true_values$b1 == 0,
      ifelse(pvals$decision, "FP", "TN"), # Negative: Detection, No Detection
      ifelse(pvals$decision, "TP", "FN") # Positive: Detection, No Detection
    )
  )

  list(
    decision=pvals,
    classification=classification,
    time=proc.time()-time0
  )

}




