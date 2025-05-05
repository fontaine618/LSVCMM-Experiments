splinectomer_wrapper = function(
    data, job, instance
){
  time0 = proc.time()

  df = instance$data
  df$group = factor(df$group)
  df$subject_id = factor(df$subject_id)

  # prepare intervals
  t0 = instance$estimated_time

  # run
  fit = splinectomeR::sliding_spliner(
    data=df,
    xvar="time",
    yvar="response",
    category="group",
    cases="subject_id",
    ints=length(t0)
  )

  # extract results

  pvals = data.frame(time=t0, pval=fit$pval_table$p_value)
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
    fit=fit,
    time=proc.time()-time0
  )

}







