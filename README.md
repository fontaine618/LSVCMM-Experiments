# LSVCMM-Experiments
Experiments accompnaying the LSVCMM paper



## Notes

- Longitudinal seems to uniformly improve on Independent, especially in estimation error
- Cross-sectional very rarely selects anything
- Longitudinal and Independent select quite a lot more than SPFDA; this results in larger error
  for null time points and larger FDR, but larger power and smaller error for non-nulls
- LSVCMM seems quite robust to the misspecification I considered
- Most trends seem natural with the x-axes (at least it seems possible to explain them all),
  except for drop in accuracy for prop observed=1 ... 
  Let's look without selection to investiggate.

- For the comparison with fixed tuning parms, I don't think lambda are really comparable
- I did the same by fixing the scale, but not the kernel only

### 12/09

- drop CS in 4th experiment
- estimation error, acc, tpr power
- time points, function
- compare with fixed tuning values

- tvreg
- new setting for sanity cehck with other function
- look into why acc drop for 2nd experiment at last x value

### 19/09

- fix kernel scale for experiments


### 26/09

- understand last time point ??
- more timepoints
- different curves

### 02/10

- Issue with lasso and translation
- you would think the adaptive penalty helps, but not really since the unpenalized estimate is also shifted

### 25/10
- Two step improves a bit, not perfect though
- Difference in estimated cov parameter: full optimization finds the true value, 2S finds something smaller

todo:
[ ] change missingness to block sqrt*qqrt
[ ] understand SPFDA issue? I think it is imputation
[ ] redo experiments + add more
[ ] check boundary weighing / implement data-adaptive

[ ] Looks like I was not setting tolerence small enough, still some movement between 1e-6 and 1e-8
[ ] Maybe need to quicken convergence


Profiling:
- small improvement possible with logdet_precision if defined per working covariance (3-4%)
- biggest improvement would be to reduce the number of total iterations
- mean update is 2x more expensive than cov updates
- if there is a way to speed up linear predictor, that'd be great
