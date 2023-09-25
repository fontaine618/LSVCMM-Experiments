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


