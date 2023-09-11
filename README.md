# LSVCMM-Experiments
Experiments accompnying the LSVCMM paper



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
