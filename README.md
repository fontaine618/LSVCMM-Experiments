# LSVCMM-Experiments
Experiments accompanying the LSVCMM paper

## Simulation experiments

### Missing data in Regular Design (3.1)

Results are produced in `sim_block.r` into the figure `sim_block.pdf`.

The four experiments are in the corresponding folders:

- `experiment_snr`: Varying noise variance
- `experiment_missing` : Varying proportion of missing data
- `experiment_re_ratio`: Varying random effect variance
- `experiment_n`: Varying sample size

Additional metrics are produced in `sim_block_supp.r` into the figure `sim_block_supp.pdf`.

### Irregular Design (3.2)

Results are produced in `sim_sparse.r` into the figure `sim_sparse.pdf`.

The four experiments are in the corresponding folders:

- `experiment_snr100`: Varying noise variance
- `experiment_missing100` : Varying proportion of missing data
- `experiment_re_ratio100`: Varying random effect variance
- `experiment_n100`: Varying sample size

Additional metrics are produced in `sim_sparse_supp.r` into the figure `sim_sparse_supp.pdf`.

### Correlation Structure (S2.2)

This experiment is contained in the `experiment_ar` folder.

### Significance experiment (S2.3)

This experiment is contained in the `experiment_pvalues` folder.

### Sparse group Lasso weights (S2.4)

This experiment is contained in the `experiment_sgl` folder.

### Inference experiment (S2.5)

Results are produced in `sim_sig.r` into the figure `sim_sig.pdf`.

The corresponding experiment are in the following folders:

- `experiment_sig10missing`
- `experiment_sig10n`
- `experiment_sig10re`

### Simulated compositional data (S2.6)


## Real data experiments

The `dmbt1` folder contains the code to reproduce the results on the DMBT1 dataset. 
