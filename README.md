# Reproducible results for "Estimation of the error structure in multivariate response linear regression models"

This codebase reproduces the simulation study section in the paper "Estimation of the error structure in multivariate response linear regression models" by Liu and Yu (2025), available open-access at [https://doi.org/10.1002/wics.70021](https://doi.org/10.1002/wics.70021).

## Directions

1. Run `generate_sim_data.R` to generate the simulation data.
2. The `run_sim_*.R` files run the simulation scenarios.
3. The tables and figures in the paper can be reproduced by running the `make_plot_*.R` and `make_tables_*.R` files.

The results related to conditional Gaussian graphical models depend on the `ANTAC` and `CAPME` R packages, which are not publicly available. Please email me if you wish to run these simulations.
When these packages are available, the files `generate_sim_ggm_data.R` and `run_sim_ggm_*.R` can be used to reproduce these results.

`methods.R` and `methods_ggm.R` contain the functions used in the simulation study.
`cov_srrr.R` is our custom implementation of the method of Chen & Huang (2016) https://doi.org/10.1007/s11222-014-9517-6.
`performance.R` contains the functions used to evaluate the performance of the methods.

A Dockerfile is also provided to run the code in a containerized environment.

## Citation

Liu, R. and Yu, G. (2025), Estimation of the Error Structure in Multivariate Response Linear Regression Models. WIREs Comput Stat, 17: e70021. https://doi.org/10.1002/wics.70021
