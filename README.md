# Reproduction of Results

The repository is composed of separate R files which contain the code for constructing examples (`covs.r`, `covsReg.r`), implementing the procedure (`precn.r`, `CovEstmtr.r`), and evaluating the results (`metric.r`).

## Numerical Simulations

For reproducing the results in **Section 8: Numerical Simulations**, run the files `Simulations - ***.r` in the respective directories in `Simulations`. They rely on the R package `doParallel` for parallel processing. We assume ~25 cores, but this can be modified in the function `AUC`.

The reason for having separate directories for every type of simulation is that for different cases the procedures for the construction of examples, their simulation, covariance estimation and parameter tuning had to be implemented slightly differently.

## Real Data Examples

For reproducing the results in **Section 9: Real Data Examples**, run the files `RealData - Section 9.1 - Infrared Absorbtion Spectroscopy.r` and `RealData - Section 9.2 - Cumulative Log Returns for Pfizer Ltd.r` in the directory `Real Data Examples`.

The data required for reproducing the results in **Section 9: Real Data Examples** of the article:
1. The file `MIR_Fruit_purees.csv` used in `RealData - Section 9.1 - Infrared Absorbtion Spectroscopy.r` can be found at [Mendeley Data](https://data.mendeley.com/datasets/frrv2yd9rg/1).
2. The file `PFIZER__EQ__NSE__NSE__MINUTE.csv` used in `RealData - Section 9.2 - Cumulative Log Returns for Pfizer Ltd.r` can be found on [Kaggle](https://www.kaggle.com/datasets/hk7797/stock-market-india).