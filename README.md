# Bayesian VAR order determination

This repository contains the code and datasets associated with the manuscript 'Bayesian inference on the order of stationary vector autoregressions'.  

## Files

The following files are included in the repository:
* multiplicative_gamma.stan - this is the associated Stan program used to fit a zero-mean vector autoregressive model with a maximum order, $p_{\max}$, using a multiplicative gamma process prior for the unconstrained partial autocorrelation matrices and an inverse Wishart prior for the error variance. We used Stan version 2.29.2 to carry out the analyses.
* functions.R - this is an R file containing functions required to simulate a new dataset. 
* run.R - this is an R file which provides illustrative R code to run the Stan program. The code allows the user to either generate a new dataset or use one of the datasets used in the simulation experiments in Section 5 of the manuscript.

## Data

In addition to the files, the data directory contains the datasets used in the simulation experiments. An example of how to run the simulation experiments for these datasets is included in the `run.R` file. Unfortunately we do not have permission to release the EEG data used in the application.
