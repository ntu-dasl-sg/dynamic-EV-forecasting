# About the project

This repository contains R code for fitting and forecasting with the dynamic extreme value model in the context of volcanic eruptions. The model and analysis are documented in the paper, "A dynamic extreme value model with applications to volcanic eruption forecasting" by Michele Nguyen, Almut E. D. Veraart, Benoit Taisne, Tan Chiou Ting and David Lallemant.

# Usage

In the repository, we have the following R scripts:
- 1_traceenv.r: This computes envelopes for the seismic traces.
- 2_covar.r: This computes covariates based on the seismic traces.
- 3_env_cov.r: This computes covariates based on trace envelopes.
- 4_threshold.r: This focuses on threshold selection for the envelope index.
- 5_POT.r: This fits the dynamic extreme value model (POT = peak over threshold) on training data.
- 6_POT_test.r: Here, we evaluate the test performance of fitted model.
- 7_EVT_value.r: This script examines the value of extreme value theory using the model fits to the Jan 2010 event (Training Set 3) with different thresholds.
- GPD_regression.r: This includes functions for conducting step-wise variable selection based on the Akaike Information Criterion (AIC) for the Generalised Pareto Distribution (GPD) regression. 

To follow the analysis in the paper, the R scripts should be run in order, excluding GPD_regression.r (this is called into 5_POT.r and 7.EVT_value.r to fit the GPD regression when required).
