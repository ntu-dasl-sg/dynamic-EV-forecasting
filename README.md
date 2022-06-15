# About the project

This repository contains R code for fitting and forecasting with the dynamic extreme value model in the context of volcanic eruptions. The model and analysis are documented in the paper, "A dynamic extreme value model with applications to volcanic eruption forecasting" by Michele Nguyen, Almut E. D. Veraart, Benoit Taisne, Tan Chiou Ting and David Lallemant.

# Code usage

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

# Data

The data used for the analysis was collected by the Institut de Physique du Globe de Paris, Observatoire Volcanologique du Piton de la Fournaise (IPGP/OVPF) and the Laboratoire de Gèophysique Interne et Tectonophysique (LGIT) within the framework of [ANR\_08\_RISK\_011/UnderVolc](https://www.fdsn.org/networks/detail/YA_2009/) project. A publicly available dataselect service can be used to obtain the seismic data in miniseed files. For Training Events 1 and 2, the Test Event and Non-events, the queries used are given below:

<details>
<summary>Training Event 1: 11/05/2009 E-L event on Flank S SE</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2009-11-04T00:00:00&endtime=2009-11-07T00:00:00&quality=B
</details>

<details>
<summary>Training Event 2: 12/14/2009 E-L event on Flank S/SW</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2009-12-13T00:00:00&endtime=2009-12-15T00:00:00&quality=B
</details>

<details>
<summary>Test Event: 10/14/2010 E-L event on Flank S</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2010-10-13T00:00:00&endtime=2010-10-15T00:00:00&quality=B
</details>

<details>
<summary>Test Non-event 1: 11/30/2009</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2009-11-30T00:00:00&endtime=2009-12-01T00:00:00&quality=B
</details>

<details>
<summary>Test Non-event 2: 12/22/2009</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2009-12-22T00:00:00&endtime=2009-12-24T00:00:00&quality=B
</details>

<details>
<summary>Test Non-event set 3: 05/08/2010</summary>
http://ws.resif.fr/fdsnws/dataselect/1/query?net=YA&sta=UV05&loc=00&cha=HHZ&starttime=2010-05-08T00:00:00&endtime=2010-05-10T00:00:00&quality=B
</details>

The data for Training Event 3 was available via text files but can also be obtained via the dataselect service.
