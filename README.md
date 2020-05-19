# P2C2M_PSMC
Implementation of P2C2M with the Pairwise Sequentially Markovian Coalescent (PSMC)

This repository includes R scripts for the Posterior Predictive Checks of Coalescent Models (P2C2M) R package. This version conducts posterior predictive checks of output from Li and Durbin's PSMC program that estimates historical changes in effective population size. We found that an identifiability issue prevents the pipeline from being able to detect model violations, but the scripts may still be of use as utility scripts for those connducting PSMC analyses. As a whole, the only dependencies are the parallel R package and Hudson's ms simulator, but their requirement is script specific.
