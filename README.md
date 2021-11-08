# Read-me

This repository contains the code used in statistical analysis presented in the manuscript 'Predicting the impacts of land management for sustainable development on depression risk in a Ugandan case study'.

This includes three scripts: 
1) ‘Data_exploration.R’ - Exploring the data, constructing the latent variables, etc. To run this script, ensure that ‘DF_anly.rds’ is in the working directory. These data can be retrieved from DOI: 10.6084/m9.figshare.16955221. This script also produces the imputed datasets used in following script. 
2) ‘Data_analysis.R’ - The Bayesian structural equation model, model diagnostics, and supplementary analysis. ‘Data_exploration.R’ has to be run to produce ‘DF_analy.list.rds’ (the multiple imputed dataset) in order to run this script. 
3) ‘function.R’ - Some functions used in the two scripts. 
