####################################
######### Primary analysis #########
####################################

# The following explores the data used in the 'Predicting the impacts of land management for sustainable development on depression risk in a Ugandan case study'  

### This script contains the following steps ###
# 1) Set the environment and load data 
# 2) Perform the model diagnostics from the WAMBS check list (Depaoli & van de Schoot 2017):
# 2.1) 'Do you understand the priors?'
# 2.2) 'Does the trace-plot exhibit convergence?' (After running the model.)
# 2.3) 'Does convergence remain after doubling the number of iterations?'
# 2.4) 'Does the histogram have enough information?'
# 2.5) 'Do the chains exhibit a strong degree of autocorrelation?'
# 2.6) 'Does the posterior distribution make substantive sense?'
# 2.7) 'Do different specifications of the multivariate variance priors influence the results?'
# 2.8) 'Is there a notable effect of the prior when compared with noninformative priors?'
# 2.9) 'Are the results stable from a sensitivity analysis?'
# 2.10) 'Is the Bayesian way of interpreting and reporting model results used?' 
# 3) Run the model on the ten imputed datasets
# 4) Extract and pool the model results 
# 5) Additional analysis associating "thinking too much" and PHQ-8 scores 

### Note: Data_exploration.R and Data_analysis.R use anonymised data (DF_anly.rds) and can be reproduced.  
### Users should ensure the DF_anly.rds data (from DOI: 10.6084/m9.figshare.16955221) is in the working directory. 
### Data_exploration.R need to be run before this one, since creates the DF_analy.list.rds imputed dataset used in the analysis. 

######### 1) Set the environment and load data #########

### Load packages ###
library(plyr)
library(tidyverse)
library(ggplot2)
library(mice)
library(lavaan)
library(semTools)
library(blavaan)
library(semPlot)
library(coda)
library(brms)

### Call the 'functions' script ###
source("functions.R")

### Seed ###
set.seed(4343)

### Load the data used in the statistical analysis ###
DF_analy.list <- readRDS("DF_analy.list.rds")

######  2.1) 'Do you understand the priors?' ###### 
# See the main text and SX for further details. The following is simply used to generate the plots used in SX. 

# Function for creating plots for SX. 
prior.plot <- function(type="normal", mean=0, variance=1, shape.a=1, shape.b=1, sec.min=-6, sec.max=6, step=.01, label=label) {
    x <- seq(sec.min, sec.max, by = step)
  
  # For a normally distributed prior 
  if (type == "normal") {
    prior.d <- dnorm(x,mean = mean, sd = sqrt(variance))
  }
  
  # For a beta distributed prior 
  if (type == "beta") {
    prior.d <- dbeta(x, shape1  = shape.a, shape2 = shape.b)
  }
  
  # Plot 
  df <- data.frame(x = x, prior.d = prior.d)  
  ggplot(data=df, aes(x=x, y=prior.d, group=1)) +
      geom_line(size = .5) +
  xlab("X") +
  ylab("Prob. den.") + 
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.1))
  
  # Save last plot 
  ggsave(filename =  paste0("C:/Users/wolf5246/Dropbox/Oxford/PhD/Chapter 3/Manuscript/Preperation/SI_prep/Prior_plots/", label, ".png"), plot = last_plot(), width = 40, height = 40, units = "mm", dpi = 400,)
  
}

# Depression ~ Food insecurity
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="1. Depression ~ Food insecurity")

# Depression ~ Economic poverty
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="2. Depression ~ Economic poverty")

# Food insecurity ~ Forest use
prior.plot(type = "beta", shape.a = 12, shape.b = 11, sec.min=0, sec.max=1, 
           label="3. Food insecurity ~ Forest use")

# Food insecurity	~	Farm size	
prior.plot(type = "normal", mean = -0.25, variance = 4, 
           label="4. Food insecurity ~ Farm size")

# Food insecurity	~	Economic poverty
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="5. Food insecurity ~ Economic poverty")

# Food insecurity	~	Distance for forest reserve	
prior.plot(type = "normal", mean = -0.25, variance = 4, 
           label="6. Food insecurity ~ Distance for forest reserve")

# Economic poverty ~ Forest use
prior.plot(type = "beta", shape.a = 12, shape.b = 11, sec.min=0, sec.max=1, 
           label="7. Economic poverty ~ Forest use")

# Economic poverty ~ Farm size
prior.plot(type = "beta", shape.a = 11, shape.b = 12, sec.min=0, sec.max=1, 
           label="8. Economic poverty ~ Farm size")

# Depression ~ Age
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="9. Depression ~ Age")

# Depression ~ Gender 
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="10. Depression ~ Gender")

# Depression ~ Education 
prior.plot(type = "normal", mean = -0.25, variance = 4, 
           label="11. Depression ~ Education")

# Depression ~ Social support 
prior.plot(type = "normal", mean = -0.25, variance = 4, 
           label="12. Depression ~ Social support")

# Depression ~ Divorced or widowed   
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="13. Depression ~ Divorced or widowed")

# Depression ~ Never married  
prior.plot(type = "normal", mean = 0, variance = 9,  
           label="14. Depression ~ Never married")

# Depression ~ General health 
prior.plot(type = "beta", shape.a = 11, shape.b = 12, sec.min=0, sec.max=1, 
           label="15. Depression ~ General health")

# Depression ~ Alcohol consumption 
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="16. Depression ~ Alcohol consumption")

# Depression ~ Smoking
prior.plot(type = "normal", mean = 0.25, variance = 4, 
           label="17. Depression ~ Smoking")

# Depression ~ Community (all)
prior.plot(type = "normal", mean = 0, variance = 9, sec.min=-9, sec.max=9,  
           label="18. Depression ~ Community (all)")

######  2.2) 'Does the trace-plot exhibit convergence?' ###### 
### Define the model and set the priors ###
# Stan uses standard deviations as the priors (rather than variance, like in Jags). See here for useful details: http://ecmerkle.github.io/blavaan/articles/prior.html 
# This means we have to take the square root of the variance to get the standard deviation, which is supplied as the prior.
# Diagnostics are performed using the first of the ten imputed datasets.

# Chains
chains = 4 

# Burn-in iterations 
burnin <- 4000 # 4000

# Post-burn-in iterations
sample <- 4000 # 4000

# Set seed
seed = 4343

# Check - standard deviation = square root of variance
sqrt(4) # = Variance = 4
sqrt(9) # = Variance = 9

# The model and priors 
model_main_1 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ prior("normal(0.25, 2)")*FIES_est + prior("normal(0.25, 2)")*Economic_poverty  
    FIES_est ~ prior("normal(-0.25, 2)")*Land_est + prior("normal(0.25, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(-0.25, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ prior("beta(11, 12)")*Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ prior("beta(11, 12)")*Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ prior("normal(0.25, 2)")*DOBB + 
    
    # Sex
    prior("normal(0.25, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(-0.25, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(-0.25, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0.25, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0.25, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0.25, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ prior("normal(0, 9)")*LOCB_EWAF + prior("normal(0, 9)")*LOCB_KADONE + prior("normal(0, 9)")*LOCB_KADTWO + 
    prior("normal(0, 9)")*LOCB_KARO + prior("normal(0, 9)")*LOCB_KYEM + prior("normal(0, 9)")*LOCB_MARA + 
    prior("normal(0, 9)")*LOCB_NONE + prior("normal(0, 9)")*LOCB_NTWO + prior("normal(0, 9)")*LOCB_NYAB + prior("normal(0, 9)")*LOCB_NYAK
    '

# Bayesian SEM
fit_main_1 <- bsem(model_main_1, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_1, "fit_main_1.rds")

# Plot estimates 
plot(fit_main_1,  plot.type = "areas", pars = 1:27, prob = 0.90,
     prob_outer = 0.95)

# Load the model (if needed)
fit_main_1 <- readRDS("fit_main_1.rds")

# Traceplots for key parameters - first ten 
key_params <- 27
plot(fit_main_1, pars = 1:10, plot.type = "trace")

# Traceplots for key parameters - secound ten 
plot(fit_main_1, pars = 11:20, plot.type = "trace")

# Traceplots for key parameters - remaining 
plot(fit_main_1, pars = 21:key_params, plot.type = "trace")

# # Geweke diagnostic
# fit_main_1_mcmc.list <- blavInspect(fit_main_1, what = "mcmc")
# geweke.plot(fit_main_1_mcmc.list) # We expect about 5% to be out more than +/- 1.96. 
# 
# # Gelman and Rubin diagnostic (rule of thumb is that everything below 1.1 is OK): https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
# gelman.diag(fit_main_1_mcmc.list) 
# gelman.plot(fit_main_1_mcmc.list)


######  2.3) 'Does convergence remain after doubling the number of iterations?' ###### 

# Double the burn-in 
burnin.d <- burnin*2

# Double the post burn-in
sample.d <- sample*2

# Bayesian SEM  
fit_main_2 <- bsem(model_main_1, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin.d, sample = sample.d, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_2, "fit_main_2.rds")

# # Load the model (if needed)
# fit_main_2 <- readRDS("fit_main_2.rds")

# Traceplots for key parameters - first ten 
plot(fit_main_2, pars = 1:10, plot.type = "trace")

# Traceplots for key parameters - secound ten 
plot(fit_main_2, pars = 11:20, plot.type = "trace")

# Traceplots for key parameters - remaining 
plot(fit_main_2, pars = 21:key_params, plot.type = "trace")

# # Geweke diagnostic
# fit_main_2_mcmc.list <- blavInspect(fit_main_2, what = "mcmc")
# geweke.plot(fit_main_2_mcmc.list) # We expect about 5% to be out more than +/- 1.96. 
# 
# # Gelman and Rubin diagnostic (rule of thumb is that everything below 1.1 is OK): https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
# gelman.diag(fit_main_2_mcmc.list) 
# gelman.plot(fit_main_2_mcmc.list)


######  2.4) 'Does the histogram have enough information?' ###### 
plot(fit_main_1 , pars = 1:key_params, plot.type = "hist")


######  2.5) 'Do the chains exhibit a strong degree of autocorrelation?'###### 
# Params 1 - 4
plot(fit_main_1, pars = 1:4, plot.type = "acf")

# Params 5 - 8
plot(fit_main_1, pars = 5:8, plot.type = "acf")

# Params 9 - 12
plot(fit_main_1, pars = 9:12, plot.type = "acf")

# Params 13 - 16
plot(fit_main_1, pars = 14:16, plot.type = "acf")

# Params 17 - 20
plot(fit_main_1, pars = 17:20, plot.type = "acf")

# Params 21 - 24
plot(fit_main_1, pars = 21:24, plot.type = "acf")

# Params 25 - 27
plot(fit_main_1, pars = 22:key_params, plot.type = "acf")

######  2.6) 'Does the posterior distribution make substantive sense?' ###### 
# Combing the results of the three chains together 
fit_1_MCMCbinded <- as.matrix(fit_main_1_mcmc.list)

# Examine the posterior for each key variable  
par(mfrow = c(4,7))
for (i in seq_along(1:key_params)) {
  plot(density(fit_1_MCMCbinded[,i]))
}
dev.off()

######  2.7) 'Do different specifications of the multivariate variance priors influence the results?' ###### 

# Examine the observed variable precision parameter (observed, because we are using plausible valued) - theta (which appears to be called 'itheta' in the manual?) - which has a default prior of gamma(1, 0.5).
dpriors(target = "stan")

# Respecify the model with a different variance associated with PHQ8_est: prior("gamma(.5, .5)")
model_main_3 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ prior("normal(0.25, 2)")*FIES_est + prior("normal(0.25, 2)")*Economic_poverty  
    FIES_est ~ prior("normal(-0.25, 2)")*Land_est + prior("normal(0.25, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(-0.25, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ prior("beta(11, 12)")*Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ prior("beta(11, 12)")*Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ prior("normal(0.25, 2)")*DOBB + 
    
    # Sex
    prior("normal(0.25, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(-0.25, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(-0.25, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0.25, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0.25, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0.25, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ prior("normal(0, 9)")*LOCB_EWAF + prior("normal(0, 9)")*LOCB_KADONE + prior("normal(0, 9)")*LOCB_KADTWO + 
    prior("normal(0, 9)")*LOCB_KARO + prior("normal(0, 9)")*LOCB_KYEM + prior("normal(0, 9)")*LOCB_MARA + 
    prior("normal(0, 9)")*LOCB_NONE + prior("normal(0, 9)")*LOCB_NTWO + prior("normal(0, 9)")*LOCB_NYAB + prior("normal(0, 9)")*LOCB_NYAK
    
    # Using a alternative to the default 
    PHQ8_est ~~ gamma(1,.05)[sd]*PHQ8_est
'

# Bayesian SEM  
fit_main_3 <- bsem(model_main_3, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_3, "fit_main_3.rds")

# # Load the model (if needed)
# fit_main_3 <- readRDS("fit_main_3.rds")

# Calculate the bias associated with using this different specification
fit_main_1_sum <- data.frame(summary(fit_main_1))[,1:2]
fit_main_2_sum <- data.frame(summary(fit_main_3))[,1:2]
round(100*(as.numeric(fit_main_1_sum$Estimate)-as.numeric(fit_main_2_sum$Estimate))/as.numeric(fit_main_2_sum$Estimate),2) # Fine, if it is only small 

######  2.8) 'Is there a notable effect of the prior when compared with noninformative priors?' ###### 
# blavaan has default very weakly informative priors, so the analysis is repeated without specified priors

# Respecify the model with default weakly informative priors
model_main_4 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ FIES_est + Economic_poverty  
    FIES_est ~ Land_est + Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ Forest_est.1 
    Economic_poverty ~~ Forest_est.1
    Economic_poverty ~~ Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ DOBB + 
    
    # Sex
    SEX +

    # Education (treated as numeric)
    EDUCAT + 
    
    # Social support 
    Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    MSTATUS_Div.wid + MSTATUS_Single +
    
    # Alchohol consumption  
    Alcohol + 
    
    # Smoking 
    SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ LOCB_EWAF + LOCB_KADONE + LOCB_KADTWO + LOCB_KARO + LOCB_KYEM + LOCB_MARA + LOCB_NONE + LOCB_NTWO + LOCB_NYAB + LOCB_NYAK
    '

# Bayesian SEM  
fit_main_4 <- bsem(model_main_4, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_4, "fit_main_4.rds")

# # Load the model (if needed)
# fit_main_4 <- readRDS("fit_main_4.rds")

# Calculate the bias associated with using weakly informative priors
fit_main_1_sum <- data.frame(summary(fit_main_1))[,1:2]
fit_main_4_sum <- data.frame(summary(fit_main_4))[,1:2]
round(100*(as.numeric(fit_main_1_sum$Estimate)-as.numeric(fit_main_4_sum$Estimate))/as.numeric(fit_main_4_sum$Estimate),2) # Fine, if it is only small 

######  2.9) 'Are the results stable from a sensitivity analysis?' ###### 
# Shifting the hyperparameters 'up' (specifically the mean and the shape of the beta distribution) 
model_main_5 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ prior("normal(0.5, 2)")*FIES_est + prior("normal(0.5, 2)")*Economic_poverty  
    FIES_est ~ prior("normal(-0, 2)")*Land_est + prior("normal(0.5, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(0, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ prior("beta(4, 2)")*Forest_est.1
    Economic_poverty ~~ prior("beta(4, 2)")*Forest_est.1
    Economic_poverty ~~ prior("beta(3, 3)")*Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ prior("beta(3, 3)")*Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ prior("normal(0.5, 2)")*DOBB + 
    
    # Sex
    prior("normal(0.5, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(0, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(0, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0.5, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0.5, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0.5, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ prior("normal(0.25, 9)")*LOCB_EWAF + prior("normal(0.25, 9)")*LOCB_KADONE + prior("normal(0.25, 9)")*LOCB_KADTWO + 
    prior("normal(0.25, 9)")*LOCB_KARO + prior("normal(0.25, 9)")*LOCB_KYEM + prior("normal(0.25, 9)")*LOCB_MARA + 
    prior("normal(0.25, 9)")*LOCB_NONE + prior("normal(0.25, 9)")*LOCB_NTWO + prior("normal(0.25, 9)")*LOCB_NYAB + prior("normal(0.25, 9)")*LOCB_NYAK
    '
# Bayesian SEM  
fit_main_5 <- bsem(model_main_5, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_5, "fit_main_5.rds")

# # Load the model (if needed)
# fit_main_5 <- readRDS("fit_main_5.rds")

# Calculate the bias associated with using 'up' shifted hyperparameters
fit_main_1_sum <- data.frame(summary(fit_main_1))[,1:2]
fit_main_5_sum <- data.frame(summary(fit_main_5))[,1:2]
round(100*(as.numeric(fit_main_1_sum$Estimate)-as.numeric(fit_main_5_sum$Estimate))/as.numeric(fit_main_5_sum$Estimate),2) # Fine, if it is only small 

# Shifting the hyperparameters 'down' (specifically the mean and the shape of the beta distribution) 
model_main_6 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ prior("normal(0, 2)")*FIES_est + prior("normal(0, 2)")*Economic_poverty  
    FIES_est ~ prior("normal(-0.5, 2)")*Land_est + prior("normal(0, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(-0.5, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ prior("beta(3, 3)")*Forest_est.1
    Economic_poverty ~~ prior("beta(3, 3)")*Forest_est.1
    Economic_poverty ~~ prior("beta(2, 4)")*Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ prior("beta(2, 4)")*Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ prior("normal(0, 2)")*DOBB + 
    
    # Sex
    prior("normal(0, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(-0.5, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(-0.5, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ prior("normal(-0.25, 9)")*LOCB_EWAF + prior("normal(-0.25, 9)")*LOCB_KADONE + prior("normal(-0.25, 9)")*LOCB_KADTWO + 
    prior("normal(-0.25, 9)")*LOCB_KARO + prior("normal(-0.25, 9)")*LOCB_KYEM + prior("normal(-0.25, 9)")*LOCB_MARA + 
    prior("normal(-0.25, 9)")*LOCB_NONE + prior("normal(-0.25, 9)")*LOCB_NTWO + prior("normal(-0.25, 9)")*LOCB_NYAB + prior("normal(-0.25, 9)")*LOCB_NYAK
    '

# Bayesian SEM  
fit_main_6 <- bsem(model_main_6, data=DF_analy.list[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_6, "fit_main_6.rds")

# # Load the model (if needed)
# fit_main_6 <- readRDS("fit_main_6.rds")

# Calculate the bias associated with using 'down' shifted hyperparameters
fit_main_1_sum <- data.frame(summary(fit_main_1))[,1:2]
fit_main_6_sum <- data.frame(summary(fit_main_6))[,1:2]
round(100*(as.numeric(fit_main_1_sum$Estimate)-as.numeric(fit_main_6_sum$Estimate))/as.numeric(fit_main_6_sum$Estimate),2) # Fine, if it is only small 

######  2.10) 'Is the Bayesian way of interpreting and reporting model results used?' ###### 
# Refer to this article: https://www.rensvandeschoot.com/wp-content/uploads/2017/02/2014-JA-RENS-Depaoli.pdf
# Credibility interval - there is a 95% probability that the true coefficient exists between the credibility interval. 
summary(fit_main_1)


#########  3) Run the model on the ten imputed datasets #########

# The model and priors 
model_main_1 <- ' 
    ### Main regression part
    # Key variables of interest
    PHQ8_est ~ f*prior("normal(0.25, 2)")*FIES_est + g*prior("normal(0.25, 2)")*Economic_poverty  
    FIES_est ~ b*prior("normal(-0.25, 2)")*Land_est + e*prior("normal(0.25, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(-0.25, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ a*prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ c*prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ d*prior("beta(11, 12)")*Land_est
    
    # Health (treated as numeric)
    PHQ8_est ~~ prior("beta(11, 12)")*Health 
    
    ### Covariates 
    # Age 
    PHQ8_est ~ prior("normal(0.25, 2)")*DOBB + 
    
    # Sex
    prior("normal(0.25, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(-0.25, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(-0.25, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0.25, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0.25, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0.52, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    PHQ8_est ~ prior("normal(0, 9)")*LOCB_EWAF + prior("normal(0, 9)")*LOCB_KADONE + prior("normal(0, 9)")*LOCB_KADTWO + 
    prior("normal(0, 9)")*LOCB_KARO + prior("normal(0, 9)")*LOCB_KYEM + prior("normal(0, 9)")*LOCB_MARA + 
    prior("normal(0, 9)")*LOCB_NONE + prior("normal(0, 9)")*LOCB_NTWO + prior("normal(0, 9)")*LOCB_NYAB + prior("normal(0, 9)")*LOCB_NYAK
    
    #### Forest use to depression
    # Indirect effects 
    af := a*f # Forest use -> food insecurity -> depression 
    cg := c*g # Forest use -> economic poverty -> depression
    cef := c*e*f # Forest use -> economic poverty -> food insecurity -> depression 
    cdbf := c*d*b*f # Forest use -> economic poverty -> farm size -> food insecurity -> depression
    
    # Total effects 
    Total.FD.Dep := af + cg + cef + cdbf

    #### Farm size to depression
    # Indirect effects 
    dg := d*g # Farm size -> economic poverty -> depression
    dcaf := d*c*a*f # Farm size -> economic poverty -> forest use -> food insecurity -> depression
    def := d*e*f # Farm size -> economic poverty -> food insecurity -> depression
    bf := b*f # Farm size -> food insecurity -> depression

    # Total effects 
    Total.FS.Dep := dg + dcaf + def + bf
'

# Model list
fit_main_list <- list()

# Number of iterations 
iters <- 10

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){

  # Bayesian SEM
  fit_main_list <- bsem(model_main_1, data=DF_analy.list[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_main_list, "fit_main_list.rds")
#fit_main_list <- readRDS("fit_main_list.rds")

# Extract the variable names 
summary_DF_names <- paste0(fit_main_list@ParTable$lhs, fit_main_list@ParTable$op , fit_main_list@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_main_list_sample <- blavInspect(fit_main_list, what="mcmc")
fit_main_list_sample <- as.matrix(fit_main_list_sample)
fit_1_MCMCbinded_DF <- data.frame(fit_main_list_sample)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_1_MCMCbinded_DF) <- make.unique(summary_DF_names)[1:length(colnames(fit_1_MCMCbinded_DF))]

# Even nicer names for key variables
colnames(fit_1_MCMCbinded_DF) <- revalue(colnames(fit_1_MCMCbinded_DF), c("PHQ8_est~FIES_est" = "Depression ~ Food insecurity", #
                                                "PHQ8_est~Economic_poverty" = "Depression ~ Economic poverty", #
                                                "FIES_est~Land_est" = "Food insecurity ~ Farm size", #
                                                "FIES_est~Economic_poverty" = "Food insecurity ~ Economic poverty", #
                                                "FIES_est~FR.dist" = "Food insecurity ~ Distance to forest reserve",#
                                                "FIES_est~~Forest_est.1" = "Food insecurity ~~ Forest use",#
                                                "Economic_poverty~~Forest_est.1" = "Economic poverty ~~ Forest use", #
                                                "Economic_poverty~~Land_est" = "Economic poverty ~~ Farm size", #
                                                "PHQ8_est~~Health" = "Depression ~ Health",#
                                                "PHQ8_est~DOBB" = "Depression ~ Age",#
                                                "PHQ8_est~SEX" = "Depression ~ Gender (RL = male)", #
                                                "PHQ8_est~EDUCAT" = "Depression ~ Education",#
                                                "PHQ8_est~Soci_est.1" = "Depression ~ Social support (family & sig. other)",#
                                                "PHQ8_est~MSTATUS_Div.wid" = "Depression ~ Divorced or widow/er", #
                                                "PHQ8_est~MSTATUS_Single" = "Depression ~ Never married", #
                                                "PHQ8_est~Alcohol" = "Depression ~ Alcohol consumption", #
                                                "PHQ8_est~SMOKING"  = "Depression ~ Daily smoker",#
                                                "PHQ8_est~LOCB_EWAF" = "Depression ~ Community (Wafala)",
                                                "PHQ8_est~LOCB_KADONE" = "Depression ~ Community (Kadukulu One)",
                                                "PHQ8_est~LOCB_KADTWO" = "Depression ~ Community (Kadukulu Two)", 
                                                "PHQ8_est~LOCB_KARO" = "Depression ~ Community (Karongo)",
                                                "PHQ8_est~LOCB_KYEM" = "Depression ~ Community (Kyempunu)",
                                                "PHQ8_est~LOCB_MARA" = "Depression ~ Community (Maramu)",
                                                "PHQ8_est~LOCB_NONE" = "Depression ~ Community (Nyabyeya One)",
                                                "PHQ8_est~LOCB_NTWO" = "Depression ~ Community (Nyabyeya Two)",
                                                "PHQ8_est~LOCB_NYAB" = "Depression ~ Community (Nyabigoma)",
                                                "PHQ8_est~LOCB_NYAK" = "Depression ~ Community (Nyakafunjo)"
                                                
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_1_MCMCbinded_DF <- fit_1_MCMCbinded_DF[c("Depression ~ Food insecurity", "Depression ~ Economic poverty","Food insecurity ~~ Forest use",
                                              "Food insecurity ~ Farm size","Food insecurity ~ Economic poverty",
                                              "Food insecurity ~ Distance to forest reserve", "Economic poverty ~~ Forest use",
                                              "Economic poverty ~~ Farm size","Depression ~ Age", "Depression ~ Gender (RL = male)",
                                              "PHQ8_est~EDUCAT" = "Depression ~ Education", "Depression ~ Social support (family & sig. other)",
                                              "Depression ~ Divorced or widow/er", "Depression ~ Never married",
                                              "Depression ~ Health","Depression ~ Alcohol consumption","Depression ~ Daily smoker", colnames(fit_1_MCMCbinded_DF)[18:length(colnames(fit_1_MCMCbinded_DF))] )]

# Save pooled samples
saveRDS(fit_1_MCMCbinded_DF, "fit_1_MCMCbinded_DF.rds")



######### 5) Additional analysis associating "thinking too much" and PHQ-8 scores #########
# Create a copy dataset for the supplementary analysis 
DF_analy_sub <- DF_analy.list
names_PHQ8 <- c("PH9INTERST", "PH9FEEL","PH9TROUBL" , "PH9TIRED","PH9APPETIT","PH9BADABT","PH9CONCEN","PH9MOVING")

# Sum PHQ-8 scores
for (i in seq_along(1:length(DF_analy_sub))){
  DF_analy_sub[[i]]$PHQ8 <- rowSums(DF_analy_sub[[i]][,c(names_PHQ8)])
}

# Initial Bayesian ordinal regression  
PHQ8_strong_model <- brm(STRONG  ~ PHQ8, data = DF_analy_sub[[1]], chains = 4, iter = 8000, warmup = 4000, family = cumulative(link = "logit", threshold = "flexible"))

# Trace plots
plot(PHQ8_strong_model, plot.type = "trace")

# Bayesian ordinal regression with all imputed data
PHQ8_strong_list <- brm_multiple(STRONG  ~ PHQ8, data = DF_analy_sub, chains = 4, iter = 8000, warmup = 4000, family = cumulative(link = "logit", threshold = "flexible"))

# Save the model 
saveRDS(PHQ8_strong_list, "PHQ8_strong_list.rds")



######### 6) Repeating the analysis with "thinking too much" as a response #########

# Converting strong thoughts to numeric, and scale and centring 
for (i in seq_along(1:length(DF_analy_sub))){
  DF_analy_sub[[i]]$Strong_num <- PHQ8.rec.num(DF_analy_sub[[i]]$STRONG) 
  DF_analy_sub[[i]]$Strong_num <- scale(DF_analy_sub[[i]]$Strong_num, center = T, scale = T)
}

# The model and priors 
model_sub_1 <- ' 
    ### Main regression part
    # Key variables of interest
    Strong_num ~ f*prior("normal(0.25, 2)")*FIES_est + g*prior("normal(0.25, 2)")*Economic_poverty  
    FIES_est ~ b*prior("normal(-0.25, 2)")*Land_est + e*prior("normal(0.25, 2)")*Economic_poverty  

    # Distance to forest reserve 
    FIES_est ~ prior("normal(-0.25, 2)")*FR.dist

    ### Covariance
    # Between socio-ecological variables
    FIES_est ~~ a*prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ c*prior("beta(12, 11)")*Forest_est.1
    Economic_poverty ~~ d*prior("beta(11, 12)")*Land_est
    
    # Health (treated as numeric)
    Strong_num ~~ prior("beta(11, 12)")*Health 
    
    ### Covariates 
    # Age 
    Strong_num ~ prior("normal(0.25, 2)")*DOBB + 
    
    # Sex
    prior("normal(0.25, 2)")*SEX +

    # Education (treated as numeric)
    prior("normal(-0.25, 2)")*EDUCAT + 
    
    # Social support 
    prior("normal(-0.25, 2)")*Soci_est.1 + 
    
    # Marital status (RL = MSTATUS_Mar.pol)
    prior("normal(0.25, 2)")*MSTATUS_Div.wid + prior("normal(0, 9)")*MSTATUS_Single +
    
    # Alchohol consumption  
    prior("normal(0.25, 2)")*Alcohol + 
    
    # Smoking 
    prior("normal(0.52, 2)")*SMOKING 
  
    ### Location (RL = LOCB_NTC)
    Strong_num ~ prior("normal(0, 9)")*LOCB_EWAF + prior("normal(0, 9)")*LOCB_KADONE + prior("normal(0, 9)")*LOCB_KADTWO + 
    prior("normal(0, 9)")*LOCB_KARO + prior("normal(0, 9)")*LOCB_KYEM + prior("normal(0, 9)")*LOCB_MARA + 
    prior("normal(0, 9)")*LOCB_NONE + prior("normal(0, 9)")*LOCB_NTWO + prior("normal(0, 9)")*LOCB_NYAB + prior("normal(0, 9)")*LOCB_NYAK
    
    #### Forest use to depression
    # Indirect effects 
    af := a*f # Forest use -> food insecurity -> depression 
    cg := c*g # Forest use -> economic poverty -> depression
    cef := c*e*f # Forest use -> economic poverty -> food insecurity -> depression 
    cdbf := c*d*b*f # Forest use -> economic poverty -> farm size -> food insecurity -> depression
    
    # Total effects 
    Total.FD.Dep := af + cg + cef + cdbf

    #### Farm size to depression
    # Indirect effects 
    dg := d*g # Farm size -> economic poverty -> depression
    dcaf := d*c*a*f # Farm size -> economic poverty -> forest use -> food insecurity -> depression
    def := d*e*f # Farm size -> economic poverty -> food insecurity -> depression
    bf := b*f # Farm size -> food insecurity -> depression

    # Total effects 
    Total.FS.Dep := dg + dcaf + def + bf
'

# Run a test model 
fit_sub_model <- bsem(model_sub_1, data=DF_analy_sub[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Trace plots
plot(fit_sub_model, plot.type = "trace", 1:10)
saveRDS(fit_sub_model, "fit_sub_model.rds")

### Imputed dataset ### 

# Model list
fit_sub_list <- list()

# Number of iterations 
iters <- 10

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){
  
  # Bayesian SEM
  fit_sub_list <- bsem(model_sub_1, data=DF_analy_sub[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_sub_list, "fit_sub_list.rds")
#fit_sub_list <- readRDS("fit_sub_list.rds")

# Extract the variable names 
summary_DF_names_sub <- paste0(fit_sub_list@ParTable$lhs, fit_sub_list@ParTable$op , fit_sub_list@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_sub_list_sample <- blavInspect(fit_sub_list, what="mcmc")
fit_sub_list_sample <- as.matrix(fit_sub_list_sample)
fit_2_MCMCbinded_DF <- data.frame(fit_sub_list_sample)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF) <- make.unique(summary_DF_names_sub)[1:length(colnames(fit_2_MCMCbinded_DF))]

# Even nicer names for key variables
colnames(fit_2_MCMCbinded_DF) <- revalue(colnames(fit_2_MCMCbinded_DF), c("Strong_num~FIES_est" = "Strong thoughts ~ Food insecurity", #
                                                                          "Strong_num~Economic_poverty" = "Strong thoughts ~ Economic poverty", #
                                                                          "FIES_est~Land_est" = "Food insecurity ~ Farm size", #
                                                                          "FIES_est~Economic_poverty" = "Food insecurity ~ Economic poverty", #
                                                                          "FIES_est~FR.dist" = "Food insecurity ~ Distance to forest reserve",#
                                                                          "FIES_est~~Forest_est.1" = "Food insecurity ~~ Forest use",#
                                                                          "Economic_poverty~~Forest_est.1" = "Economic poverty ~~ Forest use", #
                                                                          "Economic_poverty~~Land_est" = "Economic poverty ~~ Farm size", #
                                                                          "Strong_num~~Health" = "Strong thoughts ~ Health",#
                                                                          "Strong_num~DOBB" = "Strong thoughts ~ Age",#
                                                                          "Strong_num~SEX" = "Strong thoughts ~ Gender (RL = male)", #
                                                                          "Strong_num~EDUCAT" = "Strong thoughts ~ Education",#
                                                                          "Strong_num~Soci_est.1" = "Strong thoughts ~ Social support (family & sig. other)",#
                                                                          "Strong_num~MSTATUS_Div.wid" = "Strong thoughts ~ Divorced or widow/er", #
                                                                          "Strong_num~MSTATUS_Single" = "Strong thoughts ~ Never married", #
                                                                          "Strong_num~Alcohol" = "Strong thoughts ~ Alcohol consumption", #
                                                                          "Strong_num~SMOKING"  = "Strong thoughts ~ Daily smoker",#
                                                                          "Strong_num~LOCB_EWAF" = "Strong thoughts ~ Community (Wafala)",
                                                                          "Strong_num~LOCB_KADONE" = "Strong thoughts ~ Community (Kadukulu One)",
                                                                          "Strong_num~LOCB_KADTWO" = "Strong thoughts ~ Community (Kadukulu Two)", 
                                                                          "Strong_num~LOCB_KARO" = "Strong thoughts ~ Community (Karongo)",
                                                                          "Strong_num~LOCB_KYEM" = "Strong thoughts ~ Community (Kyempunu)",
                                                                          "Strong_num~LOCB_MARA" = "Strong thoughts ~ Community (Maramu)",
                                                                          "Strong_num~LOCB_NONE" = "Strong thoughts ~ Community (Nyabyeya One)",
                                                                          "Strong_num~LOCB_NTWO" = "Strong thoughts ~ Community (Nyabyeya Two)",
                                                                          "Strong_num~LOCB_NYAB" = "Strong thoughts ~ Community (Nyabigoma)",
                                                                          "Strong_num~LOCB_NYAK" = "Strong thoughts ~ Community (Nyakafunjo)"
                                                                          
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_2_MCMCbinded_DF <- fit_2_MCMCbinded_DF[c("Strong thoughts ~ Food insecurity", "Strong thoughts ~ Economic poverty","Food insecurity ~~ Forest use",
                                             "Food insecurity ~ Farm size","Food insecurity ~ Economic poverty",
                                             "Food insecurity ~ Distance to forest reserve", "Economic poverty ~~ Forest use",
                                             "Economic poverty ~~ Farm size","Strong thoughts ~ Age", "Strong thoughts ~ Gender (RL = male)",
                                             "PHQ8_est~EDUCAT" = "Strong thoughts ~ Education", "Strong thoughts ~ Social support (family & sig. other)",
                                             "Strong thoughts ~ Divorced or widow/er", "Strong thoughts ~ Never married",
                                             "Strong thoughts ~ Health","Strong thoughts ~ Alcohol consumption","Strong thoughts ~ Daily smoker", colnames(fit_2_MCMCbinded_DF)[18:length(colnames(fit_2_MCMCbinded_DF))] )]

# Save pooled samples
saveRDS(fit_2_MCMCbinded_DF, "fit_2_MCMCbinded_DF.rds")
