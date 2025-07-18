# Import the CSV file
mm.data <- read.csv("intercal_DL_metadata_111224-noCs_NEW.csv")# Import the CSV file

# remove the underway sample data from this analysis
mm.data <- mm.data[mm.data$Collection_method != "UW", ]

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(uniqiue_biorep_numeric) %>%
  slice(1) %>%
  ungroup()

# make data, index vectors and constants for nimble work
# Import the CSV file
mm.data <- read.csv("intercal_DL_metadata_111224-noCs_NEW.csv")# Import the CSV file

# remove the underway sample data from this analysis
mm.data <- mm.data[mm.data$Collection_method != "UW", ]

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(uniqiue_biorep_numeric) %>%
  slice(1) %>%
  ungroup()

# make data, index vectors and constants for nimble work
biosamp_station_index<- biosamp_dat$site.numeric
biosamp_method_index<- biosamp_dat$Collection_method_numeric
Y_biosamp_index<- mm.data$uniqiue_biorep_numeric
biosamp_Volume_filt_mL <- as.numeric(biosamp_dat$Volume_filt_mL) - mean(as.numeric(biosamp_dat$Volume_filt_mL)) #centered water volumes
biosamp_depth_Depth_m <- as.numeric(biosamp_dat$Depth_m) - mean(as.numeric(biosamp_dat$Depth_m)) #centered sample depths
Y_method_index<- mm.data$Collection_method_numeric
Y_station_index<-mm.data$site.numeric

N<-dim(mm.data)[1]
n_sites<-length(unique(mm.data$site.numeric))
n_biosamples<-length(unique(mm.data$uniqiue_biorep_numeric))
n_methods<-length(unique(mm.data$Collection_method_numeric ))

Y<-mm.data$mm_detect
# Define the NIMBLE model
edna_code_basic <- nimbleCode({
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
  }
  
  for (i in 1:n_biosamples) {
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[biosamp_method_index[i]])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) #removed [Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  prob_detection ~ dbeta(1, 1)
  
  for (i in 1:n_methods) { # probability of capture for each method in data
    prob_capture[i] ~ dbeta(1, 1)
    #prob_detection[i] ~ dbeta(1, 1)
  }
})

# Load the nimble package
library(nimble)

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(Y = Y)

# Define initial values for the parameters
inits <- list(
  prob_occurrence = 0.5,
  prob_capture = rep(0.5, n_methods),
  #prob_detection = rep(0.5, n_methods),
  site_occurrence = rep(1, n_sites),
  bio_capture = rep(1, n_biosamples)
)

# Run NIMBLE model
edna_model_basic.run <- nimbleMCMC(code = edna_code_basic, 
                          constants = constants, 
                          data = data, 
                          inits = inits,
                          niter = 10000, 
                          nburnin = 2000, 
                          thin = 10, 
                          nchains = 3,
                          summary=TRUE,
                          samplesAsCodaMCMC = TRUE,
                          WAIC = TRUE)
# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_model_basic.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_model_basic.run$samples)

# ---------------------------------------------------------------
# To do:
#  - run a model where we make occurrence probability a function of water volume
#  - run a model where we make occurrence probability a function of depth
#  - use WAIC to evaluate whether there should be different detection probabilities for each method
#  - make detection probability a  site specific random effect (account for differences in MM edna between sites)


# Lets do water volume & depth on effects capture  ######################################################################
 
# Define the NIMBLE model
edna_code_vol_depth <- nimbleCode({
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- a_logit_capture +   # categorical differences in capture method # removed [biosamp_method_index[i]]
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]       # continuous effect of volume
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) #[Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  a_logit_capture ~ dnorm(0,1.7)
  
  prob_detection ~ dbeta(1, 1)
 # for (i in 1:n_methods) {
  # a_logit_capture[i] <- dnorm(0,1.7)
#   prob_detection[i] ~ dbeta(1, 1)
 # }
})



# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_model_vol_depth.run <- nimbleMCMC(code = edna_code_vol_depth, 
                         constants = constants, 
                         data = data, 
                         inits = inits,
                         niter = 100000, 
                         nburnin = 20000, 
                         thin = 100, 
                         nchains = 3,
                         summary=TRUE,
                         samplesAsCodaMCMC = TRUE,
                         WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_model_vol_depth.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_model_vol_depth.run$samples)


# now do volume, depth and method effects on capture  #####################################

# Lets do water volume & depth
# Define the NIMBLE model
edna_code_vol_depth_meth <- nimbleCode({
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- a_logit_capture[biosamp_method_index[i]] +   # categorical differences in capture method
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]       # continuous effect of volume
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) # removed [Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection ~ dbeta(1, 1)
  
  for (i in 1:n_methods) {
    a_logit_capture[i] ~ dnorm(0,1.7)
   # prob_detection[i] ~ dbeta(1, 1)
  }
})

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_code_vol_depth_meth.run <- nimbleMCMC(code = edna_code_vol_depth_meth, 
                                   constants = constants, 
                                   data = data, 
                                   inits = inits,
                                   niter = 100000, 
                                   nburnin = 20000, 
                                   thin = 100, 
                                   nchains = 3,
                                   summary=TRUE,
                                   samplesAsCodaMCMC = TRUE,
                                   WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_vol_depth_meth.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_vol_depth_meth.run$samples)

## now do volume, depth and method effects on capture + random effects for detection by site #################
# Define the NIMBLE model
edna_code_vol_depth_meth_randDetect <- nimbleCode({
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    
    # Site-level detection probability (under assumption that + or - DNA abund at a site impacts detection)
    prob_detection_logit[i] ~ dnorm(prob_D_hat, prob_D_SD)
    logit(prob_detection[i])<-prob_detection_logit[i] #put in proportion space
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- a_logit_capture[biosamp_method_index[i]] +   # categorical differences in capture method
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]       # continuous effect of volume
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_station_index[i]])
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_D_hat ~ dnorm(0,1.7) 
  prob_D_SD~ dexp(1)
  
  for (i in 1:n_methods) {
    a_logit_capture[i] ~ dnorm(0,1.7)
  }
})

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index,
  Y_station_index = Y_station_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  prob_detection = rep(0.5, n_sites),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_code_vol_depth_meth_randDetect.run <- nimbleMCMC(code = edna_code_vol_depth_meth_randDetect, 
                                        constants = constants, 
                                        data = data, 
                                        inits = inits,
                                        niter = 100000, 
                                        nburnin = 20000, 
                                        thin = 100, 
                                        nchains = 3,
                                        summary=TRUE,
                                        samplesAsCodaMCMC = TRUE,
                                        WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_vol_depth_meth_randDetect.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_vol_depth_meth_randDetect.run$samples)

##  make capture probability a random effect of site (some sites have MMs everywhere) ##

edna_code_basic_randCap <- nimbleCode({
  
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i]~dnorm(cap_prob_hat,cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[biosamp_station_index[i]])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) #removed [Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  prob_detection ~ dbeta(1, 1)
  
 # for (i in 1:n_methods) { # probability of capture for each method in data
 #  prob_capture[i] ~ dbeta(1, 1)
    #prob_detection[i] ~ dbeta(1, 1)
 # }
})

# Load the nimble package
library(nimble)

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(Y = Y)

# Define initial values for the parameters
inits <- list(
  prob_occurrence = 0.5,
  #prob_capture = rep(0.5, n_methods),
  #prob_detection = rep(0.5, n_methods),
  site_occurrence = rep(1, n_sites),
  bio_capture = rep(1, n_biosamples)
)

# Run NIMBLE model
edna_code_basic_randCap.run <- nimbleMCMC(code = edna_code_basic_randCap, 
                                   constants = constants, 
                                   data = data, 
                                   inits = inits,
                                   niter = 10000, 
                                   nburnin = 2000, 
                                   thin = 10, 
                                   nchains = 3,
                                   summary=TRUE,
                                   samplesAsCodaMCMC = TRUE,
                                   WAIC = TRUE)
# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_basic_randCap.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_basic_randCap.run$samples)

## Lets do water volume & depth on effects capture, with site random effect  ######################################################################

# Define the NIMBLE model
edna_code_vol_depth_randCap <- nimbleCode({
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i]~dnorm(cap_prob_hat,cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- cap_prob_logit[biosamp_station_index[i]] +   # categorical differences in capture method # removed [biosamp_method_index[i]]
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]       # continuous effect of volume
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) #[Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
    # `Y_method_index[i]` -- index vector that specifies which sample method each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  
  prob_detection ~ dbeta(1, 1)
  # for (i in 1:n_methods) {
  # a_logit_capture[i] <- dnorm(0,1.7)
  #   prob_detection[i] ~ dbeta(1, 1)
  # }
})



# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  #bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_code_vol_depth_randCap.run <- nimbleMCMC(code = edna_code_vol_depth_randCap, 
                                       constants = constants, 
                                       data = data, 
                                       inits = inits,
                                       niter = 100000, 
                                       nburnin = 20000, 
                                       thin = 100, 
                                       nchains = 3,
                                       summary=TRUE,
                                       samplesAsCodaMCMC = TRUE,
                                       WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_vol_depth_randCap.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_vol_depth_randCap.run$samples)

## Lets do effects of water volume, depth & method capture, with site random effect  ######################################################################

# Define the NIMBLE model
edna_code_vol_depth_meth_randCap <- nimbleCode({
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i]~dnorm(cap_prob_hat,cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- cap_prob_logit[biosamp_station_index[i]] +   # categorical differences in capture method # removed [biosamp_method_index[i]]
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]  +    # continuous effect of volume
      b_meth[biosamp_method_index[i]]
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) #[Y_method_index[i]]
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection ~ dbeta(1, 1)
  b_meth[1] ~dnorm(0,1.7)
  b_meth[2] ~dnorm(0,1.7)
  
})

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_methods = n_methods,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_method_index = Y_method_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  #bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_code_vol_depth_meth_randCap.run <- nimbleMCMC(code = edna_code_vol_depth_meth_randCap, 
                                                   constants = constants, 
                                                   data = data, 
                                                   inits = inits,
                                                   niter = 100000, 
                                                   nburnin = 20000, 
                                                   thin = 100, 
                                                   nchains = 3,
                                                   summary=TRUE,
                                                   samplesAsCodaMCMC = TRUE,
                                                   WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_vol_depth_meth_randCap.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_vol_depth_meth_randCap.run$samples)

# compare WAICs (just a sniff test to make sure nothing is horrifically mis-specified)
edna_model_basic.run$WAIC$WAIC
edna_model_vol_depth.run$WAIC$WAIC
edna_code_vol_depth_meth.run$WAIC$WAIC
edna_code_vol_depth_meth_randDetect.run$WAIC$WAIC
edna_code_basic_randCap.run$WAIC$WAIC
edna_code_vol_depth_randCap.run$WAIC$WAIC
edna_code_vol_depth_meth_randCap.run$WAIC$WAIC

# Ok, based on WAIC exploration, and given that we care about sampling covariates (depth and vol)
# let's proceed with edna_code_vol_depth_meth_randCap.run as it is essentially the best model (lowest WAIC)
