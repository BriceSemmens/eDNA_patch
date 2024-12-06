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

## Lets do effects of water volume, depth & method capture, with site random effect  ######################################################################

# Define the NIMBLE model

edna_code_vol_depth_meth_randCap <- nimbleCode({
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i] ~ dnorm(cap_prob_hat,cap_prob_SD)
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
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection) 
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection ~ dbeta(1, 1)
  b_meth[1] <- 0
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
                                       niter = 2000000, 
                                       nburnin = 100000, 
                                       thin = 100, 
                                       nchains = 3,
                                       summary=TRUE,
                                       samplesAsCodaMCMC = TRUE,
                                       WAIC = TRUE)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
MCMCsummary(edna_code_vol_depth_meth_randCap.run$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(edna_code_vol_depth_meth_randCap.run$samples)

#put posteriors into the work environment
attach.nimble(edna_code_vol_depth_meth_randCap.run$samples)


library(ggplot2)
library(boot)

n.post<-length(b_depth) #get number of posterior draws

# Plot the linear relationship between capture probability and volume filtered

plot.stor<-matrix(data=NA,n.post,100) #make for work
# Create a range of depth values for samples (assuming standardized values)
vol_seq <- seq(min(biosamp_Volume_filt_mL), max(biosamp_Volume_filt_mL), length.out = 100)
for (i in 1:100){
  for (j in 1:n.post){
    plot.stor[j,i]<- inv.logit(b_meth[j,2]+ b_vol[j]*vol_seq[i])
  }
}
cap_prob_mean <- apply( plot.stor, 2, function(x) quantile(x, 0.5))
cap_prob_lower<- apply( plot.stor, 2, function(x) quantile(x, 0.10))
cap_prob_upper<- apply( plot.stor, 2, function(x) quantile(x, 0.90))

plot_data <- data.frame(
  vol.Filtered = (vol_seq+mean(as.numeric(biosamp_dat$Volume_filt_mL))),
  CaptureProbMean = cap_prob_mean,
  CaptureProbLower = cap_prob_lower,
  CaptureProbUpper = cap_prob_upper
)

p <- ggplot(plot_data, aes(x = vol.Filtered)) +
  geom_line(aes(y = CaptureProbMean), color = "blue") +
  geom_ribbon(aes(ymin = CaptureProbLower, ymax = CaptureProbUpper), fill = "grey", alpha = 0.5) +
  labs(x = "Volume filtered", y = "Probability of Capture") +
  theme_minimal() +
  ggtitle("Relationship Between Probability of Capture and Volume Filtered")

# Print the plot
print(p)

# Plot the linear relationship between capture probability and depth of sample

plot.stor2<-matrix(data=NA,n.post,100) #make for work
# Create a range of depth values for samples (assuming standardized values)
depth_seq <- seq(min(biosamp_depth_Depth_m), max(biosamp_depth_Depth_m), length.out = 100)
for (i in 1:100){
  for (j in 1:n.post){
    plot.stor2[j,i]<- inv.logit(b_meth[j,2]+ b_depth[j]*depth_seq[i])
  }
}
cap_prob_mean2 <- apply( plot.stor2, 2, function(x) quantile(x, 0.5))
cap_prob_lower2<- apply( plot.stor2, 2, function(x) quantile(x, 0.10))
cap_prob_upper2<- apply( plot.stor2, 2, function(x) quantile(x, 0.90))

plot_data2 <- data.frame(
  depth.samp = (depth_seq+mean(as.numeric(biosamp_dat$Depth_m))),
  CaptureProbMean = cap_prob_mean2,
  CaptureProbLower = cap_prob_lower2,
  CaptureProbUpper = cap_prob_upper2
)

p2 <- ggplot(plot_data2, aes(x = depth.samp)) +
  geom_line(aes(y = CaptureProbMean), color = "blue") +
  geom_ribbon(aes(ymin = CaptureProbLower, ymax = CaptureProbUpper), fill = "grey", alpha = 0.5) +
  labs(x = "Depth Sampled", y = "Probability of Capture") +
  theme_minimal() +
  ggtitle("Relationship Between Probability of Capture and Water Sample Depth (m)")

# Print the plot
print(p2)

# now let's plot the difference in RREAS and GEMCAP capture prob for average depth and volume
gemcap.detect<-inv.logit(cap_prob_hat)
rreas.detect<-inv.logit(b_meth[,2]+ cap_prob_hat)
boxplot(cbind(gemcap.detect,rreas.detect))
