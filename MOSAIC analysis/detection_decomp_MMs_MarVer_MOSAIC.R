library(MCMCvis)
library(boot)
library(tidyverse)
library(mcmcplots)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(nimble)

# Import the CSV file
mm.data <- read.csv("SR2408_mammals_forBXS.csv")# Import the CSV file

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(uniqiue_filter_numeric) %>%
  slice(1) %>%
  ungroup()

# make data, index vectors and constants for nimble work
biosamp_station_index<- biosamp_dat$Trans.numeric
lab_index<-biosamp_dat$Lab.numeric
biosamp_method_index<- biosamp_dat$samp.meth.num
Y_biosamp_index<- mm.data$uniqiue_filter_numeric
biosamp_Volume_filt_mL <- as.vector(scale(biosamp_dat$Volume_filt_mL)) #centered water volumes
Y_method_index<- mm.data$samp.meth.num
Y_station_index<-mm.data$Trans.numeric
lab_index<-mm.data$Lab.numeric

N<-dim(mm.data)[1]
n_sites<-length(unique(mm.data$Trans.numeric))
n_biosamples<-length(unique(mm.data$uniqiue_filter_numeric))
n_methods<-length(unique(mm.data$samp.meth.num))

Y<-mm.data$Whale.PA #HEY YOU! THIS IS BIG/RARE WHALES ONLY. IF YOU WANT ALL MMS, CHANGE.

## Lets do effects of water volume, with site random effect  ######################################################################

# Define the NIMBLE model

edna_code_vol_depth_meth_randCap <- nimbleCode({

  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i] ~ dnorm(cap_prob_hat,cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- cap_prob_logit[biosamp_station_index[i]] +   # categorical differences in capture method # removed [biosamp_method_index[i]]
      b_vol * biosamp_Volume_filt_mL[i]  +    # continuous effect of volume
      b_meth[biosamp_method_index[i]]
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
    
    #Impute missing volume(s):
    biosamp_Volume_filt_mL[i] ~ dnorm(0,1)
  }

  
  # Likelihood of detecting in each lab tech replicate
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * detect_prob[lab_index[i]]) 
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(8, 1) #informative (high likelihood occurrence)
  b_vol ~ dnorm(0, 1)
  b_meth[1] <- 0
  b_meth[2] ~dnorm(0,1.7)
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  detect_prob[1] ~ dbeta(1,1) #SIO specific detect prob
  detect_prob[2] ~ dbeta(1,1) #UW specific detect prob
  
  
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
  lab_index = lab_index
)

# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.9,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  #bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_vol = 0     # Initial value for volume coefficient
)

# Run NIMBLE model
edna_code_vol_depth_meth_randCap.run <- nimbleMCMC(code = edna_code_vol_depth_meth_randCap, 
                                       constants = constants, 
                                       data = data, 
                                       inits = inits,
                                       niter = 200000, 
                                       nburnin = 10000, 
                                       thin = 100, 
                                       nchains = 3,
                                       summary=TRUE,
                                       samplesAsCodaMCMC = TRUE,
                                       WAIC = TRUE)

# Gelman-Rubin diagnostic
MCMCsummary(edna_code_vol_depth_meth_randCap.run$samples)

# Visualize MCMC chains
mcmcplot(edna_code_vol_depth_meth_randCap.run$samples)

# --- Posterior analysis and plotting ---
attach.nimble(edna_code_vol_depth_meth_randCap.run$samples)
n.post <- length(b_depth)

# --- Define sequences for depth and volume ---

vol_min <- min(as.numeric(biosamp_dat$Volume_filt_mL),na.rm = TRUE)
vol_max <- max(as.numeric(biosamp_dat$Volume_filt_mL),na.rm = TRUE)
n_vol_steps <- 100
vol_seq <- seq(from = vol_min, to = vol_max, length.out = n_vol_steps)

# --- Calculate predicted probabilities ---
plot.stor <- matrix(data = NA, n.post, length(vol_seq))

for (i in 1:length(vol_seq)) {
  for (j in 1:n.post) {
    plot.stor[j, i] <- inv.logit(b_meth[j,2] + b_vol[j] * (vol_seq[i]-mean(as.numeric(biosamp_dat$Volume_filt_mL))))
  }
}

# --- Create plot data frames ---
plot_data_lineribbon_depth <- data.frame(
  x_value = rep(depth_seq, each = n.post),
  value = as.vector(plot.stor2)
)

plot_data_lineribbon_vol <- data.frame(
  x_value = rep(vol_seq, each = n.post),
  value = as.vector(plot.stor)
)

# --- Create plots ---
p2 <- ggplot(plot_data_lineribbon_depth, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth Sampled", y = "Probability of Capture", title = "Probability of Capture Covariates") +
  theme_minimal()

p1 <- ggplot(plot_data_lineribbon_vol, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume Sampled", y = "Probability of Capture") +
  theme_minimal()

# --- Combine plots ---
combined_data <- rbind(
  data.frame(variable = "Depth (m)", plot_data_lineribbon_depth),
  data.frame(variable = "Volume Filtered (mL)", plot_data_lineribbon_vol)
)

p <- ggplot(combined_data, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  facet_wrap(~ variable, scales = "free_x", nrow = 2) +
  labs(x = "", y = "Probability of Capture", title = "Prob. of Capture Covariates") +
  theme_minimal()

print(p)

# Let's print the probabilities of capture, detection and occurrence -- all in one plot

# --- Plot 1: Underway vs. Carboy ---
underway.detect <- as.numeric(inv.logit(cap_prob_hat))
carboy.detect <- as.numeric(inv.logit(b_meth[,2] + cap_prob_hat))

df_detect <- data.frame(underway = underway.detect, carboy = carboy.detect) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Probability")

medians_detect <- df_detect %>%
  group_by(Method) %>%
  summarize(Median = median(Probability))

viridis_cols <- viridis(2, begin = 0.3, end = 0.7)
names(viridis_cols) <- c("carboy", "underway")

p1 <- ggplot(df_detect, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  geom_vline(data = medians_detect, aes(xintercept = Median, color = Method), linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = Inf, label = names(viridis_cols), color = viridis_cols, hjust = 1.1, vjust = c(3, 4.5), size = 5) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Capture Probability", y = "Density", title = "underway vs. carboy") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
p1 <- p1 + scale_x_continuous(limits = c(0, 1)) # Set x-axis limits

# --- Plot 2: Probability of Occurrence ---
df_prob_occurrence <- data.frame(Occurrence = as.numeric(prob_occurrence))

p2 <- ggplot(df_prob_occurrence, aes(x = Occurrence)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey70", color = "black", alpha = 0.6) +
  geom_density(color = viridis(1), size = 1) +
  geom_vline(xintercept = median(prob_occurrence), linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = median(prob_occurrence), y = Inf, label = paste("Median =", round(median(prob_occurrence), 3)), vjust = 1.5, hjust = -0.1, color = "black") +
  labs(x = "Probability of Occurrence", y = "Density", title = "Probability of Occurrence") +
  theme_bw() + theme(panel.grid = element_blank())
p2 <- p2 + scale_x_continuous(limits = c(0, 1)) # Set x-axis limits

# --- Plot 3: Probability of Detection (given occurrence) ---
df_prob_detection <- data.frame(SIO = as.numeric(detect_prob[, 1]), UW = as.numeric(detect_prob[, 2])) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Probability")

medians_detection <- df_prob_detection %>%
  group_by(Method) %>%
  summarize(Median = median(Probability))

viridis_cols_detection <- viridis(2, begin = 0.3, end = 0.7)
names(viridis_cols_detection) <- c("SIO", "UW")

p3 <- ggplot(df_prob_detection, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  geom_vline(data = medians_detection, aes(xintercept = Median, color = Method), linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = Inf, label = names(viridis_cols_detection), color = viridis_cols_detection, hjust = 1.1, vjust = c(3, 4.5), size = 5) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Probability of Detection", y = "Density", title = "Probability of Detection (SIO vs. UW)") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
p3 <- p3 + scale_x_continuous(limits = c(0, 1)) # Set x-axis limits

# Combine the plots using patchwork - Probability of Detection at the bottom
combined_plot <- p1 / p3 / p2

print(combined_plot)