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

source("./Scripts/attach.nimble_v2.R")

# Import the CSV file
mm.data <- read.csv("intercal_ALL_metadata_12.30.24-noCs_NEW.csv")

# remove the underway sample data from this analysis
mm.data <- mm.data[mm.data$Collection_method != "UW", ]

# remove the passive filter sample data from this analysis
mm.data <- mm.data[mm.data$Collection_method != "PF", ]

# because I've removed some of the unique biosample reference numbers with the above 
# methods removal, I need to renumber the unique biosamples so that they are consecutive
# and start with 1. I do this using the factor trick. 
mm.data$uniqiue_biorep_numeric<-as.numeric(as.factor(mm.data$uniqiue_biorep_numeric))

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
Y_primer_index<-mm.data$Prmer_numeric

N<-dim(mm.data)[1]
n_sites<-length(unique(mm.data$site.numeric))
n_biosamples<-length(unique(mm.data$uniqiue_biorep_numeric))
n_methods<-length(unique(mm.data$Collection_method_numeric ))
n_primers<-length(unique(mm.data$Prmer_numeric))

Y<-mm.data$delphinus_all #HEY YOU! THIS IS DELPHINIDAE ONLY. IF YOU WANT ALL MMS, CHANGE.

#############
# ORIGINAL MODEL: Single occupancy probability with covariates on capture
#############

# Define the original NIMBLE model (UNCHANGED)
edna_code_vol_depth_meth_randCap <- nimbleCode({
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occupancy probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i] ~ dnorm(cap_prob_hat,cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    logit(prob_capture[i]) <- cap_prob_logit[biosamp_station_index[i]] +
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]  +    # continuous effect of volume
      b_meth[biosamp_method_index[i]]
    
    # biosample-level occupancy probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection[1] ~ dbeta(1, 1) #detection prob for dloop_original
  prob_detection[2] ~ dbeta(1, 1) #detection prob for MiFish
  b_meth[1] <- 0
  b_meth[2] ~dnorm(0,1.7)
  b_meth[3] ~dnorm(0,1.7)
})

# Define the constants, data, and initial values for original model
# UPDATED: Removed Y_method_index as it was not used in the nimble model
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_primer_index = Y_primer_index
)

data <- list(
  Y = Y,
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL,
  biosamp_depth_Depth_m = biosamp_depth_Depth_m
)

# UPDATED: Using a function for initial values to prevent errors.
# This function ensures that if Y=1 for any observation, the corresponding
# latent states (bio_capture, site_occurrence) are initialized to 1.
inits_fn_orig <- function(){
  # Determine if a biosample had any positive detections
  max_obs_per_biosample <- tapply(Y, Y_biosamp_index, max)
  bio_capture_init <- rep(0, n_biosamples)
  bio_capture_init[as.numeric(names(max_obs_per_biosample))] <- max_obs_per_biosample
  
  # Determine if a site had any positive biosamples
  site_occurrence_init_raw <- tapply(bio_capture_init, biosamp_station_index, max)
  site_occurrence_init <- rep(0, n_sites)
  site_occurrence_init[as.numeric(names(site_occurrence_init_raw))] <- site_occurrence_init_raw
  
  list(
    # Latent states (from before)
    bio_capture = bio_capture_init,
    site_occurrence = site_occurrence_init,
    
    # Model parameters (previously incomplete)
    prob_occurrence = 0.5,
    b_depth = 0,
    b_vol = 0,
    cap_prob_hat = 0,
    cap_prob_SD = 1,
    cap_prob_logit = rep(0, n_sites),
    b_meth = c(NA, 0, 0), # NA for the first fixed element, 0 for the others
    prob_detection = rep(0.5, n_primers)
  )
}


# Run original NIMBLE model
edna_code_vol_depth_meth_randCap.run <- nimbleMCMC(code = edna_code_vol_depth_meth_randCap, 
                                                   constants = constants, 
                                                   data = data, 
                                                   inits = inits_fn_orig, # Using the function here
                                                   niter = 200000, 
                                                   nburnin = 10000, 
                                                   thin = 100, 
                                                   nchains = 3,
                                                   summary=TRUE,
                                                   samplesAsCodaMCMC = TRUE,
                                                   WAIC = TRUE)


# Diagnostics
MCMCsummary(edna_code_vol_depth_meth_randCap.run$samples)
mcmcplot(edna_code_vol_depth_meth_randCap.run$samples)

# --- Original Model 4-Panel Plot ---
attach.nimble(edna_code_vol_depth_meth_randCap.run$samples)
n.post <- length(b_depth)

# Define sequences for depth and volume
depth_seq <- seq(from = min(as.numeric(biosamp_dat$Depth_m)), 
                 to = max(as.numeric(biosamp_dat$Depth_m)), 
                 length.out = 100)
vol_seq <- seq(from = min(as.numeric(biosamp_dat$Volume_filt_mL)), 
               to = max(as.numeric(biosamp_dat$Volume_filt_mL)), 
               length.out = 100)

# Calculate predicted probabilities
plot.stor_depth <- matrix(data = NA, n.post, length(depth_seq))
plot.stor_vol <- matrix(data = NA, n.post, length(vol_seq))

for (i in 1:length(depth_seq)) {
  for (j in 1:n.post) {
    plot.stor_depth[j, i] <- inv.logit(b_meth[j,2] + b_depth[j] * (depth_seq[i] - mean(as.numeric(biosamp_dat$Depth_m))))
  }
}

for (i in 1:length(vol_seq)) {
  for (j in 1:n.post) {
    plot.stor_vol[j, i] <- inv.logit(b_meth[j,2] + b_vol[j] * (vol_seq[i] - mean(as.numeric(biosamp_dat$Volume_filt_mL))))
  }
}

# Create plot data frames
plot_data_depth_orig <- data.frame(
  x_value = rep(depth_seq, each = n.post),
  value = as.vector(plot.stor_depth)
)

plot_data_vol_orig <- data.frame(
  x_value = rep(vol_seq, each = n.post),
  value = as.vector(plot.stor_vol)
)

# Create data for method and detection comparisons
gemcap.detect <- as.numeric(inv.logit(cap_prob_hat))
rreas.detect <- as.numeric(inv.logit(b_meth[,2] + cap_prob_hat))

df_detect_orig <- data.frame(GEMCAP = gemcap.detect, RREAS = rreas.detect) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Probability")

df_prob_detection_orig <- data.frame(
  Probability = c(prob_detection[, 1], prob_detection[, 2]),
  Method = rep(c("dloop", "MiFish"), each = nrow(prob_detection))
)

# Create 4-panel plot for original model
p1_orig <- ggplot(plot_data_depth_orig, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth (m)", y = "Probability of Capture", title = "Original: Depth Effect on Capture") +
  theme_minimal()

p2_orig <- ggplot(plot_data_vol_orig, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume Filtered (mL)", y = "Probability of Capture", title = "Original: Volume Effect on Capture") +
  theme_minimal()

p3_orig <- ggplot(df_detect_orig, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Capture Probability", y = "Density", title = "Original: RREAS vs. GEMCAP") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 1))

p4_orig <- ggplot(df_prob_detection_orig, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Detection Probability", y = "Density", title = "Original: Dloop vs. MiFish") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 1))

original_4panel <- (p1_orig | p2_orig) / (p3_orig | p4_orig)
print("=== ORIGINAL MODEL 4-PANEL PLOT ===")
print(original_4panel)

#############
# SITE x DEPTH MODEL: Site-depth specific occupancy states
#############

# Create site-depth combinations for occupancy states
biosamp_dat$site_depth_combo <- paste(biosamp_dat$site.numeric, 
                                      round(biosamp_dat$Depth_m, 1), sep = "_")
unique_site_depth <- unique(biosamp_dat$site_depth_combo)
n_site_depth_states <- length(unique_site_depth)

# Create mapping from site-depth combo to numeric index
site_depth_lookup <- data.frame(
  combo = unique_site_depth,
  index = 1:n_site_depth_states
)

# Extract site and depth for each site-depth state
site_depth_lookup$site <- as.numeric(sapply(strsplit(site_depth_lookup$combo, "_"), `[`, 1))
site_depth_lookup$depth <- as.numeric(sapply(strsplit(site_depth_lookup$combo, "_"), `[`, 2))

# Create index vector for biosamples to site-depth states
biosamp_dat$site_depth_index <- match(biosamp_dat$site_depth_combo, site_depth_lookup$combo)
biosamp_site_depth_index <- biosamp_dat$site_depth_index

# Center depths for site-depth states
site_depth_depths <- site_depth_lookup$depth - mean(site_depth_lookup$depth)

# Define the site x depth NIMBLE model (UNCHANGED)
edna_code_site_depth_occupancy <- nimbleCode({
  # Shared intercept for all sites
  intercept ~ dnorm(0, 1)
  
  # Site-depth state occupancy probabilities
  for (i in 1:n_site_depth_states) {
    logit(prob_site_depth_occurrence[i]) <- intercept + 
      b_depth_occ * site_depth_depths[i]
    site_depth_occurrence[i] ~ dbern(prob_site_depth_occurrence[i])
  }
  
  # Capture probability parameters
  cap_prob_hat ~ dnorm(0, 1.7)
  cap_prob_SD ~ dexp(1)
  
  # Site-level capture probability variation
  for (i in 1:n_sites) {
    cap_prob_logit_site[i] ~ dnorm(cap_prob_hat, cap_prob_SD)
  }
  
  # Biosample-level capture probabilities
  for (i in 1:n_biosamples) {
    logit(prob_capture[i]) <- cap_prob_logit_site[biosamp_station_index[i]] +
      b_vol * biosamp_Volume_filt_mL[i] +
      b_meth[biosamp_method_index[i]]
    
    bio_capture[i] ~ dbern(site_depth_occurrence[biosamp_site_depth_index[i]] *
                             prob_capture[i])
  }
  
  # Likelihood of detecting in each lab sample
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
  }
  
  # Priors
  b_depth_occ ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  b_meth[1] <- 0
  b_meth[2] ~ dnorm(0, 1.7)
  b_meth[3] ~ dnorm(0, 1.7)
})

# Define constants for the site x depth model
# UPDATED: Removed site_depth_sites as it was not used in the nimble model
constants_site_depth <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_site_depth_states = n_site_depth_states,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_method_index = biosamp_method_index,
  biosamp_site_depth_index = biosamp_site_depth_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_primer_index = Y_primer_index,
  site_depth_depths = site_depth_depths
)

data_site_depth <- list(
  Y = Y,
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL
)

# UPDATED: Using a function for initial values for the second model as well.
inits_fn_sitedepth <- function(){
  # Determine if a biosample had any positive detections
  max_obs_per_biosample <- tapply(Y, Y_biosamp_index, max)
  bio_capture_init <- rep(0, n_biosamples)
  bio_capture_init[as.numeric(names(max_obs_per_biosample))] <- max_obs_per_biosample
  
  # Determine if a site-depth combo had any positive biosamples
  site_depth_occurrence_init_raw <- tapply(bio_capture_init, biosamp_site_depth_index, max)
  site_depth_occurrence_init <- rep(0, n_site_depth_states)
  site_depth_occurrence_init[as.numeric(names(site_depth_occurrence_init_raw))] <- site_depth_occurrence_init_raw
  
  list(
    # Latent states (from before)
    bio_capture = bio_capture_init,
    site_depth_occurrence = site_depth_occurrence_init,
    
    # Model parameters (previously incomplete)
    intercept = 0,
    b_depth_occ = 0,
    b_vol = 0,
    cap_prob_hat = 0,
    cap_prob_SD = 1,
    cap_prob_logit_site = rep(0, n_sites),
    b_meth = c(NA, 0, 0), # NA for the first fixed element, 0 for the others
    prob_detection = rep(0.5, n_primers)
  )
}


# Run the site x depth NIMBLE model
edna_site_depth_occupancy.run <- nimbleMCMC(code = edna_code_site_depth_occupancy, 
                                            constants = constants_site_depth, 
                                            data = data_site_depth, 
                                            inits = inits_fn_sitedepth, # Using the function here
                                            niter = 200000, 
                                            nburnin = 10000, 
                                            thin = 100, 
                                            nchains = 3,
                                            summary = TRUE,
                                            samplesAsCodaMCMC = TRUE,
                                            WAIC = TRUE)
# Diagnostics
MCMCsummary(edna_site_depth_occupancy.run$samples)
mcmcplot(edna_site_depth_occupancy.run$samples, parms = c("b_depth_occ", "b_vol", "intercept"))

# --- Site x Depth Model 4-Panel Plot ---
attach.nimble(edna_site_depth_occupancy.run$samples)
n.post.new <- length(b_depth_occ)

# Depth effect on occupancy
depth_pred_seq <- seq(from = min(site_depth_lookup$depth), 
                      to = max(site_depth_lookup$depth), 
                      length.out = 100)
depth_pred_centered <- depth_pred_seq - mean(site_depth_lookup$depth)

occupancy_pred <- matrix(NA, n.post.new, length(depth_pred_seq))
for (i in 1:length(depth_pred_seq)) {
  for (j in 1:n.post.new) {
    occupancy_pred[j, i] <- inv.logit(intercept[j] + 
                                        b_depth_occ[j] * depth_pred_centered[i])
  }
}

plot_data_occupancy_depth <- data.frame(
  depth = rep(depth_pred_seq, each = n.post.new),
  occupancy = as.vector(occupancy_pred)
)

# Volume effect on capture
vol_pred_seq <- seq(from = min(as.numeric(biosamp_dat$Volume_filt_mL)), 
                    to = max(as.numeric(biosamp_dat$Volume_filt_mL)), 
                    length.out = 100)
vol_pred_centered <- vol_pred_seq - mean(as.numeric(biosamp_dat$Volume_filt_mL))

capture_pred_vol <- matrix(NA, n.post.new, length(vol_pred_seq))
for (i in 1:length(vol_pred_seq)) {
  for (j in 1:n.post.new) {
    capture_pred_vol[j, i] <- inv.logit(cap_prob_hat[j] + 
                                          b_vol[j] * vol_pred_centered[i])
  }
}

plot_data_capture_vol <- data.frame(
  volume = rep(vol_pred_seq, each = n.post.new),
  capture = as.vector(capture_pred_vol)
)

# Method comparison
gemcap_capture_new <- as.numeric(inv.logit(cap_prob_hat))
rreas_capture_new <- as.numeric(inv.logit(cap_prob_hat + b_meth[,2]))

df_methods_new <- data.frame(
  GEMCAP = gemcap_capture_new, 
  RREAS = rreas_capture_new
) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Probability")

# Detection comparison
df_detection_new <- data.frame(
  Probability = c(prob_detection[, 1], prob_detection[, 2]),
  Method = rep(c("dloop", "MiFish"), each = nrow(prob_detection))
)

# Create 4-panel plot for site x depth model
p1_new <- ggplot(plot_data_occupancy_depth, aes(x = depth, y = occupancy)) +
  stat_lineribbon(alpha = 0.25, fill = "#4CAF50", color = "#2E7D32", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth (m)", y = "Probability of Occupancy", title = "Site x Depth: Depth Effect on Occupancy") +
  theme_minimal()

p2_new <- ggplot(plot_data_capture_vol, aes(x = volume, y = capture)) +
  stat_lineribbon(alpha = 0.25, fill = "#FF9800", color = "#F57C00", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume Filtered (mL)", y = "Probability of Capture", title = "Site x Depth: Volume Effect on Capture") +
  theme_minimal()

p3_new <- ggplot(df_methods_new, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Capture Probability", y = "Density", title = "Site x Depth: RREAS vs. GEMCAP") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 1))

p4_new <- ggplot(df_detection_new, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Detection Probability", y = "Density", title = "Site x Depth: Dloop vs. MiFish") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 1))

sitedepth_4panel <- (p1_new | p2_new) / (p3_new | p4_new)
print("=== SITE x DEPTH MODEL 4-PANEL PLOT ===")
print(sitedepth_4panel)

#############
# MODEL COMPARISON
#############

# Model comparison summary
invisible({
  cat("\n=== MODEL COMPARISON SUMMARY ===\n")
  cat("Original Model:\n")
  cat("  - Single occupancy probability (prob_occurrence)\n")
  cat("  - Site as random effect on capture probability\n")
  cat("  - Depth as fixed effect on capture probability\n")
  cat("  - Volume as fixed effect on capture probability\n")
  cat("  - Method effects on capture probability\n")
  cat("  - Detection probability varies by primer type\n")
  cat("  - WAIC:", round(edna_code_vol_depth_meth_randCap.run$WAIC$WAIC, 2), "\n\n")
  
  cat("Site x Depth Model:\n")
  cat("  - Unique occupancy states for each site-depth combination\n")
  cat("  - Shared intercept for occupancy probability\n")
  cat("  - Depth as fixed effect on occupancy probability\n")
  cat("  - Site as random effect on capture probability (same as original)\n")
  cat("  - Volume as fixed effect on capture probability\n")
  cat("  - Method effects on capture probability\n")
  cat("  - Detection probability varies by primer type\n")
  cat("  - WAIC:", round(edna_site_depth_occupancy.run$WAIC$WAIC, 2), "\n")
  cat("  - Number of site-depth states:", n_site_depth_states, "\n")
  
  cat("Delta WAIC:", round(edna_site_depth_occupancy.run$WAIC$WAIC - edna_code_vol_depth_meth_randCap.run$WAIC$WAIC, 2), "\n")
})