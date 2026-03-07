library(nimble)
library(MCMCvis)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)

#############
# SIMULATION PARAMETERS (TRUE VALUES)
#############

# True parameter values
true_intercept <- 0.5
true_b_depth_occ <- -0.8
true_b_vol <- 0.6
true_cap_prob_hat <- -1.0
true_cap_prob_SD <- 0.5
true_b_meth <- c(0, 0.4, -0.3)  # GEMCAP (reference), RREAS, NCOG
true_prob_detection <- c(0.7, 0.5)  # D-loop, MiFish

# Data structure dimensions
n_sites <- 8
n_methods <- 3
n_primers <- 2
n_depths_per_site <- 3  # depths per site
n_biosamples_per_site_depth <- 4  # biological replicates per site-depth
n_lab_replicates_per_biosample <- 2  # lab replicates per biosample

#############
# GENERATE SIMULATED DATA
#############

# Create site-depth combinations
site_depth_data <- expand.grid(
  site = 1:n_sites,
  depth_raw = c(10, 50, 100)  # Three depth levels
)
n_site_depth_states <- nrow(site_depth_data)

# Center depths
site_depth_data$depth_centered <- site_depth_data$depth_raw - mean(site_depth_data$depth_raw)

# Generate true occupancy states for each site-depth combo
site_depth_data$prob_occurrence <- plogis(true_intercept + 
                                            true_b_depth_occ * site_depth_data$depth_centered)
site_depth_data$site_depth_occurrence <- rbinom(n_site_depth_states, 1, 
                                                site_depth_data$prob_occurrence)

# Generate site-level capture probability random effects
true_cap_prob_logit_site <- rnorm(n_sites, true_cap_prob_hat, true_cap_prob_SD)

# Create biosample data
biosample_data <- data.frame()
biosample_id <- 1

for (i in 1:n_site_depth_states) {
  site <- site_depth_data$site[i]
  depth <- site_depth_data$depth_raw[i]
  depth_centered <- site_depth_data$depth_centered[i]
  site_depth_state <- site_depth_data$site_depth_occurrence[i]
  
  for (b in 1:n_biosamples_per_site_depth) {
    # Random method and volume for this biosample
    method <- sample(1:n_methods, 1)
    volume_raw <- runif(1, 500, 2000)  # mL
    volume_centered <- volume_raw - mean(c(500, 2000))
    
    # Calculate capture probability
    logit_cap_prob <- true_cap_prob_logit_site[site] + 
      true_b_vol * volume_centered + 
      true_b_meth[method]
    prob_capture <- plogis(logit_cap_prob)
    
    # Generate bio_capture (conditional on occupancy)
    bio_capture <- rbinom(1, 1, site_depth_state * prob_capture)
    
    biosample_data <- rbind(biosample_data, data.frame(
      biosample_id = biosample_id,
      site = site,
      depth_raw = depth,
      depth_centered = depth_centered,
      site_depth_index = i,
      method = method,
      volume_raw = volume_raw,
      volume_centered = volume_centered,
      bio_capture = bio_capture
    ))
    
    biosample_id <- biosample_id + 1
  }
}

n_biosamples <- nrow(biosample_data)

# Generate lab replicate (detection) data
lab_data <- data.frame()

for (i in 1:n_biosamples) {
  bio_cap <- biosample_data$bio_capture[i]
  biosamp_id <- biosample_data$biosample_id[i]
  
  for (rep in 1:n_lab_replicates_per_biosample) {
    primer <- sample(1:n_primers, 1)
    prob_detect <- true_prob_detection[primer]
    
    # Detection conditional on bio_capture
    Y <- rbinom(1, 1, bio_cap * prob_detect)
    
    lab_data <- rbind(lab_data, data.frame(
      biosample_id = biosamp_id,
      primer = primer,
      Y = Y
    ))
  }
}

N <- nrow(lab_data)

# Recenter volume for all biosamples
mean_volume <- mean(biosample_data$volume_raw)
biosample_data$volume_centered <- biosample_data$volume_raw - mean_volume

#############
# PREPARE DATA FOR NIMBLE
#############

constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_site_depth_states = n_site_depth_states,
  N = N,
  biosamp_station_index = biosample_data$site,
  biosamp_method_index = biosample_data$method,
  biosamp_site_depth_index = biosample_data$site_depth_index,
  Y_biosamp_index = lab_data$biosample_id,
  Y_primer_index = lab_data$primer,
  site_depth_depths = site_depth_data$depth_centered
)

data <- list(
  Y = lab_data$Y,
  biosamp_Volume_filt_mL = biosample_data$volume_centered
)

# NIMBLE model code (same as original)
edna_code_site_depth_occupancy <- nimbleCode({
  intercept ~ dnorm(0, 1)
  
  for (i in 1:n_site_depth_states) {
    logit(prob_site_depth_occurrence[i]) <- intercept + 
      b_depth_occ * site_depth_depths[i]
    site_depth_occurrence[i] ~ dbern(prob_site_depth_occurrence[i])
  }
  
  cap_prob_hat ~ dnorm(0, 1.7)
  cap_prob_SD ~ dexp(1)
  
  for (i in 1:n_sites) {
    cap_prob_logit_site[i] ~ dnorm(cap_prob_hat, cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    logit(prob_capture[i]) <- cap_prob_logit_site[biosamp_station_index[i]] +
      b_vol * biosamp_Volume_filt_mL[i] +
      b_meth[biosamp_method_index[i]]
    
    bio_capture[i] ~ dbern(site_depth_occurrence[biosamp_site_depth_index[i]] *
                             prob_capture[i])
  }
  
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
  }
  
  b_depth_occ ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  b_meth[1] <- 0
  b_meth[2] ~ dnorm(0, 1.7)
  b_meth[3] ~ dnorm(0, 1.7)
})

# Initial values function
inits_fn <- function(){
  max_obs_per_biosample <- tapply(data$Y, constants$Y_biosamp_index, max)
  bio_capture_init <- rep(0, constants$n_biosamples)
  bio_capture_init[as.numeric(names(max_obs_per_biosample))] <- max_obs_per_biosample
  
  site_depth_occurrence_init_raw <- tapply(bio_capture_init, 
                                           constants$biosamp_site_depth_index, 
                                           max)
  site_depth_occurrence_init <- rep(0, constants$n_site_depth_states)
  site_depth_occurrence_init[as.numeric(names(site_depth_occurrence_init_raw))] <- site_depth_occurrence_init_raw
  
  list(
    bio_capture = bio_capture_init,
    site_depth_occurrence = site_depth_occurrence_init,
    intercept = 0,
    b_depth_occ = 0,
    b_vol = 0,
    cap_prob_hat = 0,
    cap_prob_SD = 1,
    cap_prob_logit_site = rep(0, constants$n_sites),
    b_meth = c(NA, 0, 0),
    prob_detection = rep(0.5, n_primers)
  )
}

#############
# FIT MODEL TO SIMULATED DATA
#############

cat("Fitting model to simulated data...\n")
cat("This may take several minutes...\n\n")

model_run <- nimbleMCMC(
  code = edna_code_site_depth_occupancy, 
  constants = constants, 
  data = data, 
  inits = inits_fn,
  niter = 50000,
  nburnin = 5000, 
  thin = 50, 
  nchains = 3,
  summary = TRUE,
  samplesAsCodaMCMC = TRUE,
  WAIC = FALSE
)

#############
# PARAMETER RECOVERY ASSESSMENT
#############

# Extract posterior summaries
post_summary <- MCMCsummary(model_run$samples, round = 3)

# Key parameters to check
params_to_check <- data.frame(
  Parameter = c("intercept", "b_depth_occ", "b_vol", "cap_prob_hat", 
                "cap_prob_SD", "b_meth[2]", "b_meth[3]", 
                "prob_detection[1]", "prob_detection[2]"),
  True_Value = c(true_intercept, true_b_depth_occ, true_b_vol, 
                 true_cap_prob_hat, true_cap_prob_SD, 
                 true_b_meth[2], true_b_meth[3],
                 true_prob_detection[1], true_prob_detection[2])
)

# Add posterior estimates
params_to_check$Posterior_Mean <- post_summary[params_to_check$Parameter, "mean"]
params_to_check$Posterior_SD <- post_summary[params_to_check$Parameter, "sd"]
params_to_check$CI_Lower <- post_summary[params_to_check$Parameter, "2.5%"]
params_to_check$CI_Upper <- post_summary[params_to_check$Parameter, "97.5%"]
params_to_check$Rhat <- post_summary[params_to_check$Parameter, "Rhat"]

# Check if true value is in 95% CI
params_to_check$In_CI <- (params_to_check$True_Value >= params_to_check$CI_Lower & 
                            params_to_check$True_Value <= params_to_check$CI_Upper)

# Calculate bias and relative bias
params_to_check$Bias <- params_to_check$Posterior_Mean - params_to_check$True_Value
params_to_check$Rel_Bias <- params_to_check$Bias / abs(params_to_check$True_Value)

# Print results
cat("\n=== PARAMETER RECOVERY RESULTS ===\n\n")
print(params_to_check, row.names = FALSE)

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
cat("All Rhat < 1.1:", all(params_to_check$Rhat < 1.1), "\n")
cat("Parameters with true value in 95% CI:", sum(params_to_check$In_CI), 
    "out of", nrow(params_to_check), "\n")

#############
# VISUALIZATION
#############

# Create recovery plot with parameters on x-axis
recovery_plot <- ggplot(params_to_check, aes(x = Parameter)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = Posterior_Mean, color = "Estimated"), size = 3) +
  geom_point(aes(y = True_Value, color = "True"), size = 3, shape = 17) +
  scale_color_manual(values = c("Estimated" = "darkblue", "True" = "red"),
                     name = "") +
  labs(x = "Parameter", 
       y = "Parameter Value",
       title = "Parameter Recovery: True Values vs Posterior Estimates",
       subtitle = "Error bars show 95% credible intervals") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.major.x = element_blank())

print(recovery_plot)

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Mean absolute bias:", mean(abs(params_to_check$Bias)), "\n")
cat("Mean relative bias:", mean(abs(params_to_check$Rel_Bias), na.rm = TRUE), "\n")
cat("Max |relative bias|:", max(abs(params_to_check$Rel_Bias), na.rm = TRUE), "\n")

# Data summary
cat("\n=== SIMULATED DATA SUMMARY ===\n")
cat("Total lab samples:", N, "\n")
cat("Total biosamples:", n_biosamples, "\n")
cat("Site-depth states:", n_site_depth_states, "\n")
cat("Detection rate:", mean(lab_data$Y), "\n")
cat("Biosample capture rate:", mean(biosample_data$bio_capture), "\n")
cat("Site-depth occupancy rate:", mean(site_depth_data$site_depth_occurrence), "\n")
