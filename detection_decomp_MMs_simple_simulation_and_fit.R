library(nimble)

simulate_edna_data <- function(n_sites = 80, prob_occurrence = 0.6, prob_rreas_capture = 0.7, prob_gemcap_capture = 0.5, prob_rreas_detection = 0.9, prob_gemcap_detection = 0.8) {
  # Initialize data frame to store results
  results <- data.frame(
    site = integer(),
    method = character(),
    biological_replicate = integer(),
    technical_replicate = integer(),
    dna_presence = integer()
  )
  
  # Loop through each site
  for (site in 1:n_sites) {
    # Simulate the true state of eDNA occurrence at the site
    site_occurrence <- rbinom(1, 1, prob_occurrence)
    
    # RREAS Method (3 biological replicates)
    for (bio_rep in 1:3) {
      # Simulate if mammal DNA is captured in biological replicate
      bio_capture <- ifelse(site_occurrence == 1, rbinom(1, 1, prob_rreas_capture), 0)
      
      # Loop through technical replicates
      for (tech_rep in 1:3) {
        # Simulate detection in technical replicate if biological replicate is positive
        tech_detection <- ifelse(bio_capture == 1, rbinom(1, 1, prob_rreas_detection), 0)
        
        # Append result to data frame
        results <- rbind(results, data.frame(
          site = site,
          method = "RREAS",
          biological_replicate = bio_rep,
          technical_replicate = tech_rep,
          dna_presence = tech_detection
        ))
      }
    }
    
    # GEMCAP Method (1 biological replicate)
    # Simulate if mammal DNA is captured in biological replicate
    bio_capture <- ifelse(site_occurrence == 1, rbinom(1, 1, prob_gemcap_capture), 0)
    
    # Loop through technical replicates
    for (tech_rep in 1:3) {
      # Simulate detection in technical replicate if biological replicate is positive
      tech_detection <- ifelse(bio_capture == 1, rbinom(1, 1, prob_gemcap_detection), 0)
      
      # Append result to data frame
      results <- rbind(results, data.frame(
        site = site,
        method = "GEMCAP",
        biological_replicate = 1,
        technical_replicate = tech_rep,
        dna_presence = tech_detection
      ))
    }
  }
  
  return(results)
}

# Run the simulation
#set.seed(123)  # For reproducibility
simulated_data <- simulate_edna_data()
head(simulated_data)

library(nimble)

# Define the NIMBLE model
edna_code <- nimbleCode({
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    
    for (m in 1:2) { # Loop over methods (RREAS and GEMCAP)
      for (j in 1:n_bio_reps[i, m]) {
        # Biological replicate capture probability (method dependent)
        bio_capture[i, m, j] ~ dbern(site_occurrence[i] * prob_capture[m])
      }
    }
  }
  
  # Likelihood
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[site_index[i], method_index[i], bio_index[i]] * prob_detection[method_index[i]])
  }
  
  # Priors for the parameters
  prob_occurrence ~ dbeta(1, 1)
  prob_capture[1] ~ dbeta(1, 1) # RREAS capture probability
  prob_capture[2] ~ dbeta(1, 1) # GEMCAP capture probability
  prob_detection[1] ~ dbeta(1, 1) # RREAS detection probability
  prob_detection[2] ~ dbeta(1, 1) # GEMCAP detection probability
})

# Constants for the model
n_sites <- length(unique(simulated_data$site))
n_bio_reps <- matrix(c(rep(3, n_sites), rep(1, n_sites)), nrow = n_sites, ncol = 2)  # 3 biological replicates for RREAS, 1 for GEMCAP
n_tech_reps <- array(3, dim = c(n_sites, max(n_bio_reps)))  # 3 technical replicates for each biological replicate
method <- matrix(c(rep(1, n_sites * 2), rep(2, n_sites)), nrow = n_sites, byrow = TRUE)  # 1 for RREAS, 2 for GEMCAP

# Create indices for likelihood
site_index <- simulated_data$site
method_index <- ifelse(simulated_data$method == "RREAS", 1, 2)
bio_index <- simulated_data$biological_replicate
N <- nrow(simulated_data)
Y <- simulated_data$dna_presence

# Data for the model
data <- list(n_sites = n_sites,
             n_bio_reps = n_bio_reps,
             n_tech_reps = n_tech_reps,
             method = method,
             Y = Y,
             site_index = site_index,
             method_index = method_index,
             bio_index = bio_index)

# Initial values for MCMC
inits <- list(prob_occurrence = 0.5,
              prob_capture = c(0.7, 0.5),
              prob_detection = c(0.9, 0.8))

# Create the NIMBLE model
edna_model <- nimbleModel(code = edna_code, constants = data, inits = inits)

# Compile the model
compiled_edna_model <- compileNimble(edna_model)

# Configure and run MCMC
mcmc_conf <- configureMCMC(edna_model)
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = edna_model)

# Run the MCMC
mcmc_samples <- runMCMC(compiled_mcmc, niter = 10000)

# Summary of posterior samples
summary(mcmc_samples)


# cool, seems to work great. Now let's run it on our actual data. 