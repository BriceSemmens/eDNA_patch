
---
Marine Mammal eDNA patch occupancy modeling 
---


## Modeling Approach for Marine Mammal eDNA Detection

Our eDNA monitoring involved 16 fixed sampling sites, each surveyed using two distinct water collection methods: RREAS, which collected three separate biological replicates per site, and GEMCAP, which collected a single biological replicate per site. Each biological replicate was then subsampled into three technical replicates analyzed for the presence of marine mammal DNA. Thus, the data structure included three nested levels: (1) site-level occurrence of marine mammal eDNA, (2) method- and site-specific probabilities of capturing eDNA within a given biological replicate (conditional on presence), and (3) technical replicate-level probabilities of detecting eDNA (conditional on capture).

We developed a hierarchical Bayesian model to account for the complex sampling design of our environmental DNA (eDNA) marine mammal monitoring program. The model explicitly accounts for multiple sources of variation and detection uncertainty across biological and technical replicates, sampling methods, and sites.

Data Simulation for Model Validation
Prior to fitting the model to empirical data, we developed a data simulation procedure to validate model performance and verify parameter identifiability. The simulation code, written in R, generated synthetic data consistent with our sampling design and hypothesized data-generating process. Specifically, we first specified a known occurrence probability for marine mammal eDNA across sites. For each site, we then simulated a site-level occurrence state. Given occurrence, we simulated the biological replicate capture process by assigning capture probabilities on the logit scale, incorporating fixed effects of volume filtered and depth of sampling (both standardized), as well as a site-specific random effect and a fixed effect for the sampling method. Finally, given capture within a biological replicate, we simulated the detection process at the technical replicate level using a single detection probability parameter.


### Hierarchical Model Structure

The model can be represented by the following hierarchical probability structure:

1. **Site-Level Occurrence**
   The presence of marine mammals at site $s$ is modeled as a Bernoulli random variable:

   $Z_s \sim \text{Bernoulli}(\psi)$

   where $Z_s$ is the site-level occurrence indicator, and $\psi$ is the overall occurrence probability, drawn from a Beta prior:

   $\psi \sim \text{Beta}(1,1)$

2. **Biological Replicate Capture**
   The logit-linear model for biological replicate capture probability is:

   $\text{logit}(p_{\text{capture},b}) = \beta_0 + \beta_{\text{vol}} \cdot X_{\text{vol},b} + \beta_{\text{depth}} \cdot X_{\text{depth},b} + \gamma_{s[b]} + \delta_{\text{method}[b]}$

   Where:
   - $p_{\text{capture},b}$ is the capture probability for biological replicate $b$
   - $\beta_0$ is the intercept (site-level capture probability hyperparameter)
   - $\beta_{\text{vol}}$ is the volume coefficient
   - $\beta_{\text{depth}}$ is the depth coefficient
   - $X_{\text{vol},b}$ is the centered water volume
   - $X_{\text{depth},b}$ is the centered sampling depth
   - $\gamma_{s[b]}$ is the site-specific random effect
   - $\delta_{\text{method}[b]}$ is the method-specific fixed effect, with method effects constrained such that:
   $\delta_{\text{GEMCAP}} = 0$ and 
   $\delta_{\text{RREAS}} \sim \text{Normal}(0, 1.7)$
   
   
    The $\delta_{\text{method}[b]}$ parameterization ensures that the GEMCAP method serves as the reference level, with the RREAS method's effect estimated relative to GEMCAP.

   The biological replicate capture is then modeled as:

   $Y_{\text{capture},b} \sim \text{Bernoulli}(Z_{s[b]} \cdot p_{\text{capture},b})$

3. **Technical Replicate Detection**
   Conditional on biological replicate capture, technical replicates are modeled as:

   $Y_{\text{detect},i} \sim \text{Bernoulli}(p_{\text{detect}} \cdot Y_{\text{capture},b[i]})$

   where $p_{\text{detect}}$ is the detection probability, drawn from a Beta prior:

   $p_{\text{detect}} \sim \text{Beta}(1,1)$

### Prior Distributions

We specified weakly informative priors for model parameters:

- Volume and depth coefficients: $\beta_{\text{vol}}, \beta_{\text{depth}} \sim \text{Normal}(0, 1.7)$
- Method effects: $\delta_{\text{method}} \sim \text{Normal}(0, 1.7)$
- Site-level random effects: $\gamma_s \sim \text{Normal}(\beta_0, \sigma_{\text{cap}})$
- Site-level capture probability hyperparameter: $\beta_0 \sim \text{Normal}(0, 1.7)$
- Site-level random effect standard deviation: $\sigma_{\text{cap}} \sim \text{Exponential}(1)$

### Inference

We used Markov Chain Monte Carlo (MCMC) sampling to estimate posterior distributions of model parameters:
- 3 parallel chains 
- 20,000 total iterations
- 10,000 iterations discarded as burn-in
- Thinning interval of 10

Model performance and parameter estimates were assessed using posterior summaries, trace plots, and the Widely Applicable Information Criterion (WAIC).

