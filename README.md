---
Marine mammal eDNA site occupancy, capture and detection probabilities estimated from an intensively replicated ship-based survey
---

# Abstract

Environmental DNA (eDNA) sampling holds great potential as a non-invasive method for species detection and marine biodiversity monitoring. However, eDNA monitoring strategies are limited by a paucity of quantitative studies that estimate detection probabilities for target taxa across the various levels of study design (e.g. sample sites, biological replicates, technical replicates). In this study, we applied multi-level occupancy modeling techniques to data from a ship-based marine mammal eDNA survey. Over an 8 day cruise, we surveyed 16 sites by collecting and filtering multiple Niskin bottle water samples (biological replicates) of various volumes and at various depths (including replicate biological samples at each depth) at each site. Each filter was subsequently processed using both D-loop and MiFish metabarcoding laboratory workflows, including multiple technical replicates per sample, to detect marine mammal taxa. Using the resulting sequence data, our quantitative framework compares two alternative model formulations to estimate 1) eDNA detection probabilities associated with technical replication and primer choice, 2) capture probabilities associated with biological sample replication as a function of water volume filtered, sample depth, and collection method, and 3) the probability of marine mammal eDNA site occupancy. We found weak positive effects of water volume and collection method on biological sample capture probability. The analysis compares models where depth affects either capture probability or occupancy probability directly, providing insights into the ecological interpretation of depth effects in eDNA sampling. Model comparison using WAIC indicates whether site×depth-specific occupancy states better explain the data than simpler site-level occupancy. On average, both per-sample capture and detection probabilities were below 50%, while model-estimated occurrence of marine mammal eDNA approached 100% across sampling sites. The near-ubiquitous presence of detectable eDNA across sampling sites indicates great potential for eDNA monitoring applications once survey design and replication are appropriately tuned.

# Repository Structure

The main analysis is contained in:
- **`detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo.R`**: Primary analysis script implementing both model formulations with clean NIMBLE code and comprehensive model comparison

Supporting files:
- **`./Archive/detection_decomp_MMs_real_data_with_miFish_and_NCOG.R`**: Earlier version with original model formulation (retained for reference)
- **`./Data/intercal_ALL_metadata_12.30.24-noCs_NEW.csv`**: Input data file containing eDNA detection results and sampling metadata

# Introduction

Environmental DNA (eDNA) is an emerging tool for monitoring marine biodiversity, offering a non-invasive and cost-effective approach to detect and identify species in aquatic ecosystems [@thomsen2015; @deiner2017]. This method, which analyzes genetic material shed by organisms into their environment, has shown particular promise in the field of marine mammal research, where traditional monitoring techniques often face limitations due to the elusive nature of these animals and the vast expanse of their habitats [@foote2012].

In recent years, the application of eDNA techniques has expanded rapidly, with studies demonstrating its effectiveness in detecting rare, elusive, and threatened marine mammal species [@baker2018; @parsons2018]. However, while the potential of eDNA for marine mammal monitoring has been established, quantitative estimates of occurrence rates and capture and detection probabilities have remained largely unexplored. This gap in knowledge has limited the ability of researchers to fully leverage eDNA as a tool for marine mammal population assessment and management [@goldberg2016].

## Modeling Approach for Marine Mammal eDNA Detection

Our eDNA monitoring involved 16 fixed sampling sites, each surveyed using multiple water collection methods (GEMCAP, RREAS) with varying numbers of biological replicates per site and method. Each biological replicate was collected at specific depths and volumes, then processed using two primer systems (D-loop and MiFish) with multiple technical replicates analyzed for the presence of marine mammal DNA. This nested sampling design creates a complex hierarchical structure requiring specialized analytical approaches.

We developed two alternative hierarchical Bayesian model formulations implemented in NIMBLE to account for the complex sampling design and compare different ecological hypotheses about how environmental factors affect eDNA detection:

### Model 1: Original Single Occupancy Model

This model treats each site as having a single occupancy state and models depth as affecting capture probability:

**Ecological Interpretation**: Marine mammals either occur or don't occur at a site (binary site-level occupancy). Depth affects how likely you are to detect eDNA given the animals are present (e.g., deeper water might make eDNA harder to capture due to dilution effects).

**Model Structure**:
1. **Site-Level Occurrence**: Single global occurrence probability for all sites
   ```
   site_occurrence[s] ~ Bernoulli(prob_occurrence)
   ```

2. **Biological Replicate Capture**: Depth affects capture probability
   ```
   logit(prob_capture[b]) = cap_prob_logit[site] + b_depth × depth + b_vol × volume + b_meth[method]
   bio_capture[b] ~ Bernoulli(site_occurrence[site] × prob_capture[b])
   ```

3. **Technical Replicate Detection**: Primer-specific detection probabilities
   ```
   Y[i] ~ Bernoulli(bio_capture[biosample] × prob_detection[primer])
   ```

### Model 2: Site×Depth Occupancy Model

This model creates unique occupancy states for each site-depth combination and models depth as affecting occupancy probability:

**Ecological Interpretation**: Marine mammals have different probabilities of being present at different depths within the same site. Depth affects the fundamental occupancy probability, recognizing that animals may prefer certain depth ranges and be absent from others.

**Model Structure**:
1. **Site×Depth-Level Occurrence**: Depth affects occupancy probability
   ```
   logit(prob_site_depth_occurrence[i]) = intercept + b_depth_occ × site_depth_depths[i]
   site_depth_occurrence[i] ~ Bernoulli(prob_site_depth_occurrence[i])
   ```

2. **Biological Replicate Capture**: No depth effect (depth already modeled in occupancy)
   ```
   logit(prob_capture[b]) = cap_prob_logit_site[site] + b_vol × volume + b_meth[method]
   bio_capture[b] ~ Bernoulli(site_depth_occurrence[site_depth_combo] × prob_capture[b])
   ```

3. **Technical Replicate Detection**: Same as Model 1
   ```
   Y[i] ~ Bernoulli(bio_capture[biosample] × prob_detection[primer])
   ```

### Data Simulation for Model Validation

Prior to fitting models to empirical data, we developed data simulation procedures to validate model performance and verify parameter identifiability. The simulation code generates synthetic data consistent with our sampling design and hypothesized data-generating processes, allowing us to test whether our models can recover known parameter values.

### Model Implementation and Comparison

Both models are implemented using clean NIMBLE code with proper initialization functions to eliminate warnings and ensure reliable convergence. Key improvements in the main analysis script include:

- **Robust initialization**: Functions that ensure latent states are consistent with observed data
- **Clean constants**: Removal of unused variables that caused NIMBLE warnings  
- **Complete parameterization**: All model parameters properly initialized and constrained
- **Model comparison**: WAIC-based comparison to determine which model formulation better fits the data

### Prior Distributions

We specified weakly informative priors for model parameters:

- Volume and depth coefficients: `Normal(0, 1)`
- Method effects: `Normal(0, 1.7)` with GEMCAP as reference (fixed at 0)
- Site-level random effects: `Normal(cap_prob_hat, cap_prob_SD)`
- Capture probability hyperparameters: `cap_prob_hat ~ Normal(0, 1.7)`, `cap_prob_SD ~ Exponential(1)`
- Occurrence probability: `Beta(1, 1)`
- Detection probabilities: `Beta(1, 1)` for each primer type

### Inference

MCMC sampling parameters:
- 3 parallel chains 
- 200,000 total iterations
- 10,000 iterations discarded as burn-in
- Thinning interval of 100
- Final posterior samples: 5,700 per chain (17,100 total)

Model performance assessed using:
- Posterior summaries and convergence diagnostics (R̂, effective sample size)
- Trace plots and autocorrelation diagnostics
- WAIC for model comparison
- Posterior predictive checks

### Key Research Questions

1. **Method Effects**: How do different collection methods (GEMCAP vs. RREAS) affect eDNA capture probability?

2. **Volume Effects**: Does filtering larger volumes of water increase capture probability?

3. **Depth Effects**: Does sample depth affect eDNA detectability, and is this effect better modeled as affecting occupancy or capture probability?

4. **Primer Performance**: How do D-loop and MiFish primers compare in their detection probabilities?

5. **Model Selection**: Which model structure (single occupancy vs. site×depth occupancy) better explains the observed eDNA detection patterns?

## Results

The analysis provides quantitative estimates of detection probabilities across all levels of the sampling hierarchy, with model comparison indicating whether environmental stratification (depth) is better incorporated at the occupancy or capture level. Detailed results including posterior distributions, effect sizes, and model comparison metrics are generated by the main analysis script.

## Discussion

Our study presents a novel approach to addressing knowledge gaps in eDNA monitoring by employing advanced occupancy modeling techniques to estimate occurrence rates and detection probabilities from an intensive, ship-based marine mammal eDNA sampling survey [@mackenzie2002; @schmidt2013]. The comparison of two alternative model formulations provides insights into the ecological mechanisms underlying eDNA detection patterns and optimal sampling strategies.

By providing quantitative estimates of occurrence and detection probabilities across multiple levels of sampling hierarchy, our study enhances the reliability of eDNA-based marine mammal monitoring and establishes frameworks for optimizing survey design [@bohmann2014]. The model comparison approach demonstrates how different ecological hypotheses can be formally tested using hierarchical modeling, advancing both the statistical methodology and biological understanding of eDNA detection processes.

Furthermore, our research demonstrates the power of integrating eDNA analysis with sophisticated statistical modeling, paving the way for more comprehensive and accurate assessments of marine ecosystems [@lacoursiere-roussel2016]. This approach addresses key challenges in marine mammal monitoring, such as imperfect detection and the need for cost-effective, large-scale surveillance methods [@chambert2018].

## References

@thomsen2015: Thomsen, P. F., & Willerslev, E. (2015). Environmental DNA – An emerging tool in conservation for monitoring past and present biodiversity. Biological Conservation, 183, 4-18.

@deiner2017: Deiner, K., Bik, H. M., Mächler, E., Seymour, M., Lacoursière-Roussel, A., Altermatt, F., ... & Bernatchez, L. (2017). Environmental DNA metabarcoding: Transforming how we survey animal and plant communities. Molecular Ecology, 26(21), 5872-5895.

@foote2012: Foote, A. D., Thomsen, P. F., Sveegaard, S., Wahlberg, M., Kielgast, J., Kyhn, L. A., ... & Gilbert, M. T. P. (2012). Investigating the potential use of environmental DNA (eDNA) for genetic monitoring of marine mammals. PloS one, 7(8), e41781.

@baker2018: Baker, C. S., Steel, D., Nieukirk, S., & Klinck, H. (2018). Environmental DNA (eDNA) from the wake of the whales: droplet digital PCR for detection and species identification. Frontiers in Marine Science, 5, 133.

@parsons2018: Parsons, K. M., Everett, M., Dahlheim, M., & Park, L. (2018). Water, water everywhere: environmental DNA can unlock population structure in elusive marine species. Royal Society Open Science, 5(8), 180537.

@goldberg2016: Goldberg, C. S., Turner, C. R., Deiner, K., Klymus, K. E., Thomsen, P. F., Murphy, M. A., ... & Taberlet, P. (2016). Critical considerations for the application of environmental DNA methods to detect aquatic species. Methods in Ecology and Evolution, 7(11), 1299-1307.

@mackenzie2002: MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Andrew Royle, J., & Langtimm, C. A. (2002). Estimating site occupancy rates when detection probabilities are less than one. Ecology, 83(8), 2248-2255.

@schmidt2013: Schmidt, B. R., Kéry, M., Ursenbacher, S., Hyman, O. J., & Collins, J. P. (2013). Site occupancy models in the analysis of environmental DNA presence/absence surveys: a case study of an emerging amphibian pathogen. Methods in Ecology and Evolution, 4(7), 646-653.

@ficetola2015: Ficetola, G. F., Pansu, J., Bonin, A., Coissac, E., Giguet-Covex, C., De Barba, M., ... & Taberlet, P. (2015). Replication levels, false presences and the estimation of the presence/absence from eDNA metabarcoding data. Molecular Ecology Resources, 15(3), 543-556.

@bohmann2014: Bohmann, K., Evans, A., Gilbert, M. T. P., Carvalho, G. R., Creer, S., Knapp, M., ... & de Bruyn, M. (2014). Environmental DNA for wildlife biology and biodiversity monitoring. Trends in Ecology & Evolution, 29(6), 358-367.

@sigsgaard2017: Sigsgaard, E. E., Nielsen, I. B., Bach, S. S., Lorenzen, E. D., Robinson, D. P., Knudsen, S. W., ... & Thomsen, P. F. (2017). Population characteristics of a large whale shark aggregation inferred from seawater environmental DNA. Nature Ecology & Evolution, 1(1), 0004.

@lacoursiere-roussel2016: Lacoursière-Roussel, A., Côté, G., Leclerc, V., & Bernatchez, L. (2016). Quantifying relative fish abundance with eDNA: a promising tool for fisheries management. Journal of Applied Ecology, 53(4), 1148-1157.

@chambert2018: Chambert, T., Pilliod, D. S., Goldberg, C. S., Doi, H., & Takahara, T. (2018). An analytical framework for estimating aquatic species density from environmental DNA. Ecology and Evolution, 8(6), 3468-3477.

@taberlet2018: Taberlet, P., Bonin, A., Zinger, L., & Coissac, E. (2018). Environmental DNA: For biodiversity research and monitoring. Oxford University Press.

@cristescu2018: Cristescu, M. E., & Hebert, P. D. (2018). Uses and misuses of environmental DNA in biodiversity science and conservation. Annual Review of Ecology, Evolution, and Systematics, 49, 209-230.
