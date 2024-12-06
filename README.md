
---
Marine mammal eDNA site occupancy, capture and detection probabilities estimated from an intensively replicated ship-based survey
---



# Abstract

Environmental DNA (eDNA) sampling holds great potential as a non-invasive method for species detection and marine biodiversity monitoring. However, eDNA monitoring strategies are limited by a paucity of quantitative studies that estimate detection probabilities for target taxa aross the various levels of study design (e.g. sample sites, biological replicates, technical replicates).  In this study, we appled multi-level occupancy modeling techniques to data from a ship-based marine mammal eDNA survey. Over a 8 day cruise, we surveyed 16 sites by collecting and filtering multiple Niskin bottle water samples (biological replicates) of various volumes and at various depths (including replicate biological samples at each depth) at each site. Each filter was subsequenlty processed using a Dloop metabarcoding laboratory workflow, including 3 technical replicates per sample, to detect marine mammal taxa. Using the resulting sequence data, our quantitative framework estimates 1) eDNA detection probabilities associated with technical replication, 2) capture probabilities associated with biological sample replication, and as a function of water volume filtered and sample depth, and 3) the probability of marine mammal eDNA site occupancy. We found a weak positive effect of water volume and weakly negative effect of sample depth on biological sample capture probability. On average, both per-sample capture/replicate detection probabilities were below 50%. Conversely, the model-estimated occurrence of marine mammal eDNA at any given site approached 100%. The near-ubiquitous presence of detectable eDNA across sampling sites indicates great potential for eDNA monitoring applications once survey design and replication are appropriately tuned.


# Introduction

Environmental DNA (eDNA) is an emerging tool for monitoring marine biodiversity, offering a non-invasive and cost-effective approach to detect and identify species in aquatic ecosystems [@thomsen2015; @deiner2017]. This method, which analyzes genetic material shed by organisms into their environment, has shown particular promise in the field of marine mammal research, where traditional monitoring techniques often face limitations due to the elusive nature of these animals and the vast expanse of their habitats [@foote2012].

In recent years, the application of eDNA techniques has expanded rapidly, with studies demonstrating its effectiveness in detecting rare, elusive, and threatened marine mammal species [@baker2018; @parsons2018]. However, while the potential of eDNA for marine mammal monitoring has been established, quantitative estimates of occurrence rates and capture and detection probabilities have remained largely unexplored. This gap in knowledge has limited the ability of researchers to fully leverage eDNA as a tool for marine mammal population assessment and management [@goldberg2016].

Method details and approach here. 


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

## Results

Plots from modeling code etc. go here. Plots need cleanup and fanciness first. 

## Discussion

Our study presents a novel approach to addressing this knowledge gap by employing advanced occupancy modeling techniques to estimate occurrence rates and detection probabilities associated with an intensive, ship-based marine mammal eDNA sampling survey [@mackenzie2002; @schmidt2013]. The results of our analysis reveal a striking finding: the model-estimated occurrence of marine mammal eDNA at any given site approaches 100%. This remarkable outcome suggests that eDNA sampling has the potential to detect even the rarest of marine organisms, offering unprecedented sensitivity in biodiversity monitoring [@ficetola2015].

By providing quantitative estimates of occurrence and detection probabilities, our study not only enhances the reliability of eDNA-based marine mammal monitoring but also establishes a strong case for its widespread adoption in marine conservation efforts [@bohmann2014]. The near-ubiquitous presence of detectable eDNA signals across sampling sites indicates that this method could revolutionize our ability to track and protect marine mammal populations, including those that are critically endangered or traditionally difficult to observe [@sigsgaard2017].

Furthermore, our research demonstrates the power of integrating eDNA analysis with sophisticated statistical modeling, paving the way for more comprehensive and accurate assessments of marine ecosystems [@lacoursiere-roussel2016]. This approach addresses key challenges in marine mammal monitoring, such as imperfect detection and the need for cost-effective, large-scale surveillance methods [@chambert2018].

eDNA monitoring represents a transformative tool in marine conservation, capable of providing unprecedented insights into the distribution and abundance of marine mammals [@taberlet2018]. The quantitative framework presented in this paper offers a robust foundation for future research and management strategies, potentially reshaping our approach to marine biodiversity assessment and conservation planning [@cristescu2018].

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
