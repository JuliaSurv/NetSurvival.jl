```@meta
CurrentModule = NetSurvival
```

# Theory

## General Notations

We will denote $S_O(t)$ as the overall survival probability, representing the probability that an individual is still alive regardless of the cause of death. Another key notation is $S_P(t)$ indicating the estimated probability of survival based on population data for individuals with similar attributes as the patient, often obtained from mortality rate tables. The relative survival ratio is thus defined as $S_R(t) = \frac{S_O(t)}{S_P(t)}$. This value allows comparison of the patients' survival probability to that of the general public.

We will define the overall hazard $\lambda_{O_i}$ for each individual $i$ as the sum between the hazard relative to the general population $\lambda_{P_i}$ and the one relative to the disease $\lambda_{E_i}$. Given our focus on the disease-specific hazard, we define the individual relative survival ratio as:

$$S_{E_i}(t) = exp \{-\int_0^t \lambda_{E_i}(u)du\} = \frac{exp \{-\int_0^t \lambda_{O_i}(u)du\}}{exp \{-\int_0^t \lambda_{P_i}(u)du\}}$$

Thus, the net survival of an $N$ sized cohort would be:

$$S_E(t) = \frac{1}{N} \sum_{i=1}^N S_{E_i}(t)$$

Using these notations, we can now easily define the net survival ratio as:

$$S_R(t) = \frac{\frac{1}{N}\sum_i S_{O_i}(t)}{\frac{1}{N}\sum_i S_{P_i}(t)}$$

$$S_E(t) =\frac{1}{N} \sum_i \frac{S_{O_i}(t)}{S_{P_i}(t)}$$

With these definitions in mind, net survival is interpreted as the probability that a patient remains alive in a hypothetical scenario where the disease of interest is the only potential cause of death. Estimating this from real-world data involves the assumption that the hazard $\lambda_{E_i}$ remains constant when other causes are eliminated. 

Some additional notations:

- The count number of events $dN_i(t)$ for individual $i$ at time $t$ 
- The total number of events $dN(t) = \sum_i dN_i(t)$
- The count $N_i(t) = \int_0^t dN_i(s)$ starts at 0 and skips to 1 when the individual $i$ dies
- The indicator $Y_i(t)$ for whether an individual is still at risk 
- The total number at risk $Y(t) = \sum_i Y_i(t)$ at time $t$

We can now define some key functions used in net survival analysis, such as the cumulative hazard function given by: 

$$\Lambda_{P_i}(t) = \int_0^t \lambda_{P_i}(u)du$$

As well as the population survival function given by:

$$S_{P_i}(t) = exp\{-\Lambda_{P_i}(t)\}$$

The relative survival estimator is therefore defined as such:

$$\hat{S}_R(t) = \frac{\hat{S}_O(t)}{\hat{S}_P(t)}$$

With $\hat{S}_O(t)$ being the estimated overall survival defined by the overall cumulative hazard $\Lambda_{O_i}(t) = \int_0^t \frac{dN(s)}{Y(s)}$, and $\hat{S}_P(t) = \frac{1}{n} \sum_i S_{P_i}(t)$.

## Pohar Perme

The Pohar Perme[PoharPerme2012](@cite) is a statistical method used in survival analysis to estimate net survival probabilities, particularly designed to handle situations where covariates may change over time. The net survival function is defined as:

$$\hat{S}_{E}(t) = exp(-\hat{\Lambda}_{E_i}(t))$$

The $\hat{\Lambda}_E(t)$ is the associated cumulative hazard given by : 

$$\hat{\Lambda}_E(t) = \int_0^t \frac{\sum_i \frac{dN_i(u)}{S_{P_i}(u)}}{\sum_i \frac{Y_i(u)}{S_{P_i}(u)}} - \int_0^t \frac{\sum_i \frac{Y_i(u)}{S_{P_i}(u)}d\Lambda_{P_i}(u)}{\sum_i \frac{Y_i(u)}{S_{P_i}(u)}}$$

With its respective variance estimator:

$$\hat{VAR}(\hat{\Lambda}_E(t)) = \int_0^t \frac{J(u)}{(\sum_i \frac{Y_i(u)}{S_{p_i}(u)})^2}\sum_{i=1}^n \frac{dN_i(u)}{S^2_{P_i}}$$

where $J(t) = I(Y(t) > 0)$.

## Ederer I

For the Ederer I [Ederer1961](@cite) estimator, the net survival is estimated using:

$$\hat{S}_{E}(t) = \int_0^t \frac{dN(u)}{Y(u)} - \int_0^t \frac{\sum_i S_{P_i}(u)\lambda_{P_i}(u)}{\sum_i S_{P_i}(u)}du$$

## Ederer II

As for the Ederer II estimator:

$$\hat{S}_{E}(t) = \int_0^t \frac{dN(u)}{Y(u)} - \int_0^t \frac{\sum_i Y_i(u)\lambda_{P_i}(u)}{\sum_i Y_i(u)}du$$

## Hakulinen

## Grafféo Log-Rank Test

The Grafféo Log-Rank Test [GraffeoTest](@cite) was constructed as a complement to the Pohar Perme estimator, aiming to compare the net survival functions provided by the latter. The test  is designed to compare these functions across multiple groups, including stratified covariables, and to ultimately determine, with the given p-value, which covariables are impactful to the study. 

For this test, we first define the number of deaths caused by the event studied for a time $s$ within the group $h$ noted $N_{E,h}(s)$ and the process of individuals at risk within the same group $h$ at time $s$ noted $Y_h(s)$. Both of these values are weighted with the populational estimated survival for the given patient, same as in Pohar Perme. 

The $(H_0)$ hypothesis tested 

```@bibliography
Pages = ["theory.md"]
Canonical = false
```