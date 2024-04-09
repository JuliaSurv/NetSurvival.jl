```@meta
CurrentModule = NetSurvival
```

# Getting Started

## Pohar Perme

The Pohar Perme[PoharPerme2012](@cite) is a statistical method used in survival analysis to estimate net survival probabilities, particularly designed to handle situations where covariates may change over time. The net survival function is defined as:

$$S_{E}(t) = exp(-\int_0^t\lambda_{E}(u)du)$$

The $\lambda_E$ is the associated hazard given by : 

$$\lambda_E (t) = \frac{\sum_{i=1}^{N}S_{E_i}(t)\lambda_{E_i}(t)}{\sum_{i=1}^{N}S_{E_i}(t)}$$

This weighted average is thus based on the likelihood that an individual remains alive in a hypothetical scenario where the disease is the sole cause of death. 

## Ederer II

## Grafféo Log-Rank Test

The Grafféo Log-Rank Test [GraffeoTest](@cite) was constructed as a complement to the Pohar Perme estimator, aiming to compare the net survival functions provided by the latter. The test  is designed to compare these functions across multiple groups, including stratified covariables, and to ultimately determine, with the given p-value, which covariables are impactful to the study. 

For this test, we first define the number of deaths caused by the event studied for a time $s$ within the group $h$ noted $N_{E,h}(s)$ and the process of individuals at risk within the same group $h$ at time $s$ noted $Y_h(s)$. Both of these values are weighted with the populational estimated survival for the given patient, same as in Pohar Perme. 

The $(H_0)$ hypothesis tested 

```@bibliography
Pages = ["getting_started.md"]
Canonical = false
```