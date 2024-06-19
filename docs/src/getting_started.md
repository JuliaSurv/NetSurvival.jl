```@meta
CurrentModule = NetSurvival
```

# Getting Started

In many population-based studies, the specific cause of death is unidentified, unreliable or even unavailable. Relative survival analysis addresses this scenario, previously unexplored in general survival analysis. Different methods were created with the aim to construct a consistant and reliable estimator for this purpose.

## Net survival settings

Note first that, for any positive random variable $X$, we will use extensively the following functions (that each fully characterize the distribution of $X$): 

| Function symbol and definition | Function Name |
| :--- | :---- |
| ``S_X(t) = \mathbb P(X > t)`` | Survival function |
| ``\Lambda_X(t) = -\ln S_X(t)`` | Cumulative Hazard function |
| ``\lambda_X(t) = \partial \Lambda_X(t)`` | Instantaneous hazard function |


Consider a study that consists of censored survival times from a specific cause. Such a study consists of several random objects: 

| Random variable    | Name         | Is it observed ?     |
| :----------------- | :---------------------------- | :-: |
| ``E``              | "Excess" lifetime              | ✘ |
| ``P``              | "Population" lifetime          | ✘ |
| ``O = E \wedge P`` | "Overall" lifetime             | ✘ |
| ``C``              | "Censoring" time               | ✘ |
| ``\mathbf D``      | Vector of covariates           | ✔ |
| ``T = O \wedge C`` | Event time                     | ✔ |
| ``\Delta = \mathbf{1}\{T \leq C\}`` | Event status | ✔ |

It is important to note that we do not observe a a potential indicator $\mathbf{1}\{E \geq P\}$. This is one of the key differences between net survival and standard survival. The standard Net survival analysis solves this problem by assuming that the underlying times $E$ and $P$ are independent from each other.

One central hypothesis in net survival (on top of non-informative censoring) is that $P$ and $E$ are independent. This independence can be written in terms of hazard rates $\lambda_O(t) = \lambda_P(t) + \lambda_E(t)$, or in terms of survival functions $S_O(t) = S_P(t)S_E(t)$.

The population hazard for each individual $\lambda_{P_i}$ is usally drawn from a reference life table, and may depend on covariates $\mathbf D_i$ such as age and date, sex, country, race, etc... See the [RateTables.jl](https://github.com/JuliaSurv/RateTables.jl) package for more details on the potential covariates. On the other hand, the excess mortality is assumed to be i.i.d. between individuals and not to depend on covariates at all. Thus, we mostly omit these covariates from our notations.

## Available Estimators

The estimation of net survival is usually discussed in terms of the estimation of the cumulative excess hazard $\Lambda_E(t)$ and/or the instantaneous hazard $\lambda_E = \partial\Lambda_E$. To describe the estimators, we use the following counting processes notations, similar to standard survival analysis(see e.g. [FlemingHarington2013](@cite) or [ABGK1993](@cite)). 
* The uncensored event indicatrix $\partial N_i(t)$ for individual $i$ at time $t$ 
* The total number of uncensored events process $\partial N(t) = \sum_i \partial N_i(t)$ at time $t$
* The at-risk indicatrix $Y_i(t)$, for whether an individual is still at risk 
* The total number at risk process $Y(t) = \sum_i Y_i(t)$ at time $t$

With these definitions and assumptions in mind, we will now present the four different methods implemented in this package, commonly used in literature, to estimate the excess hazard function $\partial\Lambda_E(t)$ and its variance. Recall that we estimate the variance as $\sigma_E^2(t) = \int_{0}^t \partial\sigma_E^2(s)$. 

| Name        | Reference               | Proposed (partial) Excess Hazard $\partial\hat{\Lambda}_E(s)$ | Proposed (partial) Variance $\partial\hat{\sigma}_E^2(s)$ |
| :---------- | :---------------------- | :------------------------------------ | :---------------------------------------- | 
| Ederer I    | [Ederer1961](@cite)     | $\frac{\sum_i N_i(s)}{\sum_i Y_i(s)} - \frac{\sum_i S_{P_i}(s)\partial\Lambda_{P_i}(s)}{\sum_i S_{P_i}(s)}$ | $\frac{\sum_i N_i(s)}{\left(\sum_i Y_i(s)\right)^2}$ |
| Ederer II   | [Ederer1959](@cite)     | $\frac{\sum_i N_i(s)}{\sum_i Y_i(s)} - \frac{\sum_i Y_i(s)\partial\Lambda_{P_i}(s)}{\sum_i Y_i(s)}$ | $\frac{\sum_i N_i(s)}{\left(\sum_i Y_i(s)\right)^2}$ |
| Hakulinen   | [Hakulinen1977](@cite)  | $\frac{\sum_i N_i(s)}{\sum_i Y_i(s)} - \frac{\sum_i \frac{Y_i(s)}{ S_{P_i}(s)}\partial\Lambda_{P_i}(s)}{\sum_i \frac{Y_i(s)}{ S_{P_i}(s)}}$ | $\frac{\sum_i N_i(s)}{\left(\sum_i Y_i(s)\right)^2}$ |
| Pohar Perme | [PoharPerme2012](@cite) | $\frac{\sum_i \frac{\partial N_i(s)}{S_{P_i}(s)} - \sum_i \frac{Y_i(s)}{S_{P_i}(s)}\partial\Lambda_{P_i}(s)}{\sum_i \frac{Y_i(s)}{S_{P_i}(s)}}$ | $\frac{\sum_{i=1}^n \frac{\partial N_i(s)}{S^2_{P_i}}}{\left(\sum_i \frac{Y_i(s)}{S_{p_i}(s)}\right)^2}$ |

where, in the variances, it is understood that when no more individuals are at risk $0/0$ gives $0$. 

The Pohar Perme estimator [PoharPerme2012](@cite) is the newest addition to relative survival analysis between the four methods, particularly designed to handle situations where covariates may change over time. It is trusted from the field (see e.g. [PermePavlik2018](@cite) and [CharvatBelot2021](@cite)) that only this estimator should really be used, the other ones being included mostly for historical reasons and comparisons. 


```@docs
PoharPerme
EdererI
EdererII
Hakulinen
```

## Crude Mortality

The *crude mortality rate* is the global mortality rate (from all causes of death) for a population. This measure does not use the cause of death information which, as we previously mentioned, can be unreliable and incomplete. It is given as:

$$M_E(t) = \int_0^t S_O(u-) \lambda_E(u)du$$

There exists a few estimators of this quantity, the most known one being the Cronin-Feuer estimator [cronin2000cumulative](@cite), given by:

$$\hat{M}_E(t) = \int_0^t \hat{S}_O(u-) \hat{\lambda}_E(u)du$$

where, traditionally, Kaplan-Meier is used to estimate the overall survival function $\hat{S}_O$ and Ederer II is used for the excess hazard rate $\hat{\lambda}_E$ and the population hazard rate $\hat{\lambda}_P$.

Our implementation is a bit more permissive, as any net survival estimators can be used for $\hat{\lambda}_E$. Of course, the default is still the original Ederer II which provides the original Cronin-Feuer estimator. 
 
```@docs
CrudeMortality
```

## Grafféo Log-Rank Test

The Grafféo Log-Rank Test [Graffeo2016](@cite) was constructed as a complement to the Pohar Perme estimator, aiming to compare the net survival functions provided by the latter. The test  is designed to compare these functions across multiple groups, including stratified covariables, and to ultimately determine, with the given p-value, which covariables are impactful to the study.

The null $(H_0)$ hypothesis tests the following assumption:

$$\forall t \in [0,T], \; \; \Lambda_{E,g_1}(t) = \Lambda_{E,g_2}(t) = ... = \Lambda_{E,g_k}(t),$$

where $G = \{g_1,...,g_k\}$ is a partition of $1,...,n$ consisting of disjoint groups of individuals that we wish to compare to each other. 
For all group $g \in G$, let's denote the numerator and denominator of the Pohar Perme (partial) excess hazard estimators, restricted to individuals in the group, by: 

* ``\partial N_{E,g}(s) = \sum_{i \in g} \frac{\partial N_i(s)}{S_{P_i}(s)} - \frac{Y_i(s)}{S_{P_i}(s)}\partial\Lambda_{P_i}(s)``
* ``Y_{E,g}(s) = \sum_{i \in g} \frac{Y_i(s)}{S_{P_i}(s)}``
* ``R_{g}(s) = \frac{Y_{E,g}(s)}{\sum_{g\in G} Y_{E,g}(s)}``

Then, define the vector $\mathbf Z = \left(Z_{g_r}: \; r \in 1,...,k-1 \right)$ with entries: 

$$Z_g(T) = N_{E,g}(s) - \int_{0}^T Y_{E,g}(s) \partial\hat{\Lambda}_E(s)$$

The test statistic is then given by:

$$U(T) = \mathbf Z(T)'\hat{\Sigma}_Z^{-1} \mathbf Z(T)$$

where the entries of the $\hat{\Sigma}_Z$ matrix are given by: 

$$\sigma_{g,h}(T) = \int_0^T \sum_{\ell \in G} \left(\delta_{g,\ell} - R_g(t) \right)\left(\delta_{h,\ell} - R_h(t)\right) \left(\sum_{i\in\ell} \frac{\partial N_i(s)}{S^2_{P_i}}\right)$$

Under $H_0$, the statistic $U(T)$ is asymptotically $\chi^2(k-1)$-distributed. We thus reject the $H_0$ hypothesis when the p-value obtained is under $0.05$, admitting the notable difference between the groups. 

```@docs
GraffeoTest
```

## Nessie

The Nessie function estimates the sample size by yearly intervals as well as averages an estimated lifespan left for a given group.  

This function is highly dependant on the `Life` function taken from the `RateTables.jl` package which you can find documented [here](https://juliasurv.github.io/RateTables.jl/dev/).

The sample size is thus taken by the following formula:

$$ESS(t) = \sum_i^N S_{P_i}(t)$$

While the estimated lifepsan is directly taken from the `expectation` function. 

```@docs
nessie
```

## References

```@bibliography
Pages = ["getting_started.md"]
Canonical = false
```
