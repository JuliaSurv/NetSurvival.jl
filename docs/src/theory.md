```@meta
CurrentModule = NetSurvival
```

# Theory

In many population-based studies, the specific cause of death is unidentified, unreliable or even unavailable. Relative survival analysis addresses this scenario, previously unexplored in general survival analysis. Different methods were created with the aim to construct a consistant and reliable estimator for this purpose.

## General Notations

Consider a study that includes censured data; we represent the different times with $E$,$P$ and $C$ with $E$ representing the time related to death caused by the event studied, $P$ the time related to all other causes, and $C$ the censorship. 

We will denote $O = E \land P$ as the overall time of death, regardless of the cause. $T= O \land C$ and $\Delta = \mathbf{1}\{T \leq C\}$ as well as $\mathbf{X}$, a vector of covariates, are considered observable in net survival. It is important to note that we have no information concerning the dependency structure of $(E,P)$ nor on a potential indicator $\mathbf{1}\{E \geq P\}$. This is one of the key differences between net survival and standard survival.

We will assume the overall hazard $\lambda_{O_i}$ for each individual $i$ is the sum between the hazard relative to the general population $\lambda_{P_i}$ and the one relative to the disease $\lambda_{E_i}$. $\lambda_{P_i}$ is typically derived from a reference life table alongside the covariates $\mathbf{X}$. Note that is also standard practice to assume that the distribution of $E$ and $C$ as well as the dependency structure of $(E, P, C)$ are independant from the covariates. 

Given our focus on the disease-specific hazard, we define the individual relative survival ratio as:

$$S_{E_i}(t) = exp \{-\int_0^t \lambda_{E_i}(u)du\} = \frac{exp \{-\int_0^t \lambda_{O_i}(u)du\}}{exp \{-\int_0^t \lambda_{P_i}(u)du\}}$$

With these definitions in mind, net survival is interpreted as the probability that a patient remains alive in a hypothetical scenario where the disease of interest is the only potential cause of death. Estimating this from real-world data involves the assumption that the hazard $\lambda_{E_i}$ remains constant when other causes are eliminated. 

Some additional notations:

- The count number of events $dN_i(t)$ for individual $i$ at time $t$ 
- The total number of events $dN(t) = \sum_i dN_i(t)$
- The count $N_i(t) = \int_0^t dN_i(s)$ starts at 0 and skips to 1 when the individual $i$ dies
- The indicator $Y_i(t)$ for whether an individual is still at risk 
- The total number at risk $Y(t) = \sum_i Y_i(t)$ at time $t$

As well as the population survival function given by:

$$S_{P_i}(t) = exp\{-\Lambda_{P_i}(t)\}$$

With

$$\Lambda_{P_i}(t) = \int_0^t \lambda_{P_i}(u)du$$

With these definitions and assumptions in mind, we will now present the four different methods implemented in this package commonly used to estimate $\partial\hat{\Lambda}_E(t)$.

## Ederer I

For the Ederer I [Ederer1961](@cite) estimator, the net survival is estimated using:

$$\partial\hat{\Lambda}_E(t) = \frac{dN(u)}{Y(u)} - \frac{\sum_i S_{P_i}(u)d\Lambda_{P_i}(u)}{\sum_i S_{P_i}(u)}$$

## Ederer II

As for the Ederer II [Ederer1959](@cite) estimator:

$$\partial\hat{\Lambda}_E(t) = \frac{dN(u)}{Y(u)} - \frac{\sum_i Y_i(u)d\Lambda_{P_i}(u)}{\sum_i Y_i(u)}$$

## Hakulinen

The Hakulinen method [Hakulinen1977](@cite) was later introduced and estimates the net survival using:

$$\partial\hat{\Lambda}_E(t) = \frac{dN(u)}{Y(u)} - \frac{\sum_i \frac{Y_i(u)}{ S_{P_i}(u)}d\Lambda_{P_i}(u)}{\sum_i \frac{Y_i(u)}{ S_{P_i}(u)}}$$

The variance is estimated in the same way for all three previous methods:

$$\partial\hat{\sigma}_E^2(t) =  J(u)\frac{\sum_i N_i(u)}{(\sum_i Y_i(u))^2}$$

where $J(t) = I(Y(t) > 0)$.

## Pohar Perme

The Pohar Perme[PoharPerme2012](@cite) is the newest addition to relative survival analysis between the four methods, particularly designed to handle situations where covariates may change over time. The estimation below is weighted with $S_{P_i}$, the populational survival rate.  

$$\partial\hat{\Lambda}_E(t) = \frac{\sum_i \frac{dN_i(u)}{S_{P_i}(u)} - \sum_i \frac{Y_i(u)}{S_{P_i}(u)}d\Lambda_{P_i}(u)}{\sum_i \frac{Y_i(u)}{S_{P_i}(u)}}$$

With its respective variance estimator:

$$\partial\hat{\sigma}_E^2(t) =  \frac{J(u)}{(\sum_i \frac{Y_i(u)}{S_{p_i}(u)})^2}\sum_{i=1}^n \frac{dN_i(u)}{S^2_{P_i}}$$

## Grafféo Log-Rank Test

The Grafféo Log-Rank Test [Graffeo2016](@cite) was constructed as a complement to the Pohar Perme estimator, aiming to compare the net survival functions provided by the latter. The test  is designed to compare these functions across multiple groups, including stratified covariables, and to ultimately determine, with the given p-value, which covariables are impactful to the study. The null $(H_0)$ hypothesis tests the following assumption:

$$\forall t \in [0,T], \; \; \Lambda_{E,1}(t) = \Lambda_{E,2}(t) = ... = \Lambda_{E,k}(t)$$

For this test, we first define the number of deaths caused by the event studied for a time $s$ within the group $h$ noted $N_{E,h}^w(s)$ and the process of individuals at risk within the same group $h$ at time $s$ noted $Y_h^w(s)$. Both of these values are weighted with the populational estimated survival for the given patient, same as in Pohar Perme. 

Thus, we define $\forall h \in [1;k]$:

$$Z_h^w(T) = \int_0^T \mathbf{1}(Y_.^w(s) > 0)dN_{E,h}^w(s) - \int_0^T \mathbf{1}(Y_.^w(s) > 0)\frac{Y_h^w(s)}{Y_.^w(s)}dN_{E,.}^w(s)$$

With $Y_.^w(s) = \sum_{h=i}^k Y_h^w(s)$ and $dN_{E,.}^w(s) = \sum_{h=i}^k dN_{E,h}^w(s)$.

The test statistic, defined under $H_0$, is then given by:

$$U^w(T) = (Z_1^w(T), ..., Z_{k-1}^w(T)) \hat{\Sigma}_0^w(T)^{-1} (Z_1^w(T), ..., Z_{k-1}^w(T))^t$$

Knowing that $\hat{\Sigma}_0^w(T) = (\hat{\sigma}^w_{hj}(T))_{hj}$ for $(h,j) \in [1;k-1]^2$

Finally, the asymptotic distribution of the test statistic under $H_0$:

$$U^w(t) \sim \chi^2(k-1)$$

We reject the $H_0$ hypothesis when the p-value obtained is under $0.05$, and therefore, admitting the notable difference between the groups. 

## References

```@bibliography
Pages = ["theory.md"]
Canonical = false
```