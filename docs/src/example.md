# Case study

## Introduction and datasets

We will illustrate with an example using the dataset `colrec`, which comprises $5971$ patients diagnosed with colon or rectal cancer between $1994$ and $2000$. This dataset is sourced from the Slovenia cancer registry. Given the high probability that the patients are Slovenian, we will be using the Slovenian mortality table `slopop` as reference for the population mortality rates. Subsequently, we can apply various non-parametric estimators for net survival analysis.

!!! note "N.B." 
    Mortality tables may vary in structure, with options such as the addition or removal of specific covariates. To confirm that the mortality table is in the correct format, please refer to [`RateTables.jl`'s documentation](https://JuliaSurv.github.io/RateTables.jl/), or directly extract it from there.

### Cohort details

The patients in the study are diagnosed between January 1st $1994$ and December 31st $2000$. Before we move on to the survival probabilities, it is important to be aware of how your data is distributed and of what it comprises. 

```@example 2
using NetSurvival, RateTables, DataFrames
first(colrec,10)
```

Let's explore how our data is distributed first, starting with the `age` variable.

```@example 2
println(minimum(colrec.age./365.241), " ; ",maximum(colrec.age./365.241))
```

The study can be considered diverse in terms of age seeing as the patients are between 12 and 96 years old, approximately. 

```@example 2
using Plots 
plot(histogram(colrec.age./365.241, label="Age"),
    histogram(colrec.time./365.241, label="Follow-up time"))
```

The graph above show us that although the dataset has a wide range of patients within all age groups, it is mostly centered around older adults and elderly, with the majority of the patients being between $60$ and $80$ years old. Looking at the second graph that details the distribution of the follow-up times. We notice there that the values quickly drop. Unfortunately, this is a common theme in cancer studies. 

Let's take a look at the `sex` variable now, by looking at the number of male and female patients:

```@example 2
combine(groupby(colrec, :sex), nrow)
```

There isn't too big of a difference between the two. We can say this study includes both gender relatively equally, thus, reducing bias. With these two observations, it is also worth noting that colorectal cancer is most common with men and people older than $50$.

In total, we note that we have $5971$ patients. By taking a look at the `status` variable, we can determine the deaths and censorship:

```@example 2
sum(colrec.status)
```

Out of the $5971$ patients, $4979$ have died. This is a very high mortality rate, and again, unfortunately common in cancer studies.

```@example 2
(nrow(colrec) - sum(colrec.status)) / nrow(colrec)
```

In other terms, the censorship rate is of 16.6%, meaning the event, in this case death, was not observed for only 16.6% of the patients. This is a low censorship rate and thus the quality of the signal will be pretty good. 

### Mortality table

We will be using the mortality table `slopop` taken from the `RateTables.jl` package as the study is done on Slovenian patients. 

```@example 2
slopop
``` 

The show method of the `RateTable` class shows the additional covariate `sex` that the rate table has on top of the (mandatory) `age` and `year` variables. The `sex` variable has two madalities, `:male` and `:female`. The ratetable is then three dimensional. For example, the daily hazard rate for a woman turning $45$ on the January 1st $2006$ can be accessed through the following command:

```@example 2
λ  = daily_hazard(slopop, 45*365.241, 2006*365.241; sex=:female)
``` 

Making the daily survival probability easily calculated with:

```@example 2
exp(-λ)
``` 

## Overall and expected survival

For this part, we will be using the `Survival.jl` package to apply the Kaplan Meier estimator for the overall survival.

```@example 2
using Survival 
km = fit(KaplanMeier, colrec.time./365.241, colrec.status)
plot(km.events.time, km.survival, label=false, title = "Kaplan-Meier Estimator for the Overall Survival")
```

The graph above indicates a significant dip in survival probability within the first $5$ years, and even more at $10$ years. This period in the study is then crucial for the analysis. 

## Estimated net survival

We will restrict ourselves to the first $5$ years of the study. For that, let us re-censor the dataset as follows: 

```@example 2 
for i in 1:nrow(colrec)
    if colrec.time[i] > 1826 # five years
        colrec.status[i] = false
        colrec.time[i] = 1826
    end
end
```

We can now apply the different non-parametric methods to compute the relative survival.

```@example 2
e1 = fit(EdererI, @formula(Surv(time,status)~1), colrec, slopop)
```

With the EdererI method, after $1826$ days have passed, we can say that the survival rate at this mark is around $0.456$, in the hypothetical world where patients can only die of cancer.

```@example 2
crude_e1 = CrudeMortality(e1)
println(crude_e1.Mₒ[1826], " , ", crude_e1.Mₑ[1826], " , ", crude_e1.Mₚ[1826])
```

Out of the 0.63 patients that have died, according to the EdererI method, 0.51 died because of colorectal cancer and 0.12 died of other causes.

```@example 2
e2 = fit(EdererII, @formula(Surv(time,status)~1), colrec, slopop)
```

Similarly, the EdererII method, also known as the conditional method, shows that at the $5$ year mark, the survival probability is of $0.44$ in this hypothetical world.

```@example 2
crude_e2 = CrudeMortality(e2)
println(crude_e2.Mₒ[1826], " , ", crude_e2.Mₑ[1826], " , ", crude_e2.Mₚ[1826])
```

Here, out of the 0.63 patients that have died, 0.53 are due to colorectal cancer and 0.1 due to other causes.

```@example 2
pp = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
```

We conclude for the Pohar Perme method, that in a world where cancer patients could only die due to cancer, only 41% of these patients would still be alive $5$ year after their diagnosis. The Pohar Perme estimator is the best estimator of the excess hazard under the standard hypotheses. 

```@example 2
crude_pp = CrudeMortality(pp)
println(crude_pp.Mₒ[1826], " , ", crude_pp.Mₑ[1826], " , ", crude_pp.Mₚ[1826])
```

Finally, for this estimator, we have that of the 0.64 patients that have died, 0.53 is due to colorectal cancer while 0.11 is due to other causes.

We will plot the Pohar Perme method only.

```@example 2
function mkribbon(pp)
    S = pp.Sₑ
    ci = confint(pp; level = 0.05)
    l,u = getindex.(ci, 1), getindex.(ci,2)
    rb = (S - l, u - S)
    return rb
end

p1 = plot(pp.grid, pp.Sₑ, ribbon=mkribbon(pp), xlab = "Time (days)", ylab = "Net survival", label = false)

p2 = plot(pp.grid, crude_pp.Mₑ, label = "Excess Mortality Rate")
p2 = plot!(pp.grid, crude_pp.Mₚ, label = "Population Mortality Rate")

plot(p1,p2)
```

Looking at the graph, and the rapid dip it takes, it is evident that the first $5$ years are crucial and that the survival probability is highly affected in these years. Additionally, the crude mortality graph allows us to see how much of this curve is due to the colorectal cancer studied versus other undefined causes. It is clear that the large majority is due to the cancer.

## Net survival with respect to covariates

We are now interested in comparing the different groups of patients defined by various covariates. 

```@example 2
pp_sex = fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, slopop)
pp_males = pp_sex[pp_sex.sex .== :male,:estimator][1]
pp_females = pp_sex[pp_sex.sex .== :female,:estimator][1]
```

When comparing at time $1826$, we notice that the survival probability is slightly inferior for men than for women ($0.433 < 0.449$). It is also more probable for the women to die from other causes than the men seeing as $0.0255 > 0.025$. Still, the differences are minimal. Let's confirm this with the Grafféo log-rank test:

```@example 2
test_sex = fit(GraffeoTest, @formula(Surv(time,status)~sex), colrec, slopop)
```

The p-value is indeed above $0.05$. We cannot reject the null hypothesis $H_0$ and thus we dismiss the differences between the two sexes.

As for the age, we will define two different groups: individuals aged 65 and above and those who are not.

```@example 2
colrec.age65 .= ifelse.(colrec.age .>= 65*365.241, :old, :young)
pp_age65 = fit(PoharPerme, @formula(Surv(time,status)~age65), colrec, slopop)
pp_young = pp_age65[pp_age65.age65 .== :young, :estimator][1]
pp_old = pp_age65[pp_age65.age65 .== :old, :estimator][1]
```

Here, the difference between the two is much more important. In the first group, the individuals are aged under 65 and at $5$ years time, they have a $50.1$% chance of survival. On the other hand, the individuals aged 65 and up have a $40.1$% chance of survival. 

It is also worth noting that their chances of dying from other causes is higher than the younger group, given their age. 

When applying the Grafféo test, we get the results below:

```@example 2
test_age65 = fit(GraffeoTest, @formula(Surv(time,status)~age65), colrec, slopop)
```

The p-value is well under $0.05$, meaning we reject the $H_0$ hypothesis and must admit there are differences between the individuals aged 65 and above and the others.

When plotting both we get:


```@example 2
plot1 = plot(pp_males.grid, pp_males.Sₑ, ribbon=mkribbon(pp_males), xlab = "Time (days)", ylab = "Net survival", label = "men")
plot1 = plot!(plot1, pp_females.grid, pp_females.Sₑ, ribbon=mkribbon(pp_females), label = "women")

plot2 = plot(pp_young.grid, pp_young.Sₑ, ribbon=mkribbon(pp_young), xlab = "Time (days)", ylab = "Net survival", label = "Under 65")
plot2 = plot!(pp_old.grid, pp_old.Sₑ, ribbon=mkribbon(pp_old), label = "65 and up")

plot(plot1, plot2, layout = (1, 2))
```

Visually, it is almost immediately understood that there are no worthy differences between the two sexes whereas the `age65` variable seems to play a big role.

The same kind of graph can be made on the stage: 
```@example 2
pp3 = fit(PoharPerme, @formula(Surv(time,status)~stage), colrec, slopop)
plot3 = plot(xlab = "Time (days)", ylab = "Net survival",title="Net survival per cancer stages")
for i in 1:nrow(pp3)
    e = pp3[i,:estimator]
    plot!(plot3, e.grid, e.Sₑ, ribbon=mkribbon(e), label=pp3[i,:stage])
end
plot(plot3)
```



## Estimated sample size and life expectancy 

Given that the age group plays a significant role in the study, we will now estimate the sample size by yearly intervals in order to better compare the age groups.

```@example 2
elt, ess = nessie(@formula(Surv(time,status)~age65), colrec, slopop)
elt
```

The expected life time for the younger patients is significatively higher than for older patients (24.78 years > 10.29 years).

```@example 2
hcat(ess[:,3]...)
```

Finally, the table above represents yearly expected sample sizes for both age groups under 65 and above, with the second column representing the latter. We can see that the sample size decreases for the older patients in a much more dramatic way than for the younger ages.

Unsurprisingly, we can thus conclude that age plays an important role in the study.