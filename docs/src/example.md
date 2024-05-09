# Case study

## Introduction and datasets

We will illustrate with an example using the dataset `colrec`, which comprises $5971$ patients diagnosed with colon or rectal cancer  between $1994$ and $2000$. This dataset is sourced from the Slovenia cancer registry. Given the high probability that the patients are Slovenian, we will be using the Slovenian mortality table `slopop` as reference for the populational rates. Subsequently, we can apply various non-parametric estimators for net survival analysis.

!!! note "N.B." 
    Mortality tables may vary in structure, with options such as the addition or removal of specific covariates. To confirm that the mortality table is in the correct format, please refer to [`RateTables.jl`'s documentation](https://JuliaSurv.github.io/RateTables.jl/), or directly extract it from there.

### Cohort details

The patients in the study are diagnosed between January 1st 1994 and December 31st 2000. Before we move on to the survival probabilities, it is important to be aware of how your data is distributed and of what it comprises. 

```@example 2
using NetSurvival, RateTables, DataFrames

first(colrec,10)
```

Let's explore how our data is distributed first, strating with the `age` variable.

```@example 2
println(minimum(colrec.age./365.241), " ; ",maximum(colrec.age./365.241))
```

The study can be considered diverse in terms of age seeing as the patients are between 12 and 96 years old, approximately. 

```@example 2
using Plots 

plot(
    histogram(colrec.age./365.241, label="Age"),
    histogram(colrec.time./365.241, label="Follow-up time")
)
```

The graph above show us that although it has a wide range of patients within all age groups, it is mostly centered around older adults and elderly, with the majority of the patients being between 60 and 80 years old. 

Looking at the second graph that details the distibution of the follow-up times, we notice that the values quickly drop. Unfortunately, this is a common theme in cancer studies. 

Let's take a look at the `sex` variable now: 

```@example 2
combine(groupby(colrec, :sex), nrow)
```

This dataframe shows the number of male and female patients. There isn't too big of a difference between the two. We can say this study includes both gender relatively equally, thus, reducing bias. 

With these two observations, it is also worth noting that colorectal cancer is most common with men and people older than 50.

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

By examining `slopop`, we notice it contains information regarding `age` and `year`, as expected for mortality tables. Additionally, it incorporates the covariate sex, which has two possible entries (`:male` or `:female`). The ratetable is then three dimensional, with the covariate `sex` added. For example, the daily hazard rate for a woman turning $45$ on the January 1st $2006$ can be accessed through the following command:

```@example 2
daily_hazard(slopop, 45*365.241, 2006*365.241; sex=:female)
``` 

Making the survival probability easily calculated with:

```@example 2
exp(-(daily_hazard(slopop, 45*365.241, 2006*365.241; sex=:female))*365)
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

In this part, we are interested in the first $5$ years of the study. We will thus limit the follow-up time to $5$ years, meaning we will censor all individuals with a follow-up time that is higher than this. Then, we will apply the different net survival methods.

```@example 2 
colrec.time5 .= 0.0
colrec.status5 .= Bool(true)
for i in 1:nrow(colrec)
    colrec.time5[i] = min(colrec.time[i], round(5*365.241))
    if colrec.time[i] > 5*365.241
        colrec.status5[i] = false
    end
end
```

Now that we have defined our own time and status variables according to the observations made, we can apply the different non parametric methods for relative survival.

```@example 2
e1 = fit(EdererI, @formula(Surv(time5,status5)~1), colrec, slopop)
```

With the EdererI method, after $1826$ days have passed, we can say that the survival rate at this mark is around $0.456$, in the hypothetical world where patients can only die of cancer.

```@example 2
crude_e1 = CrudeMortality(e1)
println(crude_e1.Λₒ[1826], " , ", crude_e1.Λₑ[1826], " , ", crude_e1.Λₚ[1826])
```

Out of the 0.63 patients that have died, according to the EdererI method, 0.51 died because of colorectal cancer and 0.12 died of other causes.


```@example 2
e2 = fit(EdererII, @formula(Surv(time5,status5)~1), colrec, slopop)
```

Similarily, the EdererII method, also known as the conditional method, shows that at the $5$ year mark, the survival probability is of $0.44$ in this hypothetical world.

```@example 2
crude_e2 = CrudeMortality(e2)
println(crude_e2.Λₒ[1826], " , ", crude_e2.Λₑ[1826], " , ", crude_e2.Λₚ[1826])
```

Here, out of the 0.63 patients that have dued, 0.53 are due to colorectal cancer and 0.1 due to other causes.

```@example 2
pp = fit(PoharPerme, @formula(Surv(time5,status5)~1), colrec, slopop)
```

We conclude for the Poher-Perme method, that in a world where cancer patients could only die due to cancer, only 41% of these patients would still be alive $5$ year after their diagnosis.

```@example 2
crude_pp = CrudeMortality(pp)
println(crude_pp.Λₒ[1826], " , ", crude_pp.Λₑ[1826], " , ", crude_pp.Λₚ[1826])
```

Finally, for this estimator, we have that of the 0.64 patients that have died, 0.53 is due to colorectal cancer while 0.11 is due to other causes.

We will plot the Pohar Perme method only.

```@example 2
conf_int = confint(pp; level = 0.05)
lower_bounds = [lower[1] for lower in conf_int]
upper_bounds = [upper[2] for upper in conf_int] 

p1 = plot(pp.grid, pp.Sₑ, ribbon=(pp.Sₑ - lower_bounds, upper_bounds - pp.Sₑ), xlab = "Time (days)", ylab = "Net survival", label = false)

p2 = plot(pp.grid, crude_pp.Λₑ, label = "Excess Mortality Rate")
p2 = plot!(pp.grid, crude_pp.Λₚ, label = "Population Mortality Rate")

plot(p1,p2)
```

Looking at the graph, and the rapid dip it takes, it is evident that the first $5$ years are crucial and that the survival probability is highly affected in these years. Additionnally, the crude mortality graph allows us to see how much of this curve is due to the colorectacl cancer studied versus other undefined causes. It is clear that the large majority is due to the cancer.

## Net survival with respect to covariates

We are now interested in comparing the different groups of patients defined by various covariates. 

```@example 2
pp_sex = fit(PoharPerme, @formula(Surv(time5,status5)~sex), colrec, slopop)
```

When comparing at time $1826$, we notice that the survival probability is slightly inferior for men than for women ($0.433 < 0.449$). It is also more probable for the women to die from other causes than the men seeing as $0.0255 > 0.025$. Still, the differences are minimal. Let's confirm this with the Grafféo log-rank test:

```@example 2
test_sex = fit(GraffeoTest, @formula(Surv(time5,status5)~sex), colrec, slopop)
```

The p-value is indeed above $0.05$. We cannot reject the null hypothesis $H_0$ and thus we dismiss the differences between the two sexes.

As for the age, we will define two different groups: individuals aged 65 and above and those who are not.

```@example 2
colrec.age65 .= colrec.age .>= 65 * 365.241
pp_age65 = fit(PoharPerme, @formula(Surv(time5,status5)~age65), colrec, slopop)
```

Here, the difference between the two is much more important. In the first group, the individuals are aged under 65 and at $5$ years time, they have a $50.1$% chance of survival. On the other hand, the individuals aged 65 and up have a $40.1$% chance of survival. 

It is also worth noting that their chances of dying from other causes is higher than the younger group, given their age. 

When applying the Grafféo test, we get the results below:

```@example 2
test_age65 = fit(GraffeoTest, @formula(Surv(time5,status5)~age65), colrec, slopop)
```

The p-value is well under $0.05$, meaning we reject the $H_0$ hypothesis and must admit there are differences between the individuals aged 65 and above and the others.

When plotting both we get:


```@example 2
conf_int_men = confint(pp_sex[1]; level = 0.05)
lower_bounds_men = [lower[1] for lower in conf_int_men]
upper_bounds_men = [upper[2] for upper in conf_int_men] 

conf_int_women = confint(pp_sex[2]; level = 0.05)
lower_bounds_women = [lower[1] for lower in conf_int_women]
upper_bounds_women = [upper[2] for upper in conf_int_women]

conf_int_under65 = confint(pp_age65[1]; level = 0.05)
lower_bounds_under65 = [lower[1] for lower in conf_int_under65]
upper_bounds_under65 = [upper[2] for upper in conf_int_under65] 

conf_int_65 = confint(pp_age65[2]; level = 0.05)
lower_bounds_65 = [lower[1] for lower in conf_int_65]
upper_bounds_65 = [upper[2] for upper in conf_int_65] 

plot1 = plot(pp_sex[1].grid, pp_sex[1].Sₑ, ribbon=(pp_sex[1].Sₑ - lower_bounds_men, upper_bounds_men - pp_sex[1].Sₑ), xlab = "Time (days)", ylab = "Net survival", label = "men")

plot1 = plot!(pp_sex[2].grid, pp_sex[2].Sₑ, ribbon=(pp_sex[2].Sₑ - lower_bounds_women, upper_bounds_women - pp_sex[2].Sₑ), xlab = "Time (days)", ylab = "Net survival", label = "women")

plot2 = plot(pp_age65[1].grid, pp_age65[1].Sₑ, ribbon=(pp_age65[1].Sₑ - lower_bounds_under65, upper_bounds_under65 - pp_age65[1].Sₑ), xlab = "Time (days)", ylab = "Net survival", label = "Under 65")

plot2 = plot!(pp_age65[2].grid, pp_age65[2].Sₑ, ribbon=(pp_age65[2].Sₑ - lower_bounds_65, upper_bounds_65 - pp_age65[2].Sₑ), xlab = "Time (days)", ylab = "Net survival", label = "65 and up")

plot(plot1, plot2, layout = (1, 2))
```

Visually, it is almost immediately understood that there are no worthy differences between the two sexes whereas the `age65` variable seems to play a big role.
