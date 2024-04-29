# Practice

## Data 

For standard relative survival analysis, it is required to have two different datasets: the cohort we wish to study and the population mortality table. In this practice, we will be re-using the dataset in [Getting Started](getting_started.md), `colrec`, as well as the same mortality table, `slopop`.

### Cohort

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

This dataframe shows as the number of male and female patients. There isn't too big of a difference between the two. We can say this study includes both gender relatively equally, thus, reducing bias. 

With these two observations, it is also worth noting that colorectal cancer is most common with men and people older than 50.

In total, we note that we have $5971$ patients. By taking a look at the `status` variable, we can determine the deaths and censorship:

```@example 2
sum(colrec.status)
```

Out of the $5971$ patients, $4979$ have died. This is a very high mortality rate, and again, unfortunately common in cancer studies.

```@example 2
(nrow(colrec) - sum(colrec.status)) / nrow(colrec)
```

In other terms, the censorship rate is of 16.6%, meaning the event, in this case death, was not observed for only 16.6% of the patients.

### Mortalility table

We will be using the mortality table `slopop` taken from the `RateTables.jl` package as the study is done on Slovenian patients. 

```@example 2
slopop
``` 

This command allows us to check what other covariates the mortality table has besides `age` and `year`, both expressed in days and not years. The ratetable is then three dimensional, with the covariate `sex` added. For example, the daily hazard rate for a woman turning $45$ on the January 1st $2006$ can be accessed through the following command:

```@example 2
daily_hazard(slopop, 45*365.241, 2006*365.241; sex=:female)
``` 

Making the survival probability easily calculated with:

```@example 2
exp(-(daily_hazard(slopop, 45*365.241, 2006*365.241; sex=:female))*365)
``` 

## Overall and expected survival

We noticed already in the previous graphs that there are many deaths occuring in the beginning of the study, mainly the first $5$ years.

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

With the EdererI method, and at time $1826$ days, we can say that the survival rate at this mark is around $0.456$ with only $0.017$% of patients dying from causes other than colorectal cancer.

```@example 2
e2 = fit(EdererII, @formula(Surv(time5,status5)~1), colrec, slopop)
```

The EdererII method, also known as the conditional method, shows that at the $5$ year mark, the survival probability is of $0.44$ with the same rate of patients from other causes. 

```@example 2
pp = fit(PoharPerme, @formula(Surv(time5,status5)~1), colrec, slopop)
```

We conclude, that in a world where cancer patients could only die due to cancer, only 41% of these patients would still be alive $5$ year after their diagnosis.

We will plot the Pohar Perme method only.

```@example 2
conf_int = confint(pp; level = 0.05)
lower_bounds = [lower[1] for lower in conf_int]
upper_bounds = [upper[2] for upper in conf_int] 

plot(pp.grid, pp.Sₑ, ribbon=(pp.Sₑ - lower_bounds, upper_bounds - pp.Sₑ), xlab = "Time (days)", ylab = "Net survival", label = false, title="Net survival distribution using Pohar Perme")
```

Looking at the rgaph, and its quick dip, it is evident that the first $5$ years are crucial and that the survival probability is highly affected in these years.

## Net survival with respect to covariates

We are now interested in comparing the different groups of patients defined by various covariates. We will look at individuals at age $65$ and above for now. 

```@example 2
pp_sex = fit(PoharPerme, @formula(Surv(time5,status5)~sex), colrec, slopop)
println(pp_sex[:,1826][1],pp_sex[:,1826][2])
```

```@example 2
colrec.age65 .= Bool(false)
for i in 1:nrow(colrec)
    if colrec.age[i] >= 65*365.241 
        colrec.age65[i] = true
    end
end

pp_age65 = fit(PoharPerme, @formula(Surv(time5,status5)~age65), colrec, slopop)
```

