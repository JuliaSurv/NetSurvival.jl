var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"F. Ederer. The relative survival rate: a statistical methodology. Natl. Cancer Inst. Monogr. 6, 101–121 (1961).\n\n\n\nF. Ederer and H. Heise. The effect of eliminating deaths from cancer in general survival rates, methodological notes 11. End Result Evaluation Section, National Cancer Institute (1959).\n\n\n\nT. Hakulinen. On long-term relative survival rates. Journal of Chronic Diseases 30, 431–443 (1977).\n\n\n\nM. P. Perme, J. Stare and J. Estève. On Estimation in Relative Survival. Biometrics 68, 113–120 (2011).\n\n\n\nN. Grafféo, F. Castell, A. Belot and R. Giorgi. A Log-Rank-Type Test to Compare Net Survival Distributions. Biometrics 72, 760–769 (2016).\n\n\n\nT. R. Fleming and D. P. Harrington. Counting Processes and Survival Analysis. Vol. 625 (John Wiley & Sons, 2013).\n\n\n\nP. K. Andersen, O. Borgan, R. D. Gill and N. Keiding. Statistical Models Based on Counting Processes. Springer Series in Statistics (Springer US, New York, NY, 1993).\n\n\n\nM. P. Perme and K. Pavlic. Nonparametric relative survival analysis with the R package relsurv. Journal of Statistical Software 87, 1–27 (2018).\n\n\n\nH. Charvat and A. Belot. Mexhaz: An R package for fitting flexible hazard-based regression models for overall and excess mortality with a random effect. Journal of Statistical Software 98, 1–36 (2021).\n\n\n\n","category":"page"},{"location":"references/","page":"References","title":"References","text":"T. Hakulinen and K. H. Abeywickrama. A computer program package for relative survival analysis. Computer programs in biomedicine 19, 197–207 (1985).\n\n\n\nT. Hakulinen and L. Tenkanen. Regression analysis of relative survival rates. Journal of the Royal Statistical Society Series C: Applied Statistics 36, 309–317 (1987).\n\n\n\n","category":"page"},{"location":"benches/","page":"Benchmarking results","title":"Benchmarking results","text":"CurrentModule = NetSurvival","category":"page"},{"location":"benches/#Benchmarking-results","page":"Benchmarking results","title":"Benchmarking results","text":"","category":"section"},{"location":"benches/","page":"Benchmarking results","title":"Benchmarking results","text":"The following benchmarks are run on github actions continuous integration platform, which is a very slow computing engine. Local experiments suggests performances that are twice as fast on correct hardware – note that we do not use multithreading at all, but underlying BLAS calls might. ","category":"page"},{"location":"benches/","page":"Benchmarking results","title":"Benchmarking results","text":"using RCall\nusing NetSurvival, RateTables, BenchmarkTools\n\nR_bench = @benchmark R\"\"\"\nrelsurv::rs.surv(\n    survival::Surv(time, stat) ~1, \n    rmap=list(age = age, sex = sex, year = diag), \n    data = relsurv::colrec, \n    ratetable = relsurv::slopop, \n    method = \"pohar-perme\", \n    add.times=1:8149)\n\"\"\"\n\njl_bench = @benchmark fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)\n\nratio = time(minimum(R_bench)) / time(minimum(jl_bench))","category":"page"},{"location":"benches/#Benchmarking-across-time","page":"Benchmarking results","title":"Benchmarking across time","text":"","category":"section"},{"location":"benches/","page":"Benchmarking results","title":"Benchmarking results","text":"The folloiwng charts provide a glimpse of NetSurvival.jl's performance along time: ","category":"page"},{"location":"benches/","page":"Benchmarking results","title":"Benchmarking results","text":"<iframe src=\"../../benchmarks/\" style=\"height:500px;width:100%;\"></iframe>","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"CurrentModule = NetSurvival","category":"page"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"In many population-based studies, the specific cause of death is unidentified, unreliable or even unavailable. Relative survival analysis addresses this scenario, previously unexplored in general survival analysis. Different methods were created with the aim to construct a consistant and reliable estimator for this purpose.","category":"page"},{"location":"theory/#General-Notations","page":"Theory","title":"General Notations","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"Note first that, for any positive random variable X, we will use extensively the following functions (that each fully characterize the distribution of X): ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Function symbol and definition Function Name\nS_X(t) = mathbb P(X  t) Survival function\nLambda_X(t) = -ln S_X(t) Cumulative Hazard function\nlambda_X(t) = partial lambda_X(t) Instantaneous hazard function","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Consider a study that consists of censored survival times from a specific cause. Such a study consists of several random objects: ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Random variable Name Is it observed ?\nE \"Excess\" lifetime ✘\nP \"Population\" lifetime ✘\nO = E wedge P \"Overall\" lifetime ✘\nC \"Censoring\" time ✘\nmathbf D Vector of covariates ✔\nT = O wedge C Event time ✔\nDelta = mathbf1T leq C Event status ✔","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"It is important to note that we do not observe a a potential indicator mathbf1E geq P. This is one of the key differences between net survival and standard survival. The standard Net survival analysis solves this problem by assuming that the underlying times E and P are independent from each other.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"One central hypothesis in net survival (on top of non-informative censoring) is that P and E are independent. This independence can be written in terms of hazard rates lambda_O(t) = lambda_P(t) + lambda_E(t), or in terms of survival functions S_O(t) = S_P(t)S_E(t).","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The population hazard for each individual lambda_P_i is usally drawn from a reference life table, and may depend on covariates mathbf D_i such as age and date, sex, country, race, etc... See the RateTables.jl package for more details on the potential covariates. On the other hand, the excess mortality is assumed to be i.i.d. between individuals and not to depend on covariates at all. Thus, we mostly omit these covariates from our notations.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The estimation of net survival is usually discussed in terms of the estimation of the cumulative excess hazard Lambda_E(t) and/or the instantaneous hazard lambda_E = partialLambda_E. To describe the estimators, we use the following counting processes notations, similar to standard survival analysis(see e.g. [6] or [7]). ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The uncensored event indicatrix partial N_i(t) for individual i at time t \nThe total number of uncensored events process partial N(t) = sum_i partial N_i(t) at time t\nThe at-risk indicatrix Y_i(t), for whether an individual is still at risk \nThe total number at risk process Y(t) = sum_i Y_i(t) at time t","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"With these definitions and assumptions in mind, we will now present the four different methods implemented in this package, commonly used in literature, to estimate the excess hazard function partialLambda_E(t) and its variance. Recall that we estimate the variance as sigma_E^2(t) = int_0^t partialsigma_E^2(s). ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Name Reference Proposed (partial) Excess Hazard partialhatLambda_E(s) Proposed (partial) Variance partialhatsigma_E^2(s)\nEderer I [1] fracsum_i N_i(s)sum_i Y_i(s) - fracsum_i S_P_i(s)partialLambda_P_i(s)sum_i S_P_i(s) fracsum_i N_i(s)left(sum_i Y_i(s)right)^2\nEderer II [2] fracsum_i N_i(s)sum_i Y_i(s) - fracsum_i Y_i(s)partialLambda_P_i(s)sum_i Y_i(s) fracsum_i N_i(s)left(sum_i Y_i(s)right)^2\nHakulinen [3] fracsum_i N_i(s)sum_i Y_i(s) - fracsum_i fracY_i(s) S_P_i(s)partialLambda_P_i(s)sum_i fracY_i(s) S_P_i(s) fracsum_i N_i(s)left(sum_i Y_i(s)right)^2\nPohar Perme [4] fracsum_i fracpartial N_i(s)S_P_i(s) - sum_i fracY_i(s)S_P_i(s)partialLambda_P_i(s)sum_i fracY_i(s)S_P_i(s) fracsum_i=1^n fracpartial N_i(s)S^2_P_ileft(sum_i fracY_i(s)S_p_i(s)right)^2","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where, in the variances, it is understood that when no more individuals are at risk 00 gives 0. ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The Pohar Perme estimator [4] is the newest addition to relative survival analysis between the four methods, particularly designed to handle situations where covariates may change over time. It is trusted from the field (see e.g. [8] and [9]) that only this estimator should really be used, the other ones being included mostly for historical reasons and comparisons. ","category":"page"},{"location":"theory/#Grafféo-Log-Rank-Test","page":"Theory","title":"Grafféo Log-Rank Test","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"The Grafféo Log-Rank Test [5] was constructed as a complement to the Pohar Perme estimator, aiming to compare the net survival functions provided by the latter. The test  is designed to compare these functions across multiple groups, including stratified covariables, and to ultimately determine, with the given p-value, which covariables are impactful to the study.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The null (H_0) hypothesis tests the following assumption:","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"forall t in 0T   Lambda_Eg_1(t) = Lambda_Eg_2(t) =  = Lambda_Eg_k(t)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where G = g_1g_k is a partition of 1n consisting of disjoint groups of individuals that we wish to compare to each other.  For all group g in G, let's denote the numerator and denominator of the Pohar Perme (partial) excess hazard estimators, restricted to individuals in the group, by: ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"partial N_Eg(s) = sum_i in g fracpartial N_i(s)S_P_i(s) - fracY_i(s)S_P_i(s)partialLambda_P_i(s)\nY_Eg(s) = sum_i in g fracY_i(s)S_P_i(s)\nR_g(s) = fracY_Eg(s)sum_gin G Y_Eg(s)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Then, define the vector mathbf Z = left(Z_g_r  r in 1k-1 right) with entries: ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Z_g(T) = N_Eg(s) - int_0^T Y_Eg(s) partialhatLambda_E(s)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The test statistic is then given by:","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"U(T) = mathbf Z(T)hatSigma_Z^-1 mathbf Z(T)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where the entries of the hatSigma_Z matrix are given by: ","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"sigma_gh(T) = int_0^T sum_ell in G left(delta_gell - R_g(t) right)left(delta_hell - R_h(t)right) left(sum_iinell fracpartial N_i(s)S^2_P_iright)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Under H_0, the statistic U(T) is asymptotically chi^2(k-1)-distributed. We thus reject the H_0 hypothesis when the p-value obtained is under 005, admitting the notable difference between the groups. ","category":"page"},{"location":"theory/#References","page":"Theory","title":"References","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"F. Ederer. The relative survival rate: a statistical methodology. Natl. Cancer Inst. Monogr. 6, 101–121 (1961).\n\n\n\nF. Ederer and H. Heise. The effect of eliminating deaths from cancer in general survival rates, methodological notes 11. End Result Evaluation Section, National Cancer Institute (1959).\n\n\n\nT. Hakulinen. On long-term relative survival rates. Journal of Chronic Diseases 30, 431–443 (1977).\n\n\n\nM. P. Perme, J. Stare and J. Estève. On Estimation in Relative Survival. Biometrics 68, 113–120 (2011).\n\n\n\nN. Grafféo, F. Castell, A. Belot and R. Giorgi. A Log-Rank-Type Test to Compare Net Survival Distributions. Biometrics 72, 760–769 (2016).\n\n\n\nT. R. Fleming and D. P. Harrington. Counting Processes and Survival Analysis. Vol. 625 (John Wiley & Sons, 2013).\n\n\n\nP. K. Andersen, O. Borgan, R. D. Gill and N. Keiding. Statistical Models Based on Counting Processes. Springer Series in Statistics (Springer US, New York, NY, 1993).\n\n\n\nM. P. Perme and K. Pavlic. Nonparametric relative survival analysis with the R package relsurv. Journal of Statistical Software 87, 1–27 (2018).\n\n\n\nH. Charvat and A. Belot. Mexhaz: An R package for fitting flexible hazard-based regression models for overall and excess mortality with a random effect. Journal of Statistical Software 98, 1–36 (2021).\n\n\n\n","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"CurrentModule = NetSurvival","category":"page"},{"location":"getting_started/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"getting_started/#Fitting-the-non-parametric-estimators","page":"Getting Started","title":"Fitting the non parametric estimators","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"PoharPerme\nEdererI\nEdererII\nHakulinen","category":"page"},{"location":"getting_started/#NetSurvival.PoharPerme","page":"Getting Started","title":"NetSurvival.PoharPerme","text":"PoharPerme\n\nThis method estimates net survival probabilities by applying the following estimation:\n\npartialhatLambda_E(t) = fracsum_i fracdN_i(u)S_P_i(u) - sum_i fracY_i(u)S_P_i(u)dLambda_P_i(u)sum_i fracY_i(u)S_P_i(u)\n\nTo fit the Pohar Perme to your data based on a certain rate table, apply the example below to your code : \n\nfit(PoharPerme, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)\n\nReferences: \n\n[4] Perme, Maja Pohar and Stare, Janez and Estève, Jacques (2012). On Estimation in Relative Survival.\n\n\n\n\n\n","category":"type"},{"location":"getting_started/#NetSurvival.EdererI","page":"Getting Started","title":"NetSurvival.EdererI","text":"EdererI\n\nTo call this function: \n\nfit(EdererI, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)\n\n\n\n\n\n","category":"type"},{"location":"getting_started/#NetSurvival.EdererII","page":"Getting Started","title":"NetSurvival.EdererII","text":"EdererII\n\nTo call this function: \n\nfit(EdererII, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)\n\n\n\n\n\n","category":"type"},{"location":"getting_started/#NetSurvival.Hakulinen","page":"Getting Started","title":"NetSurvival.Hakulinen","text":"Hakulinen\n\nTo call this function: \n\nfit(Hakulinen, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)\n\n\n\n\n\n","category":"type"},{"location":"getting_started/#Example","page":"Getting Started","title":"Example","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"We will illustrate with an example using the dataset colrec, which comprises 5971 patients diagnosed with colon or rectal cancer  between 1994 and 2000. This dataset is sourced from the Slovenia cancer registry. Given the high probability that the patients are Slovenian, we will be using the Slovenian mortality table slopop as reference for the populational rates. Subsequently, we can apply various non-parametric estimators for net survival analysis.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"note: N.B.\nMortality tables may vary in structure, with options such as the addition or removal of specific covariates. To confirm that the mortality table is in the correct format, please refer to the documentation of RateTables.jl, or directly extract it from there.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"By examining slopop, we notice it contains information regarding age and year, as expected for mortality tables. Additionally, it incorporates the covariate sex, which has two possible entries (:male or :female).","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Pohar Perme","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"using NetSurvival, RateTables\npp1 = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)","category":"page"},{"location":"getting_started/#Applying-the-Grafféo-log-rank-test","page":"Getting Started","title":"Applying the Grafféo log-rank test","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"GraffeoTest","category":"page"},{"location":"getting_started/#NetSurvival.GraffeoTest","page":"Getting Started","title":"NetSurvival.GraffeoTest","text":"GraffeoTest\n\nThe Grafféo test is a log-rank type test and is typically used in net survival analysis to determine the impact of certain covariates in the study.\n\nThe null (H_0) hypothesis tests the following assumption:\n\nforall t in 0T Lambda_Eg_1(t) = Lambda_Eg_2(t) =  = Lambda_Eg_k(t)\n\nwhere G = g_1g_k is a partition of 1n consisting of disjoint groups of individuals that we wish to compare to each other.  For all group g in G, let's denote the numerator and denominator of the Pohar Perme (partial) excess hazard estimators, restricted to individuals in the group, by: \n\npartial N_Eg(s) = sum_i in g fracpartial N_i(s)S_P_i(s) - fracY_i(s)S_P_i(s)partialLambda_P_i(s)\nY_Eg(s) = sum_i in g fracY_i(s)S_P_i(s)\nR_g(s) = fracY_Eg(s)sum_gin G Y_Eg(s)\n\nThen, define the vector mathbf Z = left(Z_g_r r in 1k-1 right) with entries: \n\nZ_g(T) = N_Eg(s) - int_0^T Y_Eg(s) partialhatLambda_E(s)\n\nThe test statistic is then given by:\n\nU(T) = mathbf Z(T)hatSigma_Z^-1 mathbf Z(T)\n\nwhere the entries of the hatSigma_Z matrix are given by: \n\nsigma_gh(T) = int_0^T sum_ell in G left(delta_gell - R_g(t) right)left(delta_hell - R_h(t)right) left(sum_iinell fracpartial N_i(s)S^2_P_iright)\n\nUnder H_0, the statistic U(T) is asymptotically chi^2(k-1)-distributed.\n\nTo apply the test to your data based on a certain rate table, apply the example below to your code : \n\nfit(GraffeoTest, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)\n\nIf you wish to stratify a covariate:\n\nfit(GraffeoTest, @formula(Surv(time,status)~covariable1 + Strata(covariable2)), data, ratetable)\n\nReferences: \n\n[5] Grafféo, Nathalie and Castell, Fabienne and Belot, Aurélien and Giorgi, Roch (2016). A Log-Rank-Type Test to Compare Net Survival Distributions.  \n\n\n\n\n\n","category":"type"},{"location":"getting_started/#Example-2","page":"Getting Started","title":"Example","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"When applying the test to the same data as before, we get:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"test1 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"The p-value is well under 005, meaning that the different groups identified by the stage variable have different survival probabilities. Thus, it should be taken into consideration in the study.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"test2 = fit(GraffeoTest, @formula(Surv(time,status)~sex), colrec, slopop)","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"For the sex variable, we notice that the p-value is above 005 indicating that there isn't a difference between male and female patients.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NetSurvival","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The NetSurvival.jl package provides the necessary tools to perform estimations and analysis in the Net Survival field. This specialized branch of Survival Analysis focuses on estimating the probability of survival from a specific event of interest, for example a given cancer, without considering other causes of death. This is especially relevant in the (unfortunately quite common) case where the cause of death indicatrix is either unavailable or untrustworthy. Consequently, the so-called missing indicatrix issue forbids the use of standard competitive risks survival analysis methods on these datasets.  For that, a few standard estimators were established in the last 50 years, backed by a wide literature.","category":"page"},{"location":"","page":"Home","title":"Home","text":"By integrating observed data from the target population with historical population mortality data (usually sourced from national census datasets), Net Survival allows the extraction of the specific mortality hazard associated with the particular disease, even under the missing indicatrix issue. The concept of relative survival analysis dates back several decades to the seminal article by Ederer, Axtell, and Cutler in 1961 [1] and the one by Ederer and Heise in 1959 [2].","category":"page"},{"location":"","page":"Home","title":"Home","text":"For years, the Hakulinen estimator (1977) [3] and the Ederer I and II estimators were widely regarded as the gold standard for non-parametric survival curve estimation. However, the introduction of the Pohar-Perme, Stare, and Estève estimator in 2012 [4] resolved several issues inherent in previous estimators, providing a reliable and consistent non-parametric estimator for net survival analysis.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Standard tools nowadays are composed of R packages, with underlying C and C++ routines, that are hard to read, maintain, and use. This package is an attempt to bring standard relative survival analysis modeling routines to Julia, while providing an interface that is close to the relsurv standard, albeit significantly faster and easier to maintain in the future. Our hope is that the junction with classical modeling API in Julia will allow later extensions of the existing modeling methods, with a simple interface for the practitioners.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some key features in NetSurvival.jl are:","category":"page"},{"location":"","page":"Home","title":"Home","text":"A panel of different non-parametric net survival estimators (Ederer I [1], Ederer II [2], Hakulinen [3], Pohar Perme [4]) with an interface compliant with Julia's standards. \nGrafféo's log-rank test [5] to compare net survival curves accross groups, including stratified testing.\nA compact, readable and efficient codebase (up to 1000x less LOC than relsurv for the same functionalities), ensuring long-term maintenability.\nSignificant performance improvements (up to 50x) compared to the R package relsurv.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not yet available on Julia's general registry, and thus can be installed through the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"https://github.com/JuliaSurv/NetSurvival.jl.git\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the rest of this documentation to have a glimpse of the functionalities!","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"F. Ederer. The relative survival rate: a statistical methodology. Natl. Cancer Inst. Monogr. 6, 101–121 (1961).\n\n\n\nF. Ederer and H. Heise. The effect of eliminating deaths from cancer in general survival rates, methodological notes 11. End Result Evaluation Section, National Cancer Institute (1959).\n\n\n\nT. Hakulinen. On long-term relative survival rates. Journal of Chronic Diseases 30, 431–443 (1977).\n\n\n\nM. P. Perme, J. Stare and J. Estève. On Estimation in Relative Survival. Biometrics 68, 113–120 (2011).\n\n\n\nN. Grafféo, F. Castell, A. Belot and R. Giorgi. A Log-Rank-Type Test to Compare Net Survival Distributions. Biometrics 72, 760–769 (2016).\n\n\n\n","category":"page"}]
}