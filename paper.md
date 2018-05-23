---
title: Environmental Justice Analysis of Chlorpyrifos Use in California's Central Valley
author: Daniel J. Hicks
bibliography: spatial_project.bib
---

# Introduction #

This paper conducts an distributional environmental justice (dEJ) analysis of chlorpyrifos use in California's Central Valley.  Chlorpyrifos is one of the most widely-used pesticides in the world.  It is known to be a moderately powerful neurotoxicant, and was banned for residential use in the US in 2001.  During the Obama administration, US EPA had developed plans to remove all tolerances for chlorpyrifos residue on food, which effectively would have banned all US use of chlorpyrifos.  However, one of Scott Pruitt's first actions as Administrator of US EPA was to halt this process.  [Need cites, actual lit review]

*["distributional"]

An dEJ analysis examines the way the distribution of environmental risks intersect with race, ethnicity, class, gender, and other systems of structural oppression.  *[UCC study and since]  While spatial methods are frequently used in dEJ analysis, they are often not statistically sophisticated; @ChakrabortyRevisitingToblerFirst2011 notes that EJ regression analyses often use ordinary least squares (OLS), and fail to account for spatial autocorrelation.  OLS may also be inappropriate for some left-bounded data, such as pounds of a chemical used or disease rates.  This lack of statistical sophistication makes spatial EJ analyses vulnerable to technical criticism by "merchants of doubt" [@OreskesMerchantsDoubtHow2011], which can limit their effectiveness in pushing for policy change or legal remediation.  

*[real lit review] A quick literature search identified one previous spatial analysis of chlorpyrifos exposure in California's Central Valley [@LuoSpatiallydistributedpesticide2010].  This analysis focused on modeling fate and transport of the pesticide using a physical-chemical simulation, and did not examine the population being exposed.  It was therefore not a dEJ analysis.  In contrast, my analysis will focus on demographic covariates for chlorpyrifos use in the Central Valley.  
*[Lovasi, G. S., Quinn, J. W., Rauh, V. A., Perera, F. P., Andrews, H. F., Garfinkel, R., ... & Rundle, A. (2011). Chlorpyrifos exposure and urban residential environment characteristics as determinants of early childhood neurodevelopment. American journal of public health, 101(1), 63-70.; 
CHAMACOS study]


# Data #

## Pesticide Use ##

This paper combines two datasets.  First, California's Department of Pesticide Regulation (DPR) releases annual public datasets for pesticide use across the state, known Pesticide Use Reports (PUR) [^DPR].  These datasets provide data on pesticide use at the township and section level, including pounds of an active chemical used on particular days and times.  The most recent data release covers 2015.  These data exclusively report agriculture use; this limitation is acceptable for studying chlorpyrifos in the agriculture-heavy Central Valley.  

[^DPR]: <http://www.cdpr.ca.gov/docs/pur/purmain.htm>

I retrieved the full datasets for 2011-2015.  Combining data over several years can help smooth annual variation in use due to variation in insect pest prevalence or other factors.  These datasets were combined and filtered on chemical name for chlorpyrifos.  Data for all counties across the state were used, not just the Central Valley, to avoid edge effects.  For example, while Sacramento County was not included in the study area, there is some chlorpyrifos use in the southeastern part of Sacramento County, very close to areas of San Joaquin County that are included in the study area.  

DPR conducts automated quality checks, and flags errors and changed values in the data.  I closely examined entries with "potential duplicate" and "outlier" error codes.  For the first few years of the dataset (2011 through 2013), it appears that DPR flagged but did not remove potential duplicates.  In later years, it appears that DPR removes potential duplicates in the released data.  I removed potential duplicates from the earlier years.  

DPR automatically identifies potential outliers as uses with *[def'n].  Inspecting the entries with "potential outlier" flags, it appears that the main datafile includes estimated or modified values, not the original potential outliers.  I used modified values without further adjustment.  

Chlorpyrifos uses are associated with townships and sections in the data, and annual totals of active ingredient were calculated at the centroid of each 1 mile-square ($1.6 km \times 1.6 km$) section.  Because these centroids do not match the actual use locations (i.e., farm fields), the centroid totals might be unreliable for the smallest CTD, 1km (see discussion and table \ref{tab.ctd} below).  However, this error should be negligible for the other CTDs.  


## Demographics ##

The second primary dataset comprises ACS five-year estimates, from 2011-2015, for all census tracts and places in 17 Central Valley Counties.  For the purposes of this study, I included Solano County (which has substantial chlorpyrifos use in the Fairfield-Vacaville area) but excluded Sacramento County (which is primarily urban, and has only a relatively small agricultural area in the southwest).  

*[control variables:  population density; percent employed in agriculture; housing values]

For each region, I retrieved estimates and margin of error (MOE) values for three categories of demographic variables:  *race and ethnicity* (Hispanic, non-Hispanic White, non-Hispanic Black, Indigenous, and Asian residents), *foreign-born noncitizens*, *children under 5* (who may be especially sensitive to chlorpyrifos exposure; *[cite]), and *poverty* (individuals with an income-poverty ratio below 1).  (See also table \ref{tab.dv}.)  I calculated population densities and proportions (e.g., the fraction of all residents who are Hispanic) for the total population and each of these demographic variables for every tract and place, using Census-recommended methods to calculate MOEs for these derived variables *[cite].  A few *[how many?] tracts and places with an estimated total population of 0 were excluded from all analyses.  

Because much of the study area is rural, tract size and population density varied over four orders of magnitude, from *[low to high].  Places cover *[85%] of the population, including *[non-Hispanic Whites and Hispanic residents], and have *[more consistent densities?].  However, places are geographically separated from each other; in constructing contiguity-based spatial weights, *[62%] of places had no neighbors.  

The Modifiable Areal Unit Problem (MAUP) has been used to criticize dEJ projects [@SteelEnvironmentalJusticeValues2012; TaylorToxiccommunitiesenvironmental2014, 41ff].  "Egocentric neighborhood" methods have been used to address the MAUP in other fields *[cite].  However, this method assumes that populations are spread evenly over each discrete region (e.g., each census tract).  This assumption is inappropriate for this project, which includes many  rural regions.  In addition, the large geographic scale of this project would require millions — if not tens of millions — of egocentric neighborhoods, and so literally trillions of distance calculations between neighborhoods and chlorpyrifos uses.  The census tract-level distance calculations already pressed the limits of the computing power available for this project.  More fine-grained regions (e.g., census block groups or blocks) would have multiplied uncertainty in the ACS estimates, and also likely would have exceeded the available computing power.  

As a compromise, I used block population counts from the 2010 Census to calculate weighted centroids for each tract *[and place?].  These centroids more accurately represent the "average location" of the population in each tract, without requiring more computing power for the distance calculations.  *[Lievanos?? uses census tract centroids] ^[NLCD]

^[NLCD]: An alternative approach would be to use the National Land Cover Database (NLCD), which was most recently released in 2011.  <https://www.mrlc.gov/nlcd2011.php>  Specifically, tract populations could be interpolated across the four categories of developed land identified in the NLCD.  Unfortunately, this public data source did not come to my attention until too late in the process of preparing this manuscript.  

Chlorpyrifos use centroids, census tracts, and places included in the study area are shown in figure \ref{fig.chlor_use}.  

![Data used in this study.  Red points are chlorpyrifos use totals, added across 2011-2015 and shown on a $log_10$ pounds scale.  Blue regions are census tracts included in the study area; yellow regions are included places.  All California counties are shown for context; counties included in the Central Valley for the purposes of this study are named. \label{fig.chlor_use}](03_chlor_use.png)


## Linking Pesticide Use to Tracts and Places ##

Chlorpyrifos use is linked to tracts and places using the concept of *characteristic travel distance* (CTD).  CTD is defined as the distance "at which 63% of the original mass of volatilized [chemical] is degraded or deposited" [@MackayFateenvironmentlongrange2014].  In the simplest sense, CTD models assume exponential decay:  

\begin{align}
    q(d) &= q_0 \beta^d \label{eqn.decay}
\end{align}

where $q(d)$ is the quantity of a chemical at a distance $d$ from its source; $q_0$ is the quantity volatilized at the source; and $\beta$ is a constant related to the CTD:  

\begin{align} 
    (1-0.63) &= \beta^{CTD}\\
    \beta &= 0.37^{1/CTD}
\end{align}

I will refer to $\beta^d$ as the *decay coefficient*.  Given a CTD, the decay coefficient can be used in equation \ref{eqn.decay} to produce a crude estimate of potential exposure at a distance $d$ from a chemical source, e.g., a use of chlorpyrifos associated with a section centroid.  Then, for each location (tract or place) $i$, the total potential exposure at the location centroid can be estimated as 
\begin{align}
    q_i &= \sum_u q_u \beta^{d_{iu}} \ref{eqn.tot_pot_exp}\\
        &= \sum_u q_u 0.37^{d_{iu}/CTD}
\end{align}
where $q_u$ is the total use at section $u$ and $d_{iu}$ is the Euclidean distance between the centroid of $i$ and the centroid of $u$.  To slightly account for the fact that the residents of a location are not located at its centroid, the decay coefficient is set to 1 whenever the centroid associated with use $u$ is within $i$, regardless of $d_{iu}$.  

It is important to stress that equations \ref{eqn.decay} and \ref{eqn.tot_pot_exp} are, at best, crude estimate of potential exposure.  They do not take into account prevailing or occurrent winds, or other chemical transport mechanisms.  They represent the chemical moving uniformly away from the source along a one-dimensional path, not dispersing over a two- or three-dimensional space surrounding the source at various speeds in various directions at various times.  It represents diverse and variable processes of application, fixation, and chemical transformation as simple exponential decay.  By comparison, @LuoSpatiallydistributedpesticide2010 use a physio-chemical model and PUR data to produce much more sophisticated estimates of chlorpyrifos loading.  

Since physio-chemical modeling is beyond our capacities here, we standardize the total use values $q_u$ (across all sections with non-zero use values).  That is, we work with relative rather than absolute potential exposures, and assume that on average over five years fate and transport processes are symmetrical around each use.  

@MackayFateenvironmentlongrange2014 estimate the CTD for chlorpyrifos to be 62 km.  For robustness, this paper considers five CTD values, ranging from 1 km to 90 km.  See table \ref{tab.ctd}.  

| CTD (km) | $\beta$ |
|:---------|:--------|
| 1        | 0.370   |
| 10       | 0.905   |
| 30       | 0.967   |
| 60       | 0.984   |
| 90       | 0.989   |

Table: Characteristic Travel Distance (CTD) values used in this study, and corresponding scaling values $\beta$. \label{tab.ctd}


# Methods #

The primary analysis of this study is a regression of relative potential chlorpyrifos exposure against Census demographic data.  We construct separate models for each of the five CTD values listed in table \ref{tab.ctd}, as well as separate models for tracts and places.  These two methodological choices give $5 \times 2 = 10$ models.  The software language R was used to clean and analyze all data *[cite], and the Bayesian regression models were written in Stan *[cite].  Complete cleaning and analysis code is available at *[ref].  

In the terminology of Bogen and Woodward [-@BogenSavingPhenomena1988], this project aims to characterize a phenomenon rather than confirm or disconfirm a theory.  Bogen and Woodward argue that phenomena are intermediate between data and theories.  Data can be thought of as the more-or-less direct outputs of instruments — in the case of this project, the PUR forms filled out by farmers and survey instruments used by the Census.  Phenomena are more-or-less stable patterns or trends in the world, and are presumed to reflect "signal" rather than "noise" in the data.  But phenomena are not explanations of the patterns or trends.  That is the role of theories, which offer models of the causal processes that produce phenomena.  The distinction between phenomena and theories is close to the distinction between correlation and causation:  phenomena describe correlations, while theories aim at causation.  Thus, this study aims to describe phenomenal relationships between chlorpyrifos use and population demographics, but not to explain them.  

Note that the theory-phenomenon distinction does not mean that phenomena are theory-independent.  I have already introduced theoretical assumptions about, for instance, how chlorpyrifos decays over space.  More generally, separating "signal" and "noise" requires background assumptions about what kinds of processes might be at work in the data-generating process.  But this project merely relies on those assumptions; it does not intend to confirm or deny them.  

Arguably, the theory-phenomenon distinction aligns with the difference between procedural and distributive aspects of EJ.  The unjust distribution of environmental hazards is the outcome of social processes, which are often unjust in themselves.  

Non-spatial exploratory data analysis indicated that, for almost all (>80%) tracts and places, almost all residents (>80%) were either Hispanic or non-Hispanic White.  Thus, Hispanic and non-Hispanic White proportions are strongly negatively correlated ($r = -.9$), and we drop non-Hispanic White rate from the independent variable list.  Similarly, there are strong correlations between Hispanic and non-citizen proportions ($r = .8$), and so we also drop non-citizen proportion from the independent variable list.  (Alternatively, factor analysis methods might have been used to construct composite variables combining Hispanic and non-citizen proportions; compare *[Lievandos??].  This approach was not taken for simplicity.)  There are moderate correlations ($r = .4-.6$) between Hispanic, young children, and poverty proportions; so, while we use all three of these variables, we expect to see greater uncertainty in the Bayesian posterior estimates for these three variables.  

Intersectionality theory *[cite] suggests that interlocking systems of oppression, such as race and class, can produce "synergistic" effects, with greater disparities than either system would produce "in isolation."  We therefore considered including a Hispanic poverty proportion in the independent variable list.  However, taking the proportion with respect to the total population (Hispanic poverty / total population) yielded a very strong correlation with the Hispanic proportion (Hispanic / total population; $r = .9$); and taking the proportion with respect to the Hispanic population (Hispanic poverty / Hispanic) yielded a very strong correlation with general poverty (poverty / total population; $r = .8$).  We therefore did not include this intersectional variable here.  

We also include two control variables.  Because PUR data come only from agricultural uses, and chlorpyrifos is not registered for domestic use, we expect chlorpyrifos use to be negatively correlated with population density and positively correlated with agricultural employment.  Density was not substantially correlated with any variables of interest (except for a low-moderate correlation negative correlation with White; $r = -.4$).  Agricultural employment has a moderate-strong correlation with Hispanic ($r = .67$), so *[we consider models both with and without this control].  

All independent variables are given in table \re{tab.dv}.  

| Category       | Independent Variable |
|                | (proportion of total population) |
|:---------------|:-------------------|
| Race-Ethnicity | Hispanic           |
|                | Black              |
|                | Indigenous         |
|                | Asian              |
| Children       | Children under 5   |
| Class          | Income-Poverty Ratio < 1.0 |
| Controls       | Population Density |
|                | Employed in Agriculture | 

Table: Independent variables (IVs) used in this study. All race-ethnicity groups other than Hispanic are non-Hispanic.  All (first-order) IVs are proportion of total population in the tract or place. \label{tab.dv}

Preliminary plots indicated that potential exposure values were multimodal and highly non-Gaussian.  However, plotting separate distributions for each county suggested that this was due to very different county-level baselines; for counties with many tracts or places, the distribution of values appeared to be sufficiently normal for Gaussian regression.  We include a county-level random intercept to account for these different baselines.  

Spatial exploratory data analysis of both independent and dependent variables suggested substantial degrees of spatial autocorrelation on both sides of the model.  We therefore use a spatial Durbin model.  

Finally, ACS estimates can have large margins of error, especially for subpopulations of rural tracts.  We therefore incorporate an error-in-variables approach for the dependent variables.  However, I was unable to find a R package that produced maximum likelihood estimates for combined spatial Durbin-error-in-variables models.  I therefore take a Bayesian approach, and use the flexibility of Stan to specify such a model.  

All together, the regression model is specified as follows:  

*[prior on Xtrue]

\begin{align*}
    X^{obs} &\sim N(X^{true}, S^2 I_n)\\
    Y &\sim N(\hat Y, \sigma^2 I_n)\\
    \hat Y &= X^{true}\beta + \gamma_{county} + W_{ind}X_{true}\rho + W_{dep}Y\theta\\
    \gamma_{county} &\sim N(mu_{county}, \sigma^2_{county})\\
    \\
    X^{true} &\sim
    \beta_j, \theta, \rho &\sim N(0, 5^2)
    \mu_{county} &\sim N(0, 5^2)\\
    \sigma, \sigma_{county} &\sim HCauchy(0, 5)\\
\end{align*}
where 

- $Y$ is the length-$n$ dependent variable column vector, $\hat Y$ is the fitted dependent variable, $\sigma$ is the residual standard deviation and $I_n$ is the $n\times n$ identity matrix; 
- $X^{obs} = \left< x_{i,j} \right>$ is the matrix of observed values of the $j$th predictor (demographic variable) at location $i$,  $X^{true}$ is the corresponding matrix of true predictor values, and $S$ is the matrix of standard deviations for these observed values, calculated from ACS margins of error; 
- $\beta$ is a length-$j$ column vector of regression coefficients; 
- $\gamma_{county}$ are the county-level random intercepts, which are drawn from a Gaussian distribution with second-level parameters $\mu_{county}$ and $\sigma_{county}$; 
- $W_{ind}$ and $W_{dep}$ are the spatial weights matrices for the independent and dependent variables, respectively; and 
- $\rho$ and $\theta$ are the spatial lag coefficients for the independent and dependent variables, respectively.  

The final three rows in the specification give weakly informative priors on the top-level parameters.  The regression coefficients $\beta, \theta, \rho$ have (independent) Gaussian priors, centered at 0 with standard deviation 10.  These coefficients are in units of standard deviations of relative chlorpyrifos potential exposure per change of rate from 0% to 100%.  $\mu_{county}$ also has a Gaussian prior centered at 0 with standard deviation 10.  These random intercepts are in units of standard deviations of relative chlorpyrifos potential exposure.  Finally, the variance parameters $\sigma, \sigma_{county}$ are given (independent) half-Cauchy priors (where "half" means that it is constrained to positive values) with location 0 and scale 5 [@GelmanPriordistributionsvariance2006].  

Besides fitting the model, our Stan code also calculates estimates for $R^2$ and Moran's $I$ on the residuals, to produce Bayesian posterior distributions for the fit and residual spatial autocorrelation.  


# References #
