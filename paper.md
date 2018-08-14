---
title: Environmental Justice Analysis of Chlorpyrifos Use in California's Central Valley
author: Daniel J. Hicks
bibliography: spatial_project.bib
---

# Introduction #

Chlorpyrifos is one of the most widely-used pesticides in the world. In California in 2016, it was the 29th most heavily used pesticide active ingredient, with over 900,000 pounds of active ingredient applied over 640,000 acres [@PesticideUseReporting2017].  It is generally recognized to be a moderate neurotoxin [@WorldHealthOrganizationWHORecommendedClassification2010 67; @USEPAChlorpyrifosRevisedHuman2016; citations in @TrasandeWhenenoughdata2017], and was banned for residential use in the US in 2001.  Based on the toxicological and epidemiological evidence, in 2007 the environmental organizations Pesticide Action Network North America (PANNA) and Natural Resources Defense Council (NRDC) filed a petition with US EPA, calling on the agency to revoke all tolerances for chlorpyrifos.  In 2017, US EPA rejected this petition [@USEPAChlorpyrifosOrderDenying2017], and as of this writing it remains in widespread use in California and elsewhere. 

This paper conducts an distributional environmental justice (dEJ) analysis of chlorpyrifos use in California's Central Valley.  A dEJ analysis examines the way the distribution of environmental risks intersect with race, ethnicity, class, gender, and other systems of structural oppression.  Since the landmark report "Toxic Wastes and Race in the United States"  [@CommissionforRacialJusticeToxicWastesRace1987], a significant dEJ literature has emerged [@PulidoEnvironmentalismEconomicJustice1996; @Shrader-FrechetteEnvironmentaljusticecreating2002; @BrownToxicExposuresContested2007; @OttingerTechnoscienceenvironmentaljustice2011; @TaylorToxiccommunitiesenvironmental2014].  

This paper uses methods of spatial analysis developed in econometrics, especially spatial regression modeling.  While spatial methods are frequently used in dEJ analysis, they are often not statistically sophisticated.  @ChakrabortyRevisitingToblerFirst2011 notes that dEJ regression analyses often use ordinary least squares (OLS), and fail to account for spatial autocorrelation.  This lack of statistical sophistication makes spatial EJ analyses vulnerable to technical criticism by "merchants of doubt" [@OreskesMerchantsDoubtHow2011], which can limit their effectiveness in pushing for policy change or legal remediation.  

It is important to recognize that environmental justice issues are not exhausted by the distribution of environmental hazards.  @SchlosbergDefiningEnvironmentalJustice2007, drawing on previous work by @YoungJusticePoliticsDifference1990 and @Shrader-FrechetteEnvironmentaljusticecreating2002, argues that environmental justice also includes procedural justice and appropriate recognition and respect for community identity.  For example, racialized communities that are outside of an administrative district — and so formally excluded from land-use decisions within the district — might be exposed to pollution emitted as a result of those land-use decisions [@LondonStruggleWaterJustice2018]; this is a form of procedural injustice.  Or, communities' claims and arguments might be ignored because they are racialized or lack formal scientific credentials [@OttingerBucketsResistanceStandards2010].  This is an example of misrecognition and disrespect.  

However, dEJ remains an important aspect of EJ, and the kinds of quantitative methods deployed in this project can be especially useful for identifying distributive environmental injustices.  

Previous spatial analyses of chlorpyrifos use and exposure in California fall into two categories.  First, physical-chemical simulation methods have been used to develop fate-and-transport estimates of chlorpyrifos presence across the entire state.  For example, @LuoSpatiallydistributedpesticide2010 use public chlorpyrifos use data and a fate-and-transport simulation to estimate how the chemical moves through space.  However, this study did not examine the population exposed to the pesticide, and therefore was not a dEJ analysis.  In contrast, this study focuses on demographic covariates for chlorpyrifos use in the Central Valley, and thus is more focused on the population that is potentially exposed than on the chemical itself.  

The other category of studies use relatively small-scale epidemiological methods to examine the public health impacts of chlorpyrifos exposure.  Several studies in this category have been conducted as part of the Center for the Health Assessment of Mothers and Children of Salinas (CHAMACOS) Study, based at University of California, Berkeley.  Over the past 20 years, the CHAMACOS Study has followed roughly 800 children in a farmworker community in California's Salinas Valley, a major agricultural region south of the San Francisco Bay [@CenterforEnvironmentalResearchandChildrensHealthCHAMACOSStudy].  @GunierPrenatalResidentialProximity2017 compare public pesticide use records to Wechsler Intelligence Scale for Children (WISC) scores for 255 7-year-old children.  Examining a 1 km buffer around the residence of pregnant women participants, they find that a 1 standard deviation increase in organophosphate use (including chlorpyrifos) in this buffer during pregnancy is associated with a 1-4 point decrease in WISC scores.  

The Columbia Center for Children's Environmental Health (CCEH) at Columbia University has also conducted epidemiological studies of the health impacts of chlorpyrifos that incorporate spatial data.  @LovasiChlorpyrifosExposureUrban2011 use neighborhood poverty and built environment indicators as controls in their analysis of chlorpyrifos impacts on cognitive and motor development in 266 3-year-old children.  This study focuses on children of African-American and Dominican mothers. 

Because these last two studies focus on populations that are likely to be highly exposed or socially vulnerable to chlorpyrifos impacts, they can be considered dEJ studies.  Both of these studies focus on estimating the health impacts of chlorpyrifos exposure, rather than the (relative or absolute) degree of exposure per se.  In contrast, the current paper considers social vulnerability as a predictor for (a proxy of) chlorpyrifos exposure.  The current study also works at a much larger scale than these two studies:  rather than analyzing data for a few hundred individuals, this study analyzes data for more than a thousand Census tracts and tens of thousands of uses of chlorpyrifos.  

Methodologically, the current study closely resembles two studies that used spatial methods to identify demographic predictors of potential exposure to other kinds of environmental health hazards.   @LievanosRacedeprivationimmigrant2015 uses data from across the continental US and spatial methods to construct clusters of high lifetime cancer risk (LCR) due to air pollution, then (nonspatially) regresses these clusters against composite Census tract demographic variables.  This study concludes that "isolated Latino immigrant-economic deprivation is the strongest positive demographic predictor of tract presence in air-toxic LCR clusters, followed by black-economic deprivation and isolated Asian/Pacific Islander immigrant-economic deprivation" [@LievanosRacedeprivationimmigrant2015 50], a significant dEJ finding.  

@SilvaSpatialModelingIdentify2018 examine the presence of hydraulic fracturing wells in Ohio Census block groups.  They compare spatial and nonspatial regression models, that is, models that do and do not account for correlation between neighboring block groups.  They observe that the nonspatial model is susceptible to spatial confounding in this case; and find that increased median income and population density are associated with an increased presence of wells.  



# Data #

## Study Area ##

The study area for this project is California's Central Valley.  The Central Valley can be defined in a number of different ways.  Since the units of analysis for this project are tracts and places designated by the US Census, a county-based definition was judged to be most appropriate.  Sacramento County was excluded because, unlike the rest of the region, most of its area is urban.  17 other counties were used to define the Central Valley; see table \ref{tab.counties} and figure \ref{fig.chlor_use}.  

||
|:-------|:------------|
| Shasta | San Joaquin |
| Tehama | Stanislaus  |
| Glenn  | Merced |
| Butte  | Madera |
| Colusa | Fresno |
| Yuba   | Kings  |
| Sutter | Tulare |
| Yolo   | Kern   |
| Solano ||
Table: California counties comprising the Central Valley for the purposes of this study.  Counties are listed roughly in north-south order.  \label{tab.counties}

## Pesticide Use ##

This paper combines two datasets.  First, California's Department of Pesticide Regulation (DPR) releases annual public datasets for agricultural pesticide use across the state, known as Pesticide Use Reports (PUR) [^DPR].  These datasets provide data on pesticide use at the township and section level, including pounds of active chemical used on particular days and times.  The most recent data release covers 2015.  These data exclusively report agriculture use; this limitation is acceptable for studying chlorpyrifos, which has been banned for residential use, in the agriculture-heavy Central Valley.  

[^DPR]: <http://www.cdpr.ca.gov/docs/pur/purmain.htm>

I retrieved the full datasets for 2011-2015.  Combining data over several years can help smooth annual variation in use due to variation in insect pest prevalence or other factors.  These datasets were combined and filtered on chemical name for chlorpyrifos.  Data for all counties across the state were used, not just the Central Valley, to avoid edge effects.  For example, while Sacramento County was not included in the study area, there is some chlorpyrifos use in the southeastern part of Sacramento County, very close to areas of San Joaquin County that are included in the study area.  

DPR conducts automated quality checks, and flags errors and changed values in the data.  I closely examined entries with "potential duplicate" and "outlier" error codes.  For the first few years of the dataset (2011 through 2013), it appears that DPR flagged but did not remove potential duplicates.  In later years, it appears that DPR removes potential duplicates in the released data.  I removed potential duplicates from the earlier years.  

DPR automatically identifies potential outliers as entries satisfying at least one of two criteria.  The first criterion is 200 lbs (solids and sprays) or 1000 lbs (fumigants) or more of active ingredient used per acre.  The second criterion is "50 times the median rate for all uses with the same pesticide product, crop treated, unit treated, and record type" (from the documentation file included in the data release).  Inspecting the entries with "potential outlier" flags, it appears that the main datafiles include corrected/estimated/modified values, not the original potential outliers.  I used the corrected/estimated/modified values given in the data without further adjustment.  

Chlorpyrifos uses are associated with townships and sections in the data, and I calculated annual and all-study-period totals of active ingredient at the centroid of each 1 mile-square ($1.6 km \times 1.6 km$) section.  Because these centroids do not match the actual use locations (i.e., farm fields), the centroid totals might be unreliable for the smallest CTD, 1km (see discussion and table \ref{tab.ctd} below).  However, this error should be negligible for the other CTDs.  

*[total num. observations, range of total use values]


## Demographics ##

The second primary dataset comprises American Community Survey (ACS) five-year estimates, from 2011-2015, for all Census tracts and places in the 17 Central Valley Counties.  

*[total num. tracts, places]

For each tract and place, I retrieved estimates and margin of error (MOE) values for three categories of demographic variables:  *race and ethnicity* (Hispanic, non-Hispanic White, non-Hispanic Black, Indigenous, and Asian residents), *foreign-born noncitizens*, *children under 5* (who may be especially sensitive to chlorpyrifos exposure due to small body weight and critical neurodevelopmental processes), and *poverty* (individuals with an Census-determined income-poverty ratio below 1).  Intersectionality theory [CollinsIntersectionality2013] suggests that interlocking systems of oppression, such as race-ethnicity and class, can produce synergistic effects, with greater disparities than either system would produce "on its own."  I therefore retrieved an intersectional Hispanic poverty variable (Hispanic individuals below the poverty line) calculated by the Census.  Because PUR data come only from agricultural uses, I predicted that chlorpyrifos use would be positively correlated with agricultural employment, and therefore retrieved an agricultural employment variable as a potential control.  (See also table \ref{tab.iv} and section \ref{sec.iv_selection}.)  

I calculated population densities and proportions (e.g., the fraction of all residents who are Hispanic) for the total population and each of these ACS variables for every tract and place, using Census-recommended methods to calculate MOEs for these derived variables [@USCensusBureauAmericanCommunitySurvey, 11ff].  9 tracts and 15 places with estimated total population or total employment of 0 were excluded from all analyses.  

Because much of the study area is rural, tract size and population density varied over four orders of magnitude, from 0.8 to 5600 residents per $km^2$.  Places cover 87% of the population, including 83% of non-Hispanic White and 88% of Hispanic residents, with population densities between 1.4 and 4400 residents per $km^2$.  However, places are geographically separated from each other, covering only 8% of the area of the tracts; in constructing contiguity-based spatial weights, 62% of places had no neighbors.  

The Modifiable Areal Unit Problem (MAUP) has been used to criticize dEJ projects [@SteelEnvironmentalJusticeValues2012;  @TaylorToxiccommunitiesenvironmental2014, 41ff].  "Egocentric neighborhood" methods have been used to address the MAUP in segregation research [@ReardonMeasuresSpatialSegregation2004].  However, this method assumes that populations are distributed evenly within each discrete region (e.g., each Census tract).  This assumption is inappropriate for this project, which includes many spatially heterogeneous rural regions.  In addition, the large geographic scale of this project would require millions — if not tens of millions — of egocentric neighborhoods, and so literally trillions of distance calculations between neighborhoods and chlorpyrifos uses.  The Census tract-level distance calculations already pressed the limits of the computing power available for this phase of the project.  More fine-grained regions (e.g., Census block groups or blocks) would have multiplied uncertainty in the ACS estimates, and also likely would have exceeded the available computing power.  

As a compromise, I used block population counts from the 2010 Census to calculate weighted centroids for each tract and place.  These centroids more accurately represent the "average location" of the population in each tract, without requiring more computing power in the distance calculation step.  By comparison, @LievanosRacedeprivationimmigrant2015 conducts a dEJ study of air toxic exposure across the US using Census tract centroids alone.[^NLCD]

[^NLCD]: An alternative approach would be to use the National Land Cover Database (NLCD), which was most recently released in 2011.  (<https://www.mrlc.gov/nlcd2011.php>)  Specifically, tract populations could be interpolated across the four categories of developed land identified in the NLCD.  Unfortunately, this public data source did not come to my attention until too late in the process of preparing this manuscript.  

Chlorpyrifos use section centroids, Census tracts, and places included in the study area are shown in figure \ref{fig.chlor_use}.  

![Data used in this study.  Red points are chlorpyrifos use totals, shown on a $log_{10}$ pounds scale and for this map 2015 only.  Blue regions are Census tracts included in the study area; yellow regions are included places.  All California counties are shown for context. \label{fig.chlor_use}](03_chlor_use.png)



## Linking Pesticide Use to Tracts and Places ##

Chlorpyrifos use is linked to tracts and places using the concept of *characteristic travel distance* (CTD).  CTD is defined as the distance "at which 63% of the original mass of [a] volatilized [chemical] is degraded or deposited" [@MackayFateenvironmentlongrange2014].  In the simplest sense, CTD models assume exponential decay:  

\begin{align}
    q(d) &= q_0 \beta^d \label{eqn.decay}
\end{align}

where $q(d)$ is the quantity of a chemical at a distance $d$ from its source; $q_0$ is the quantity volatilized at the source; and $\beta$ is a constant related to the CTD:  

\begin{align} 
    (1-0.63) &= \beta^{CTD}\\
    \beta &= 0.37^{1/CTD}
\end{align}

I will refer to $\beta^d$ as the *decay coefficient*.  Given a CTD, the decay coefficient can be used in equation \ref{eqn.decay} to calculate the *distance-weighted local use* of chlorpyrifos at each location (tract or place centroid) $i$:  
\begin{align}
    q_i &= \sum_u q_u \beta^{d_{iu}} \label{eqn.tot_pot_exp}\\
        &= \sum_u q_u 0.37^{d_{iu}/CTD}
\end{align}
where $q_u$ is the total use at section $u$ and $d_{iu}$ is the Euclidean distance between the centroid of $i$ and the centroid of $u$.  To slightly account for the fact that the residents of a location are not located at its centroid, the decay coefficient is set to 1 whenever the centroid associated with use $u$ is within $i$, regardless of $d_{iu}$.  

It is important to stress that equations \ref{eqn.decay} and \ref{eqn.tot_pot_exp} are, at best, crude estimate of potential exposure.  They do not take into account prevailing or occurrent winds, or other chemical transport mechanisms.  Read as exposure estimates, they represent the chemical moving uniformly away from the source along a one-dimensional path, not dispersing over a two- or three-dimensional space surrounding the source at various speeds in various directions at various times.  diverse and variable processes of application, fixation, and chemical transformation are represented as simple exponential decay.  I therefore refer to these aggregate statistics as "local use" rather than "exposure."  By comparison, @LuoSpatiallydistributedpesticide2010 use a physical-chemical model and PUR data to produce much more sophisticated estimates of chlorpyrifos loading.  

@MackayFateenvironmentlongrange2014 estimate the CTD for chlorpyrifos to be 62 km.  For robustness, this paper considers five CTD values, ranging from 1 km to 90 km.  See table \ref{tab.ctd} and \ref{fig.ctd}.  

| CTD (km) | $\beta$ |
|:---------|:--------|
| 1        | 0.370   |
| 10       | 0.905   |
| 30       | 0.967   |
| 60       | 0.984   |
| 90       | 0.989   |

Table: Characteristic Travel Distance (CTD) values used in this study, and corresponding decay-rate values $\beta$. \label{tab.ctd}

![Impact of Characteristic Travel Distance (CTD) value on weighted local use values in Census tracts. Panels correspond to the different CTD values used in this study. Weighted local use is the log (base 10) of aggregate chlorpyrifos use around each tract, scaled using the decay coefficient, 2011-15. Color scales are only roughly consistent between panels, with the scale midpoint set at 5 ($10^5 = 100,000$ lbs).](07_ctd.png)


# Methods #

The primary analysis of this study is a regression of relative potential chlorpyrifos exposure against Census demographic data.  Separate models are constructed for each of the five CTD values listed in table \ref{tab.ctd}, as well as for tracts and places.  These two methodological choices give $5 \times 2 = 10$ models.  The software language R was used to clean and analyze all data, with significance use of the `tidyverse`, `tidycensus`, `sf`, `spdep`, and `tmap` packages [@WickhamtidyverseEasilyInstall2017; @WalkertidycensusLoadUS2018; @PebesmasfSimpleFeatures2018; @BivandspdepSpatialDependence2018; @TennekestmapThematicMaps2018].  Complete cleaning and analysis code is available at <https://github.com/dhicks/spatial_project>.  

In the terminology of @BogenSavingPhenomena1988, this project aims to characterize a phenomenon rather than confirm or disconfirm a theory.  Bogen and Woodward argue that phenomena are intermediate between data and theories.  Data can be thought of as the more-or-less direct outputs of instruments — in the case of this project, the PUR forms filled out by farmers and survey instruments used by the Census.  Phenomena are more-or-less stable patterns or trends in the world, and are presumed to reflect "signal" rather than "noise" in the data.  But phenomena are not explanations of the patterns or trends.  Rather, explanations are given by theories (in their terminology), which offer models of the causal processes that produce phenomena.  The distinction between phenomena and theories is close to the distinction between correlation and causation:  phenomena describe correlations, while theories aim at causation.  Thus, this study aims to describe phenomenal relationships between chlorpyrifos use and population demographics, but not to explain them.  

Note that the theory-phenomenon distinction does not mean that phenomena are theory-independent.  I have already introduced theoretical assumptions about, for instance, how chlorpyrifos decays over space, and which Census variables might be associated with chlorpyrifos use.  More generally, separating "signal" and "noise" in statistical analysis requires background assumptions about what kinds of processes might be at work in the data-generating process.  But this project merely relies on those assumptions; it does not intend to confirm or deny them.  

Arguably, the theory-phenomenon distinction aligns with the difference between procedural and distributive aspects of EJ.  The unjust distribution of environmental hazards is the outcome of social processes, which are often unjust in themselves.  


## Independent Variable Selection ## \label{sec.iv_selection}

Non-spatial exploratory data analysis indicated that, for almost all (>80%) tracts and places, almost all residents (>80%) were either Hispanic or non-Hispanic White.  Thus, Hispanic and non-Hispanic White proportions are strongly negatively correlated ($r = -.9$), and I dropped non-Hispanic White proportion from the independent variable list.  Similarly, there are strong correlations between Hispanic and non-citizen proportions ($r = .8$), and so we also drop non-citizen proportion from the independent variable list.  (Alternatively, factor analysis methods might have been used to construct composite variables combining Hispanic and non-citizen proportions; compare @LievanosRacedeprivationimmigrant2015.  This approach was not taken here for simplicity.)  There are moderate correlations ($r = .4-.6$) between Hispanic, young children, and poverty proportions; so, while we use all three of these variables, we expect to see greater uncertainty in their coefficient estimates, due to the difficulty of disentangling their statistical effects [@MorrisseyMultipleRegressionNot2018].  

I also considered including an intersectional, Hispanic poverty proportion in the independent variable list.  However, taking the proportion with respect to the total population (Hispanic poverty / total population) yielded a very strong correlation with the Hispanic proportion (Hispanic / total population; $r = .9$); and taking the proportion with respect to the Hispanic population (Hispanic poverty / Hispanic) yielded a very strong correlation with general poverty (poverty / total population; $r = .8$).  I therefore did not include this intersectional variable here.  

Also for intersectional reasons, I considered including interaction terms between Hispanic, young children, and poverty proportions.  I examined the effect of including various interaction terms on coefficient variance inflation factors (VIF) for a non-spatial linear regression model on Census-designated places.  Standard statistical practice suggests that VIFs greater than 10 are concerning.  All of these interaction terms produced VIF > 10 for at least 1 coefficient in either tracts or places (or both), and so were not included.  

Finally, I also include two control variables.  Because PUR data come only from agricultural uses, we expect chlorpyrifos use to be negatively correlated with population density and positively correlated with agricultural employment.  Density was not substantially correlated with any variables of interest (except for a low-moderate negative correlation with White; $r = -.4$).  Because density is left-bounded at 0, we use log density to simplify the Bayesian prior on its true value.  Agricultural employment has a moderate-strong correlation with Hispanic ($r = .67$).  

All independent variables are given in table \ref{tab.iv}.  

| Category       | Independent Variable |
|:---------------|:-------------------|
| Race-Ethnicity | Hispanic           |
|                | Black              |
|                | Indigenous         |
|                | Asian              |
| Children       | Children under 5   |
| Class          | Income-Poverty Ratio < 1.0 |
| Controls       | log Population Density  |
|                | Employed in Agriculture | 

Table: Independent variables (IVs) used in this study. All race-ethnicity groups other than Hispanic are non-Hispanic.  All IVs other than population density are proportion of total population in the tract or place. \label{tab.iv}

I took an effects estimation approach to these variables, rather than a hypothesis testing approach.  That is, I did not test the null hypothesis that the effect of any given variable on local use is exactly 0.  In line with the aphorism that "everything is influenced by everything else," in this kind of observational/field research setting it is a priori highly implausible that the effect is exactly 0.  Hypothesis testing therefore merely tells us whether we have enough data to detect the non-zero effect, which is not informative.  In addition, "stepwise" regression model selection methods (which include or exclude variables based on whether their coefficients are statistically significant) are known to introduce bias and interpretability problems for both hypothesis testing and estimation [@gungAlgorithmsautomaticmodel2012].  

Instead of testing hypotheses, I aimed to *[estimate effects as accurately possibly]


## Spatial Weights ##

Exploratory data analysis considered 8 spatial weight constructions for both tracts and places:  contiguity, inverse distance weights (with an outer limit of 50 km and a decay of $\frac{1}{d}$), and $k$-nearest-neighbors (KNN) with $k$ ranging from 3 to 8.  All weights were row-normalized.  As noted above, 62% of places had no contiguity-based neighbors.  

To examine the impact of different spatial weights constructions, Moran's $I$ was calculated for population densities corresponding to the independent variables identified above, e.g., density of Hispanic population, calculated as the number of Hispanic residents per square kilometer.  

Among tracts, contiguity produced the highest values of Moran's $I$, followed by KNN, and distance weights.  For example, for Hispanic population density, Moran's $I$ was slightly less than .7 for contiguity weights; was between .45-.55 for KNN; and was slightly greater than .3 for distance weights.  

Among places, KNN and distance-based weights were similar, especially for larger values of K, while contiguity-based weights were much less for most variables.  For example for Hispanic population density, Moran's I was between .55-.65 for KNN; about .53 for distance weights; but only .3 for contiguity weights.  

KNN weights were selected for use in further analysis.  They produced reasonably high and consistent values of Moran's $I$ for both tracts and places.  For the block bootstrap procedure used to account for measurement error (discussed below), KNN blocks were dichotomous (for a given location $i$, every other location $j$ is either one of $i$'s $k$ nearest neighbors or not), had uniform size $k+1$, and could be characterized exactly by two values (namely, $k$ and the location $i$).  

Initially, the block bootstrap procedure was written to allow exploring several values of $k$.  However, runtime issues required reducing the number of free parameters explored in this stage.  Again, EDA indicated that values of $k$ from 3 to 8 produced roughly the same spatial correlation values, while CTD values produced large differences between dependent variable values.  I therefore judged that different CTD values were likely to be a more important source of variations in effects than different spatial weights.  (See also @LeSageBiggestMythSpatial2014.)  I selected $k=3$ as the sole spatial weights construction for all analyses.  


## Regression Specification ##

Exploratory data analysis indicated that local chlorpyrifos use values are non-Gaussian.  For lower CTD values (1, 10), the left bound at 0 creates highly asymmetrical distributions.  A log transformation of the DV appeared to mitigate the left bound for the lower CTD values without qualitatively changing the distributions for the higher CTD values.  I therefore use a log-transformed dependent variable.  Base 10 was used for interpretability, so that dependent variable units are locally-weighted pound orders of magnitude.  

For higher CTD values (30, 60, 90), distributions are bi- or trimodal.  Plotting separate distributions for each county suggested that this was due to very different county-level baselines.  For counties with many tracts or places, the distribution of values appeared to be sufficiently Gaussian for standard regression.  County-level dummy variables were considered, and were found to improve homoscedasticity,  but also introduced perfect multicollinearity in lagged models.  However, inspection of residuals in the spatial Durbin models suggested that lagged dependent variables substantially improved homoscedasticity without county-level dummies.  Still, even the spatial Durbin models exhibit some heteroscedasticity.  

Spatial exploratory data analysis of both independent and dependent variables suggested substantial degrees of spatial autocorrelation on both sides of the regression formula.  A sequence of three model specifications was considered:  "vanilla" linear regression, without any spatial component (referred to below as "regression" below); spatial regression with lagged independent variables, or "spatial lag X"; and spatial Durbin regression, which incorporates lags for both dependent and independent variables [@LeSageIntroductionSpatialEconometrics2009].  If the regression model is specified as 
\begin{align}
    Y &= \alpha 1_n + X \beta + \varepsilon,
    \varepsilon &\sim N(0, \sigma^2)
\end{align}
where $\alpha$ is a scalar parameter, $1_n$ is a $n$ column vector of 1s, $X$ is a $n \times p$ design matrix, and $\beta$ is a length-$p$ column vector of regression coefficients, then the spatial lag X model is specified as
\begin{align}
    Y &= \alpha 1_n + X \beta + WX \theta + \varepsilon
\end{align}
where $W$ is the $n \times n$ spatial weights matrix and $\theta$ is a length-$p$ column vector of coefficients.  The spatial Durbin model is further specified as
\begin{align}
    Y &= \alpha 1_n + X\beta + WX\theta + \rho WY + \varepsilon\\
    Y &= (I_n - \rho W)^{-1} (\alpha 1_n + X\beta + WX\theta)
\end{align}
where $\rho$ is a scalar parameter and $I_n$ is the $n \times n$ identity matrix.  

KNN spatial weights, with $k=3$, were used for all of these models.  This sequence of models was fitted for both places and tracts, for CTD values of 10 and 60.  Residual plots were examined for indications of heteroscedasticity, and $R^2$, AIC, and Moran's $I$ of the residuals were compared within model-dataset combinations.  


## Measurement Error and Block Bootstrap ##

All three of these models assume perfect measurements in the independent variables.  However, ACS estimates can have large margins of error, especially for subpopulations of rural tracts \*[cite].  This kind of measurement error induces *attenuation bias* or *regression dilution*, in which correlation estimates are shrunk towards 0 [@FrostCorrectingregressiondilution2000].  In the context of dEJ analysis, attenuation bias is potentially a serious problem, insofar as it leads to the underestimation of environmental disparities.  That is, measurement error can make distributive environmental injustices seem less serious than they actually are.  

A block bootstrap approach was developed to explore the effect of measurement error on the impact estimates.  In a "basic" or unparameterized bootstrap, independent samples are taken from the observations in the dataset (with replacement), forming a "resampled dataset" of the same size as the dataset.  These samples approximate drawing a new sample from the original population.  By calculating model statistics on a set of resampled datasets (1000 resampled datasets is common), we can estimate features of the sampling distribution of the population statistic.  For example, if $\hat\beta^1, \hat\beta^2, \ldots, \hat\beta^{1000}$ are regression coefficient estimates calculated on 1,000 resampled datasets, $\mathrm{sd}_l(\hat\beta^l)$ estimates the standard error of the coefficient estimate. 

However, because spatial and temporal observations are not independent, this simple sampling method cannot be used in these cases.  The *block bootstrap* approach has been developed for these kinds of data [@LahiriResamplingMethodsDependent2003].  Observations are arranged into blocks, and blocks are sampled rather than individual observations to form the resampled datasets.  

Specifically, a moving block bootstrap approach designed for raster data was modified for the present case [@LahiriResamplingMethodsDependent2003, ch 12].  Each block was defined as a given location (place or tract) and its $k$ nearest neighbors (with $k=3$ selected as above).  Since each block contained the same number of observations (namely, $k+1=4=n_b$),  $\lfloor{\frac{n}{n_b}}\rfloor$ blocks were drawn to produce a resampled dataset of approximately the same size as the original.  The spatial weights matrix for the resample was constructed by subsetting the original spatial weights (repeating rows/columns as individual locations appeared multiple times in the data) and re-row normalizing.  

500 such resampled datasets were constructed for each combination of dataset (places or tracts) and CTD value.  With two datasets and 5 CTD values, 10 top-level "models" were examined.  For each resampled dataset, a spatial Durbin model was fit and the following statistics were calculated: IV impacts for all covariates; $\rho$, the coefficient on the lagged dependent variables; and Moran's $I$ on the model residuals.  Impacts were calculated using a Monte Carlo method, with 100 samples for each resampled dataset.  These Monte Carlo draws were then combined at the dataset-CTD level, producing $100 \times 500 = 50,000$ samples for the impacts for each IV.  

It is important to note that this block bootstrap approach has not been validated, either formally or with simulated values.  It is based on an approach designed for raster data that are used to model a smoothly-varying spatial field, not regression models fit to lattice data.  I include it here as the best available way to attempt to analyze the effects of measurement error and attenuation bias.  


# Results #

## Model Selection ##

For each of 4 dataset-CTD combinations examined in model selection, spatial Durbin models consistently outperformed both regression and spatial lagged X models.  AIC and Moran's $I$ on the residuals were lowest with spatial Durbin models in every case.  $R^2$ was generally .6-.65 for the spatial Durbin models on the places (except for CTD = 1, where $R^2 > .75), compared to values around .4 for regression and .55 for spatial lagged X models.  For tracts, $R^2$ was near .5 for spatial Durbin models, in each case slightly greater than approximately .45 for the spatial lagged X models.  ($R^2$ for the spatial Durbin models was calculated using "trend" predictions, which ignore the lagged dependent variable.  When the lagged dependent variable was included, $R^2$ was 1.00 after rounding to two decimal places.)  These three statistics all indicate that spatial Durbin models provide a better fit.  However, also in every case, Moran's $I$ was substantially (and statistically significantly) greater than 0 even for spatial Durbin models, with values of approximately $.07$ for tracts and $.15$ for places.  This suggests that there may be spatial non-stationarity; that is, the effects of the independent variables may vary across different sub-regions in the study area.  This possibility was not explored further in the current study.  

Residual plots had strong signs of heteroscedasticity for both regression and spatial lagged X.  In general, residual variances were larger at small fitted values.  County-level heterogeneity appeared to be responsible for this pattern.  In tracts, with spatial Durbin models, any lingering indication of heteroscedasticity was due to 3 (out of 1044) locations with small fitted values and large (negative) residuals.  These locations can plausibly be treated as outliers, and so there do not appear to be any issues of heteroscedasticity at the tract level.  Places showed more heteroscedasticity, with approximately 15 (out of 391) locations with small fitted values and large (negative and positive) residuals in the spatial Durbin models.  This heteroscedasticity primarily means that the coefficient estimates are even more uncertain (e.g., larger standard errors) than are reported below.  

## Moran's $I$ and $\rho$ ##

Violin plots showing the block bootstrap distribution of Moran's $I$ and $\rho$, the coefficient on the lagged local use values, are shown in figures \ref{fig.moran} and \ref{fig.rho}, respectively.  

![Block bootstrap distributions of Moran's $I$ for the residuals for each dataset-CTD combination.  Distribution width corresponds to the number of resamples at that value of $I$.  Horizontal lines indicate 5%, 50% (median), and 95% quantiles. Underlying data are 500 block bootstrap resamples at each dataset-CTD combination.  Points indicate estimated values for spatial Durbin models fitted on the complete original dataset.  \label{fig.moran}]{11_moran.png}
![Block bootstrap distributions of $\rho$, coefficient on lagged dependent variables, for each dataset-CTD combination.  Distribution width corresponds to the number of resamples at that value of $\rho$.  Horizontal lines indicate 5%, 50% (median), and 95% quantiles. Underlying data are 500 block bootstrap resamples at each dataset-CTD combination. Points indicate estimated values for spatial Durbin models fitted on the complete original dataset.  \label{fig.rho}]{11_rho.png}

For both statistics, tracts and places are comparable given a CTD value.  Moran's $I$ value tends to increase with CTD, while $\rho$ tends to decrease.  Scatterplots comparing $\rho$ to Moran's $I$ indicated that, for most dataset-CTD combinations, there was very strong negative correlation between these two statistics.  (This trend was weaker for places with CTD 1km, and the association was weakly positive for tracts with CTD 1km.)  

Points in figures \ref{fig.moran} and \ref{fig.rho} show the point estimates calculated for the spatial Durbin models fitted on the complete original dataset.  Notably, for larger CTD values, the estimates for the complete models differ dramatically from those for the block bootstrap resamples.  For example, $\rho$ is high for these models, near 1, while the block bootstrap remodels tend to have low values of $\rho$.  

At a theoretical level, spatial patterns in the dependent variable that are not associated with the independent variables need to be partitioned between the autoregression coefficient $\rho$ and the autocorrelation structure of the residuals measured by Moran's $I$.  The block bootstrap resamples, as it were, explore different ways of making this partition.  It may be that there is no theoretically correct way of making this partition, and that the complete models have simply happened to settle on especially extreme partitions, with no deeper implications to be drawn.  Or it may be that theoretical or coding problems with the block bootstrap approach produce the "exploration" process as an artifact or error.  Again, the particular block bootstrap approach used here has not been validated.  

## IV Impacts ##

Unlike "vanilla" linear regression models, which treat observations as formally and mathematically independent, spatial models treat observations — locations — as connected.  This means that changes in an IV at one location can influence the DV at another location, corresponding to the term $WX\theta$.  Further, the spatial Durbin model's lagged dependent variable term, $WY\rho$, introduces the possibility for feedback loops:  a change IV $\Delta x_i$ at location $l$ induces a change $\Delta y_{l'} = w_{l'l}\theta \Delta x_i$ in $y$ at neighbor $l'$, which feeds back to location $l$ as $w_{ll'} \Delta y_{l'} \rho$.  (This and the next paragraph generally follow @LeSageIntroductionSpatialEconometrics2009, §2.7.)  

Spatial econometricians have introduced the notion of *impacts* for the interpretation of regression coefficients under spatial feedbacks.  In "vanilla" linear regression models without interaction, the coefficient $\beta_i$ for IV $x_i$ is identical to the partial derivative $\partial y / \partial x_i$.  $\beta_i$ can therefore be interpreted directly as the marginal effect of $x_i$ on $y$.  But in the spatial Durbin model, the partial derivative 
\[ \frac{\partial y}{\partial x_i} = (I_n - W\rho)^{-1}(I_n \beta_i + W \theta_i = S_i(W) \]
depends not just on the coefficients $\beta_i$ and $\theta_i$, but also the autoregression coefficient $\rho$.  And the value of this partial derivative at location $l$ depends on its connection to other locations, as encoded in $W$.  The *total impacts* for IV $x_i$ are formally defined as the mean row sum of $S_i(W)$, which corresponds to averaging $S_i(W)$ across locations.  These impacts are efficiently calculated using the \texttt{impacts()} function in the R package \texttt{spdep}, which uses a Monte Carlo method based on the traces of the powers of $W$.  The Monte Carlo method produces a sequence of estimates for total impact, which can be analyzed as an estimate of the sampling distribution for the true total impact.  

Figures \ref{fig.impacts_1}, \ref{fig.impacts_10}, \ref{fig.impacts_369} show total impact Monte Carlo estimates for IVs for CTD values of 1, 10, and 30-60-90 respectively. In the present study, total impact estimates were made for both the spatial Durbin models fit on the entire complete dataset and for each block bootstrap resample.  Resample estimates were then combined into a single estimated sampling distribution for each dataset-CTD-IV combination.  

![Total IV impacts, CTD = 1 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_1}](11_impacts_1.png)

![Total IV impacts, CTD = 10 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_10}](11_impacts_10.png)

![Total IV impacts, CTD = 30, 60, and 90 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_369}](11_impacts_369.png)

Comparing the three plots, the estimates tends to be much more uncertain for smaller CTD values.  Indeed, for CTD = 1 km, several estimates are uncertain across more than 10 orders of magnitude (evaluating uncertainty by the width of the 5%-95% quantile interval).  By contrast, for CTD values of 30 or greater, uncertainty for almost all variables is less than about 2 orders of magnitude, and in some cases less than 1 order of magnitude.  In short, whatever physical validity different CTD values have, higher CTD values lead to more precise effects estimates.  The largest uncertainties are for the smallest subpopulations in the study area, namely, Asian, Black, Indigenous, and children proportions.  

The figures show both the block bootstrap resamples (violin plots) and estimates from the complete original dataset (point-range plots).  For larger CTD values, the estimates for the complete model can be quite different from the resample estimates.  For example, for tracts, the complete model estimates for the effects of agricultural employment and population density are well above the 90% quantiles from the resamples; while the estimates for Hispanic proportion are well below those of the resamples.  In many cases the resample distributions are more conservative (closer to zero); but for Hispanic and Indigenous proportion the complete model estimates are more conservative.  It may be that, just as the resamples effectively explored ways of partitioning spatial patterns in the dependent variable between autoregression and autocorrelation, they also effectively explored ways of partitioning effects estimates between partially correlated independent variables.  In any case, the block bootstrap approach was adopted to examine the potential effect of attenuation bias, but it is not clear whether this was successful.  For this reason, the discussion below focuses on the estimated impacts from the full model.  

![Total impact estimates from the spatial Durbin models for CTD = 60 km. Total impact estimates $\hat\zeta$ are transformed as $\zeta_{trans} = 10^{\zeta/10}$.  $\zeta_{trans}$ can be interpreted as the relative change in local use when the corresponding IV increases by 10%.  For example, if $\zeta_{trans} = 1.5$, then local use is 50% greater when the corresponding IV is 10% greater.  Transformed values $>1$ therefore correspond to increases; values $<1$ correspond to decreases.  Line ranges give 90% confidence intervals. All estimates are based on the complete original datasets, not the block resamples.  \label{fig.impacts_backtrans}](11_impacts_backtrans.png)

|geography |CTD |IV                 | est. (transformed)| estimate| 95% CI|      |
|:---------|:---|:------------------|------------------:|--------:|------:|-----:|
|places    |60  |ag. employment     |               0.99|    -0.04|  -0.54|  0.49|
|tracts    |60  |ag. employment     |               1.14|     0.56|   0.27|  0.86|
|places    |60  |Asian              |               0.99|    -0.05|  -2.21|  1.10|
|tracts    |60  |Asian              |               1.01|     0.05|  -0.33|  0.36|
|places    |60  |Black              |               0.67|    -1.75|  -3.69| -0.45|
|tracts    |60  |Black              |               0.89|    -0.52|  -0.97| -0.14|
|places    |60  |children           |               0.95|    -0.24|  -3.14|  2.16|
|tracts    |60  |children           |               1.14|     0.56|  -1.22|  1.86|
|places    |60  |Hispanic           |               1.08|     0.34|  -0.03|  0.69|
|tracts    |60  |Hispanic           |               1.07|     0.30|   0.10|  0.50|
|places    |60  |Indigenous         |               0.66|    -1.82|  -4.30|  0.06|
|tracts    |60  |Indigenous         |               0.68|    -1.70|  -4.51|  1.00|
|places    |60  |pop. density (log) |               1.05|     0.20|   0.08|  0.33|
|tracts    |60  |pop. density (log) |               1.06|     0.27|   0.23|  0.30|
|places    |60  |poverty            |               1.06|     0.25|  -0.28|  0.76|
|tracts    |60  |poverty            |               0.94|    -0.26|  -0.52|  0.02|

Table: Estimates of total impact (direct + indirect) from spatial Durbin models for CTD = 60km.  *IV*:  Independent variable.  *Est. (transformed)*: Point estimate for total impact, transformed as $10^{\zeta/10}$ to aid interpretation. \label{tab.impacts}


# Discussion #

I focus this discussion on the total impact estimates for Hispanic, poverty, and agricultural employment proportion, and population density, for a CTD of 60 km.  These estimates are reported in figure \ref{fig.impacts_backtrans} and table \ref{tab.impacts}.  Transformed estimates are reported to aid interpretation:  when the estimated value of $\zeta$ is transformed as $\zeta_{trans} = 10^{\zeta/10}$, $\zeta_{trans}$ can be interpreted as the multiplicative change in local chlorpyrifos use when the corresponding IV increases by 10%.  For example, if $\zeta_{trans} = 1.5$ for Hispanic proportion, then a 10-point increase in Hispanic proportion is associated with a 50% in local chlorpyrifos use.  

*[write discussion]

- why CTD = 60?  
    - good precision in estimates - comparable to 90, better than 30
    - more conservative (closer to 0) for these variables
    - supported by literature

- ag employment and poverty: ambiguous
- population density and Hispanic: evidence of positive effect












Limitations
- PUR limitations
    - agricultural use only, not industrial, commercial, or state use (or residential, but there's no legal residential use of chlorpyrifos)
- local use, not exposure
    - but if practices are similar across amounts used, it's plausible that residential exposure will vary linearly with local use
    - and just one chemical; mixtures
    - roughly models residential exposure, not occupational, children at school, take-home occupational
    - only transport through air
- spatial distribution of population
    - land use data could be used to refine centroid locations
    - more computing power would be needed for "egocentric neighborhoods"
    - or raster methods?  
- correlated IVs
- block bootstrap approach - validate? getting at something other than attenuation bias? 
- heteroscedasticity suggests spatial heterogeneity
    - crop species
    - suggestions from Noli
- phenomena, not theory
- only distributive


Policy implications
- buffers and label guidelines
- 


# References #
