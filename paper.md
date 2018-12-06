---
title: Environmental Justice Analysis of Chlorpyrifos Use in California's Central Valley
author: Daniel J. Hicks
bibliography: spatial_project.bib
numbersections: true
---

*[limitations]*

# Abstract #

## Background ##

Chlorpyrifos is one of the most widely-used pesticides in the world, and is generally recognized to be a moderate neurotoxin.  

## Objectives ##

This paper conducts an distributional environmental justice (dEJ) analysis of chlorpyrifos use in California's Central Valley.  A dEJ analysis examines the way the distribution of environmental risks are associated with race, ethnicity, class, gender, and other systems of structural oppression.

## Methods ##

Data on chlorpyrifos use in townships and sections are retrieved from California's Department of Pesticide Registration public pesticide use records for 2011-2015.  These data are combined with demographic data on Census tracts and Census-designated places in the Central Valley from the American Community Survey (ACS).  Spatial and non-spatial regression models are used to estimate effects of demographic covariates on local chlorpyrifos use.  A bootstrap method is used to account for measurement error in the ACS estimates.  

## Results ##

Model selection statistics consistently favored spatial Durbin models, which are used in further analysis.  *[The block bootstrap results are difficult to interpret, and it is unclear whether this method has avoided attenuation bias.]* 

## Discussion ## 

*[rewrite; move to Results and replace w/ summary of Discussion?]*
There is consistent evidence that a 10-point increase in Hispanic population is associated with a 7-8% increase in potential chlorpyrifos exposure.  The effects of agricultural employment and poverty proportions are ambiguous.  Due to significant uncertainty in effects estimates, it is unclear if there is any association between potential chlorpyrifos exposure and Asian, Indigenous, and young children population proportion. 


# Introduction #

Chlorpyrifos is one of the most widely-used pesticides in the world. In California in 2016, it was the 29th most heavily used pesticide active ingredient, with over 900,000 pounds of active ingredient applied over 640,000 acres [@PesticideUseReporting2017].  It is an organophosphate (OP) pesticide, is generally recognized to be a moderate neurotoxin [@WorldHealthOrganizationWHORecommendedClassification2010 67; @USEPAChlorpyrifosRevisedHuman2016; citations in @TrasandeWhenenoughdata2017], and was banned for residential use in the US in 2001.  Based on the toxicological and epidemiological evidence, in 2007 the environmental organizations Pesticide Action Network North America (PANNA) and Natural Resources Defense Council (NRDC) filed a petition with US EPA, calling on the agency to revoke all tolerances for chlorpyrifos.  In 2017, US EPA rejected this petition [@USEPAChlorpyrifosOrderDenying2017], and as of this writing it remains in widespread use in California and elsewhere. 

This paper conducts an distributional environmental justice (dEJ) analysis of chlorpyrifos use in California's Central Valley.  A dEJ analysis examines the way the distribution of environmental risks intersect with race, ethnicity, class, gender, and other systems of structural oppression.  Since the landmark report "Toxic Wastes and Race in the United States"  [@CommissionforRacialJusticeToxicWastesRace1987], a significant dEJ literature has emerged [@PulidoEnvironmentalismEconomicJustice1996; @Shrader-FrechetteEnvironmentaljusticecreating2002; @BrownToxicExposuresContested2007; @OttingerTechnoscienceenvironmentaljustice2011; @TaylorToxiccommunitiesenvironmental2014].  

This paper uses methods of spatial analysis developed in econometrics, especially spatial regression modeling.  While spatial methods are frequently used in dEJ analysis, they are often not statistically sophisticated.  @ChakrabortyRevisitingToblerFirst2011 notes that dEJ regression analyses often use non-spatial methods that fail to account for spatial autocorrelation.  This lack of statistical sophistication makes spatial EJ analyses vulnerable to technical criticism by "merchants of doubt" [@OreskesMerchantsDoubtHow2011], which can limit their effectiveness as tools for policy change or legal remediation.  

It is important to recognize that environmental justice issues are not exhausted by the distribution of environmental hazards.  @SchlosbergDefiningEnvironmentalJustice2007, drawing on previous work by @YoungJusticePoliticsDifference1990 and @Shrader-FrechetteEnvironmentaljusticecreating2002, argues that environmental justice also includes procedural justice and appropriate recognition and respect for community identity.  For example, racialized communities that are outside of an administrative district — and so formally excluded from land-use decisions within the district — might be exposed to pollution emitted as a result of those land-use decisions [@LondonStruggleWaterJustice2018]; this is a form of procedural injustice.  Or, communities' claims and arguments might be ignored because they are racialized or lack formal scientific credentials [@OttingerBucketsResistanceStandards2010].  This is an example of misrecognition and disrespect.  

However, dEJ remains an important aspect of EJ, and the kinds of quantitative methods deployed in this project can be especially useful for identifying distributive environmental injustices.  

Previous spatial analyses of chlorpyrifos use and exposure in California fall into two categories.  First, physical-chemical simulation methods have been used to develop fate-and-transport estimates of chlorpyrifos presence across the entire state.  For example, @LuoSpatiallydistributedpesticide2010 use public chlorpyrifos use data and a fate-and-transport simulation to estimate how the chemical moves through space.  However, this study did not examine the population exposed to the pesticide, and therefore was not a dEJ analysis.  In contrast, this study focuses on demographic covariates for chlorpyrifos use in the Central Valley, and thus is more focused on the population that is potentially exposed than on the chemical itself.  

The other category of studies use comparatively small-scale epidemiological methods to examine the public health impacts of chlorpyrifos exposure.  Several studies in this category have been conducted as part of the Center for the Health Assessment of Mothers and Children of Salinas (CHAMACOS) Study, based at University of California, Berkeley.  Over the past 20 years, the CHAMACOS Study has followed roughly 800 children in a farmworker community in California's Salinas Valley, a major agricultural region south of the San Francisco Bay [@CenterforEnvironmentalResearchandChildrensHealthCHAMACOSStudy].  @GunierPrenatalResidentialProximity2017 compare public pesticide use records to Wechsler Intelligence Scale for Children (WISC) scores for 255 7-year-old children.  Examining a 1 km buffer around the residence of pregnant women participants, they find that a 1 standard deviation increase in organophosphate use (including chlorpyrifos) in this buffer during pregnancy is associated with a 1-4 point decrease in WISC scores.  

The Columbia Center for Children's Environmental Health (CCEH) at Columbia University has also conducted epidemiological studies of the health impacts of chlorpyrifos that incorporate spatial data.  @LovasiChlorpyrifosExposureUrban2011 use neighborhood poverty and built environment indicators as controls in their analysis of chlorpyrifos impacts on cognitive and motor development in 266 3-year-old children.  This study focuses on children of African-American and Dominican mothers. 

Because these last two studies focus on populations that are likely to be highly exposed or socially vulnerable to chlorpyrifos impacts, they can be considered dEJ studies.  Both of these studies focus on estimating the health impacts of chlorpyrifos exposure, rather than the (relative or absolute) degree of exposure per se.  In contrast, the current paper considers social vulnerability as a predictor for potential chlorpyrifos exposure.  The current study also works at a much larger scale than these two studies:  rather than analyzing data for a few hundred individuals, this study analyzes data for more than a thousand Census tracts and tens of thousands of uses of chlorpyrifos across more than 10,000 square miles.  

Methodologically, the current study closely resembles a number of other studies that use spatial methods to identify demographic predictors of potential exposure to other kinds of environmental health hazards.  Often these studies are framed explicitly in terms of environmental justice.  @LievanosRacedeprivationimmigrant2015 uses data from across the continental US and spatial methods to construct clusters of high lifetime cancer risk (LCR) due to air pollution, then (non-spatially) regresses these clusters against composite Census tract demographic variables.  This study concludes that "isolated Latino immigrant-economic deprivation is the strongest positive demographic predictor of tract presence in air-toxic LCR clusters, followed by black-economic deprivation and isolated Asian/Pacific Islander immigrant-economic deprivation" [@LievanosRacedeprivationimmigrant2015 50], a significant dEJ finding.  (See also @LievanosSociospatialDimensionsWater2017; @LievanosRetoolingCalEnviroScreenCumulative2018.)  

@GrineskiAsianAmericansdisproportionate2017 also focus on demographic predictors of air pollution exposure.  They define non-spatial clusters of Census tracts based on median year of housing construction, then apply generalized estimating equations to model air pollution hazard as a function of demographics.  This study focuses on residents of Asian background, and estimates effects for 7 different Asian ancestry groups.  They find that "an increase in proportion of Asian residents in census tracts is associated with significantly greater cancer risk from HAPs" [@GrineskiAsianAmericansdisproportionate2017 71].  

@BakhtsiyaravaEnvironmentalinequalitypollution2017 construct a spatial regression of immigrant population against toxic chemical releases, identified using the Toxics Release Inventory.  In this case, spatial regression is used as a control, to account for spatial autocorrelation in both migration and chemical releases.  They find that "immigrants from Europe and Latin America show a pollution advantage and tend to migrate to less polluted areas" relative to non-immigrants (defined as individuals who lived in the US the previous year); at the same time "for Mexican, Chinese, and Asian immigrants we found no relationship between toxin release and their concentration in PUMAs, indicating exposure patterns comparable to their native-born [sic] counterparts" [@BakhtsiyaravaEnvironmentalinequalitypollution2017 65]. 

@SilvaSpatialModelingIdentify2018 examine the presence of hydraulic fracturing wells in Ohio Census block groups.  They compare spatial and non-spatial regression models, that is, models that do and do not account for correlation between neighboring block groups.  They observe that the non-spatial model is susceptible to spatial confounding in this case; and find that increased median income and population density are associated with an increased presence of wells.  



# Data #

## Study Area ##

The study area for this project is California's Central Valley.  This study area was chosen for a number of reasons.  First, California is a  major US agricultural producer, producing over 13% of US agricultural value *[Ag report 2]*.  And the Central Valley is the largest center of California's agricultural production; 7 of the 10 most agriculturally productive counties in California are in the Central Valley *[Ag report 5]*.  Consequently, the Central Valley is also a major user of pesticides, including chlorpyrifos.  Second, demographically, the Central Valley is home to substantial populations of both Hispanic and non-Hispanic White residents, which raises the possibility of inequitable distributions of pesticide exposure, i.e., distributive environmental injustice.  Third, California's Department of Pesticide Regulation (DPR) makes available public, detailed, geocoded data on pesticide use in the state.  Combined with public data from the US Census, this makes it straightforward to retrieve data for a pesticide-related dEJ study.  

The primary units of analysis for this study are Census tracts and places.  *[tracts are exhaustive, but end up w/ large uncertainty for rural communities.  places are designed to pick out locally-identifiable communities, and so can have smaller uncertainty but also miss a fraction of the population]*

The Central Valley can be defined in a number of different ways.  Since the units of analysis for this project are tracts and places designated by the US Census, a county-based definition was judged to be most appropriate.  Sacramento County was excluded because, unlike the rest of the region, most of its area is urban.  17 other counties were used to define the Central Valley; see table \ref{tab.counties} and figure \ref{fig.chlor_use}.  

|Sacramento Valley| San Joaquin Valley |
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
Table: California counties comprising the Central Valley for the purposes of this study.  Counties are listed roughly in north-south order.  Counties in the left column are in the Sacramento Valley (northern half of the Central Valley); counties in the right column are in the San Joaquin Valley (southern half).\label{tab.counties}

## Pesticide Use ##

This study combines two datasets.  First, California's Department of Pesticide Regulation (DPR) releases annual public datasets for agricultural pesticide use across the state, known as Pesticide Use Reports (PUR) [^DPR].  These datasets provide data on pesticide use at the township and section level, including pounds of active chemical used on particular days and times.  As of this writing, the most recent data release covers 2015.  These data exclusively report agriculture use; this limitation is acceptable for studying chlorpyrifos, which has been banned for residential use, in the agriculture-heavy Central Valley.  

[^DPR]: <http://www.cdpr.ca.gov/docs/pur/purmain.htm>

Full datasets for 2011-2015 were retrieved.  Combining data over several years can help smooth annual variation in use due to variation in insect pest prevalence or other factors.  These datasets were combined and filtered on chemical name for chlorpyrifos.  Data for all counties across the state were used, not just the Central Valley, to avoid edge effects.  For example, while Sacramento County was not included in the study area, there is some chlorpyrifos use in the southeastern part of Sacramento County, very close to areas of San Joaquin County that are included in the study area.  These Sacramento County uses are incorporated into the analysis.  

DPR conducts automated quality checks, and flags errors and changed values in the data.  Entries with "potential duplicate" and "outlier" error codes were examined.  For the first few years of the dataset (2011 through 2013), it appears that DPR flagged but did not remove potential duplicates.  In later years, it appears that DPR removes potential duplicates in the released data.  Potential duplicates from the earlier years were therefore removed. 

DPR automatically identifies potential outliers as entries satisfying at least one of two criteria.  The first criterion is 200 lbs (solids and sprays) or 1000 lbs (fumigants) or more of active ingredient used per acre.  The second criterion is "50 times the median rate for all uses with the same pesticide product, crop treated, unit treated, and record type" (from the documentation file included in the data release).  Inspecting the entries with "potential outlier" flags, it appears that the main datafiles include corrected/estimated/modified values, not the original potential outliers.  The corrected/estimated/modified values given in the data were used without further adjustment.  

Chlorpyrifos uses are associated with townships and sections in the data; annual and all-study-period totals of active ingredient at the centroid of each 1 mile-square ($1.6 km \times 1.6 km$) section were calculated.  Because these centroids do not match the actual use locations (i.e., farm fields), the centroid totals might be unreliable for the smallest CTD, 1km (see discussion and table \ref{tab.ctd} below).  However, this error should be negligible for the other CTDs.  

All together, 1,113,398 use records for chlorpyrifos were identified in the DPR datasets for 2011-15.  After aggregating by sections and years, there were 31,789 records, with annualized use values ranging from $10^-2$ to $10^4$ lbs of active ingredient.  


## Demographics ##

The second primary dataset comprises American Community Survey (ACS) five-year estimates, from 2011-2015, for Census tracts and places in the 17 Central Valley Counties.  

For each tract and place, estimates and margin of error (MOE) values were retrieved for four categories of demographic variables:  *race and ethnicity* (Hispanic, non-Hispanic White, non-Hispanic Black, Indigenous, and Asian residents), *foreign-born noncitizens*, *children under 5* (who may be especially sensitive to chlorpyrifos exposure due to small body weight and critical neurodevelopmental stages), and *poverty* (individuals with an Census-determined income-poverty ratio below 1).  Because PUR data come only from agricultural uses, agricultural employment was also retrieved as a potential control.  (See also table \ref{tab.iv} and section \ref{sec.iv_selection}.)  

Population densities and proportions (e.g., the fraction of all residents who are Hispanic) were calculated for the total population and each of these ACS variables for every tract and place, using Census-recommended methods to calculate MOEs for these derived variables [@USCensusBureauAmericanCommunitySurvey, 11ff].  9 tracts and 15 places with estimated total population or total employment of 0 were excluded; 391 places and 1,044 tracts were used in all further analysis.  

Because much of the study area is rural, tract size and population density varied over four orders of magnitude, from 0.8 to 5600 residents per $km^2$.  Places cover 87% of the population, including 83% of non-Hispanic White and 88% of Hispanic residents, with population densities between 1.4 and 4400 residents per $km^2$.  However, places are geographically separated from each other, covering only 8% of the area of the tracts; in constructing contiguity-based spatial weights, 62% of places had no neighbors.  

The Modifiable Areal Unit Problem (MAUP) has been used to criticize dEJ projects [@SteelEnvironmentalJusticeValues2012;  @TaylorToxiccommunitiesenvironmental2014, 41ff].  "Egocentric neighborhood" methods have been used to address the MAUP in segregation research [@ReardonMeasuresSpatialSegregation2004].  However, these methods assume that populations are distributed evenly within each discrete region (e.g., each Census tract).  This assumption is inappropriate for this project, which includes many spatially heterogeneous rural regions.  In addition, the large geographic scale of this project would require millions of egocentric neighborhoods, and so trillions of distance calculations between neighborhoods and chlorpyrifos uses.  The Census tract-level distance calculations already pressed the limits of the computing power available for this phase of the project.  More fine-grained regions (e.g., Census block groups or blocks) would have multiplied uncertainty in the ACS estimates, and also likely would have exceeded the available computing power.  

As a compromise, block population counts from the 2010 Census were used to calculate weighted centroids for each tract and place.  These centroids more accurately represent the "average location" of the population in each tract, without requiring more computing power in the distance calculation step.  

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

$\beta^d$ is referred to here as the *decay coefficient*.  Given a CTD, the decay coefficient can be used in equation \ref{eqn.decay} to calculate the *distance-weighted local use* of chlorpyrifos at a location (tract or place centroid) $i$:  
\begin{align}
    q_i &= \sum_u q_u \beta^{d_{iu}} \label{eqn.tot_pot_exp}\\
        &= \sum_u q_u 0.37^{d_{iu}/CTD}
\end{align}
where $q_u$ is the total use at section $u$ and $d_{iu}$ is the Euclidean distance between the centroid of $i$ and the centroid of $u$.  To slightly account for the fact that the residents of a location are not located at its centroid, the decay coefficient is set to 1 whenever the centroid associated with use $u$ is within $i$, regardless of $d_{iu}$.  

It is important to stress that equations \ref{eqn.decay} and \ref{eqn.tot_pot_exp} provide, at best, rough estimates of potential exposure.  They do not take into account prevailing or occurrent winds, or other chemical transport media such as water.  Read as exposure estimates, they represent the chemical moving uniformly away from the source along a one-dimensional path, not dispersing over a two- or three-dimensional space surrounding the source at various speeds in various directions at various times.  Diverse and variable processes of application, fixation, and chemical transformation are represented as simple exponential decay.  I therefore refer to these aggregate statistics as "local use" rather than "exposure."  By comparison, @LuoSpatiallydistributedpesticide2010 use a physical-chemical model and PUR data to produce much more sophisticated estimates of chlorpyrifos loading.  

@MackayFateenvironmentlongrange2014 estimate the CTD for chlorpyrifos to be 62 km.  For robustness, this paper considers five CTD values, ranging from 1 km to 90 km.  See table \ref{tab.ctd} and figure \ref{fig.ctd}.  

| CTD (km) | $\beta$ |
|:---------|:--------|
| 1        | 0.370   |
| 10       | 0.905   |
| 30       | 0.967   |
| 60       | 0.984   |
| 90       | 0.989   |

Table: Characteristic Travel Distance (CTD) values used in this study, and corresponding decay-rate values $\beta$. \label{tab.ctd}

![Impact of Characteristic Travel Distance (CTD) value on weighted local use values in Census tracts. Panels correspond to the different CTD values used in this study. Weighted local use is the log (base 10) of aggregate chlorpyrifos use around each tract, scaled using the decay coefficient, 2011-15. Color scales are only roughly consistent between panels, with the scale midpoint set at 5 ($10^5 = 100,000$ lbs). \label{fig.ctd}](07_ctd.png)

It may be objected that this approach effectively assumes that, when chlorpyrifos is used, all of the chemical is volatilized.  This assumption is plainly false, but fortunately the current study only requires a much weaker assumption.  Decay coefficients are multiplicative, so that if $q_u$ is the amount used at section $u$ and $\pi$ is the volatilization rate, then $q_i = \sum_{u} \pi q_u \beta^{d_{iu}} = \pi \sum q_u \beta^{d_{iu}}$ is weighted local use at location $i$.  The second equality assumes that the volatilization rate $\pi$ is constant across the dataset.  This entails that the local use values can all be scaled by $\pi$.  Then, since the log of weighted local use is the dependent variable in the regression models below, the scaling coefficient $\pi$ becomes a constant term:  $\log q_i = \log \hat q_i + \log \pi$, where $\hat q_i$ is the un-scaled estimate above.  Thus, while the value of $\pi$ matters for estimating the intercept of the regression models, it does not matter for estimating the coefficients on the covariates.  The assumption of a constant volatilization rate is still unlikely, of course; it likely varies with the weather at the time chlorpyrifos is applied, for example.  However, insofar as variation in the volatilization rate is uncorrelated with the section $u$ and the demographic covariates of interest — which both seem plausible — the argument above can be modified to conclude that this variation will be bundled into the error term of the regression models, and so will not bias the coefficient estimates of interest.  


# Methods #

The primary analysis of this study is a regression of potential chlorpyrifos exposure against Census demographic data.  Separate models are constructed for each of the five CTD values listed in table \ref{tab.ctd}, as well as for tracts and places.  These two methodological choices give $5 \times 2 = 10$ models.  The software language R was used to clean and analyze all data, with significance use of the `tidyverse`, `tidycensus`, `sf`, `spdep`, and `tmap` packages [@WickhamtidyverseEasilyInstall2017; @WalkertidycensusLoadUS2018; @PebesmasfSimpleFeatures2018; @BivandspdepSpatialDependence2018; @TennekestmapThematicMaps2018].  Complete cleaning and analysis code is available at <https://github.com/dhicks/spatial_project>.  

This study uses an effects estimation approach, rather than a hypothesis testing approach [@GelmanFailureNullHypothesis2017].  That is, there are no tests of any null hypothesis that the effect of a given independent variable on local chlorpyrifos use is exactly 0.  In line with the aphorism that "everything is related to everything else" (sometimes called Tobler's first law of geography), in the kind of observational/field research setting of this study it is a priori highly implausible that any such effect is exactly 0.  Hypothesis testing therefore merely tells us whether we have enough data to detect the non-zero effect; this conclusion is not informative.  In addition, "stepwise" regression model selection methods (including or excluding variables based on whether their coefficients are statistically significant) are known to introduce bias and interpretability problems for both hypothesis testing and estimation [@gungAlgorithmsautomaticmodel2012 and citations therein].  

Instead of testing hypotheses, the aim of this study is to estimate effects as accurately possible, including estimates of the uncertainty of those first-order estimates.  Specifically, the modeling approach discussed below is designed to avoid or mitigate bias due to spatial correlation and errors in the independent variable measures.  In addition, two major "researcher degrees of freedom" [@SimmonsFalsePositivePsychology2011] — namely, geographic unit of analysis (tracts vs. places) and choice of CTD value — are tracked by analyzing 10 combinations of geography and CTD value in parallel.  The data collection and analysis process is carried out in a series of publicly available computer scripts, including two scripts devoted to exploratory data analysis (EDA) as well as interstitial comments on analytical decisions.  This research process enables future researchers to reflect on the extent to which particular decisions made in the course of study might have influenced the findings [@LeonelliReThinkingReproducibilityCriterion2018]. 


## Independent Variable Selection ## 
\label{sec.iv_selection}

Non-spatial exploratory data analysis indicated that, for almost all (>80%) tracts and places, almost all residents (>80%) were either Hispanic or non-Hispanic White.  Thus, Hispanic and non-Hispanic White proportions are strongly negatively correlated ($r = -.9$), and so non-Hispanic White proportion was dropped from the independent variable list.  There are moderate to large correlations ($r = .4-.8$) between Hispanic, noncitizen, young children, and poverty proportions; but these were not large enough to introduce large variance inflation factors that would have made the estimates uninterpretable.  

Four control variables were also used.  Because PUR data come only from agricultural uses, we expect chlorpyrifos use to be negatively correlated with population density and positively correlated with agricultural employment.  Density was not substantially correlated with any variables of interest (except for a low-moderate negative correlation with White; $r = -.4$).  Because density is left-bounded at 0 and ranges over several orders of magnitude, log density was used.  Agricultural employment has a moderate-strong correlation with Hispanic ($r = .67$).  To account for residual spatial heterogeneity, discussed below, county-level population density and agricultural employment variables were also included as controls.  

All independent variables are given in table \ref{tab.iv}, and descriptive statistics for both places and tracts are given in table \ref{tab.desc_stats} 

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


|geography |variable           | mean|   sd|   min|  max| Moran's I|
|:---------|:------------------|----:|----:|-----:|----:|---------:|
|places    |ag. employment     | 0.17| 0.19|  0.00| 1.00|      0.54|
|places    |Asian              | 0.03| 0.05|  0.00| 0.42|      0.29|
|places    |Black              | 0.02| 0.05|  0.00| 0.61|      0.13|
|places    |children           | 0.07| 0.04|  0.00| 0.21|      0.18|
|places    |pop. density (log) | 2.40| 0.73|  0.14| 3.64|      0.52|
|places    |Hispanic           | 0.43| 0.32|  0.00| 1.00|      0.74|
|places    |Indigenous         | 0.02| 0.04|  0.00| 0.45|     -0.01|
|places    |poverty            | 0.23| 0.17|  0.00| 1.00|      0.35|
|tracts    |ag. employment     | 0.10| 0.13|  0.00| 0.65|      0.65|
|tracts    |Asian              | 0.08| 0.08|  0.00| 0.53|      0.62|
|tracts    |Black              | 0.05| 0.06|  0.00| 0.46|      0.68|
|tracts    |children           | 0.07| 0.03|  0.00| 0.20|      0.28|
|tracts    |pop. density (log) | 2.75| 0.87| -0.11| 3.75|      0.46|
|tracts    |Hispanic           | 0.42| 0.24|  0.02| 0.98|      0.79|
|tracts    |Indigenous         | 0.01| 0.02|  0.00| 0.19|      0.16|
|tracts    |poverty            | 0.22| 0.13|  0.00| 0.64|      0.54|

Table:  Descriptive statistics for independent variables used in this study. All race-ethnicity groups other than Hispanic are non-Hispanic.  All IVs other than population density are proportion of total population in the tract or place. \label{tab.desc_stats}


## Spatial Weights ##

As part of the exploratory data analysis step (i.e., prior to regression analysis), 8 spatial weight constructions were considered for both tracts and places:  contiguity, inverse distance weights (with an outer limit of 50 km and a decay of $\frac{1}{d}$), and $k$-nearest-neighbors (KNN) with $k$ ranging from 3 to 8.  All weights were row-normalized.  As noted above, 62% of places had no contiguity-based neighbors.  

To examine the impact of different spatial weights constructions, Moran's $I$ was calculated for population densities corresponding to the independent variables identified above, e.g., density of Hispanic population, calculated as the number of Hispanic residents per square kilometer.  Moran's $I$ is a measure of spatial correlation, and is interpreted as a correlation coefficient:  values range from -1 to 1, with 1 indicating perfect spatial correlation (values at locations are identical to values at neighbors) and 0 indicating no spatial correlation.  

Among tracts, contiguity produced the highest values of Moran's $I$, followed by KNN, and finally distance weights.  For example, for Hispanic population density, Moran's $I$ was slightly less than .7 for contiguity weights; was between .45-.55 for KNN; and was slightly greater than .3 for distance weights.  

Among places, KNN and distance-based weights were similar, especially for larger values of K, while contiguity-based weights were much less for most variables.  For example for Hispanic population density, Moran's I was between .55-.65 for KNN; about .53 for distance weights; but only .3 for contiguity weights. 

KNN weights were selected for use in further analysis.  They produced reasonably high and consistent values of Moran's $I$ for both tracts and places.  

Initially, the bootstrap procedure described below was written to allow exploring several values of $k$ in the regression modeling stage.  However, runtime issues required reducing the number of free parameters explored in this stage.  Again, EDA indicated that values of $k$ from 3 to 8 for spatial weights produced roughly the same spatial correlation values for the independent variables, while CTD values produced multiple-order-of-magnitude differences across dependent variable values.  I therefore judged that different CTD values were likely to be a more important source of variation in effects than different spatial weights.  (See also @LeSageBiggestMythSpatial2014.)  $k=3$ was used as the sole spatial weights construction for all regression models.  


## Regression Specification ##

Exploratory data analysis indicated that local chlorpyrifos use values are non-Gaussian and range over multiple orders of magnitude across the entire study area.  In addition, for lower CTD values (1, 10), the left bound at 0 creates highly asymmetrical distributions.  A log transformation of the DV appeared to mitigate the left bound for the lower CTD values without qualitatively changing the distributions for the higher CTD values.  I therefore use a log-transformed dependent variable.  Base 10 was used for interpretability, so that dependent variable units are locally-weighted pound orders of magnitude.  

For higher CTD values (30, 60, 90), distributions are bi- or trimodal.  Plotting separate distributions for each county suggested that this was due to very different county-level baselines.  For counties with many tracts or places, the distribution of values appeared to be sufficiently Gaussian for standard regression.  County-level dummy variables were considered, and were found to improve homoscedasticity in non-spatial models,  but also introduced multicollinearity in spatial lagged models.  However, inspection of residuals in the spatial Durbin models suggested that lagged dependent variables substantially improved homoscedasticity without county-level dummies.  Still, even the spatial Durbin models exhibit some heteroscedasticity.  Separate county-level regression models were therefore constructed.  However, because many counties have few geographies (e.g., only 11 places), these models could only be fit for counties with 30 or more places or 50 or more tracts, and many of the estimates had high uncertainty.  Full-data models are therefore used for primary analysis, with county-level models used as a robustness/heterogeneity check.  

Spatial exploratory data analysis of both independent and dependent variables suggested substantial degrees of spatial autocorrelation on both sides of the regression formula.  A sequence of three model specifications was considered:  "standard" linear regression, without any spatial component (referred to simply as "regression" below); spatial regression with lagged independent variables, or "spatial lag X"; and spatial Durbin regression, which incorporates lags for both dependent and independent variables [@LeSageIntroductionSpatialEconometrics2009].  If the regression model is specified as 
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

KNN spatial weights, with $k=3$, were used for all of these models.  This sequence of models was fitted for both places and tracts and for CTD values 10 through 90.  Residual plots were examined for indications of heteroscedasticity, and $R^2$, AIC, and Moran's $I$ of the residuals were compared within model-dataset combinations.  


## Measurement Error and Bootstrap ##

All these regression model specifications assume perfect measurements in the independent variables.  However, ACS estimates can have large margins of error, especially for subpopulations of difficult-to-survey rural tracts [@SpielmanReducingUncertaintyAmerican2015].  This kind of measurement error can induce *attenuation bias* or *regression dilution*, in which correlation estimates are shrunk towards 0 [@FrostCorrectingregressiondilution2000].  In the context of dEJ analysis, attenuation bias is potentially a serious problem, insofar as it leads to the underestimation of environmental disparities.  That is, measurement error can make distributive environmental injustices seem less serious than they actually are.  

A bootstrap approach was developed to explore the effect of measurement error on the impact estimates.  In a "basic" or unparameterized bootstrap, independent samples are taken from the observations in the dataset (with replacement), forming a "resampled dataset" of the same size as the dataset.  These samples approximate drawing a new sample from the original population.  By calculating model statistics on a set of resampled datasets (1000 resampled datasets is common), we can estimate features of the sampling distribution of the population statistic.  For example, if $\hat\beta^1, \hat\beta^2, \ldots, \hat\beta^{1000}$ are regression coefficient estimates calculated on 1,000 resampled datasets, $\mathrm{sd}_l(\hat\beta^l)$ estimates the standard error of the coefficient estimate, and quantiles of the $\hat\beta^l$ distribution estimate quantiles of the parameter's sampling distribution.  

In the present study, the uncertainty of concern is with the ACS estimates for the IV values.  Assuming normality for these estimates, resamples can be drawn from the Gaussian distribution centered at the ACS point estimate with standard deviation the ACS margin of error.  Because spatial observations are not independent, the simple sampling method described above was not applicable to the current data.  Instead, each resampled dataset comprised the actual observed response values (potential chlorpyrifos exposure) at each actual location, with resampled IV values.  

500 such resampled datasets were constructed for each combination of geography (places or tracts) and CTD value (total $500 \times 2 \times 5 = 5000$ resampled datasets).  For each resampled dataset, a spatial Durbin model was fit and the following statistics were calculated: IV impacts for all covariates; $\rho$, the coefficient on the lagged dependent variables; and Moran's $I$ on the model residuals.  Impacts were calculated using a Monte Carlo method, with 100 draws for each resampled dataset.  These Monte Carlo draws were then combined at the geography-CTD level, producing $100 \times 500 = 50,000$ draws for the impacts for each IV-geography-CTD combination.  *[Means?]* and quantiles for these resampled distributions are reported below.  No resampling was done for county-level regression models.  



# Results #

Except when noted, the analysis in this section focuses on the full, non-resampled dataset and its models.  

## Model Selection and Evaluation ##

Model selection considered a "standard" non-spatial linear regression, spatial lagged X, and spatial Durbin models across each of the $10 = 2 \times 5$ geography-CTD combinations.  This can be interpreted as taking the "general-to-specific" approach to spatial model selection [@ElhorstAppliedSpatialEconometrics2010].  

Spatial Durbin models consistently outperformed both regression and spatial lagged X models.  Specifically, AIC and Moran's $I$ on the residuals were lowest with spatial Durbin models in every case.  For the spatial Durbin models, on both tracts and places $R^2>.8$ for all CTD values, and nearly 1 for CTD $\geq 30$, compared to values around .4-.5 for regression and .45-.55 for spatial lagged X.  These extremely high $R^2$ values appeared to be due to the way higher CTD values smoothed the response variable over space; see figure \ref{fig.ctd}.  

Discounting these high $R^2$ values, AIC and Moran's I statistics both still indicate that spatial Durbin models provide the best fit, and these models were selected for further analysis.  

However, Moran's $I$ was still substantially greater than 0 for all of the spatial Durbin models, with values of approximately $.07$ for tracts and $.15$ for places.  This suggests that there may be spatial non-stationarity; that is, the effects of the independent variables may vary across different sub-regions in the study area.  

County-level regression models were used to examine this possibility further.  Moran's $I$ was consistently close to 0 ($\pm .05$) for some county-CTD combinations; but was still greater than $.05$ for tracts in several counties, especially with larger CTD values.  In contrast, Moran's $I$ was consistently substantially negative (less than $-.05$) for places in Fresno and Kern county.  Thus, even at the county level, there are indications of spatial non-stationarity.  Non-stationary models were not explored in the current study, but may be an important direction for future work.  

Residual plots had strong signs of heteroscedasticity for both regression and spatial lagged X.  In general, residual variances were larger at small fitted values.  In tracts, with spatial Durbin models, any lingering indication of heteroscedasticity was due to 3 (out of 1044) locations with small fitted values and large (negative) residuals.  These locations can plausibly be treated as outliers, and so there do not appear to be any issues of heteroscedasticity at the tract level.  Places showed more heteroscedasticity, with approximately 15 (out of 391) locations with small fitted values and large (negative and positive) residuals in the spatial Durbin models.  This heteroscedasticity primarily means that the coefficient estimates are more uncertain (e.g., should have larger standard errors) than are reported below.  

Examination of residual scatterplots of the full dataset suggested that county-level heterogeneity might be responsible for this heteroscedasticity.  Supporting this assessment, residual scatterplots for county-level models generally did not show signs of worrisome heteroscedasticity.  The exceptions here were Butte (51 tracts) and Stanislaus (30 places, 94 tracts) counties, which may have been susceptible to small $n$ problems.  



## IV Impacts ##

Unlike "standard" non-spatial linear regression models, which treat observations as formally and mathematically independent, spatial models treat observations — locations or geographic units — as connected.  This means that changes in an IV at one location can influence the DV at another location, corresponding to the term $WX\theta$.  Further, the spatial Durbin model's lagged dependent variable term, $WY\rho$, introduces the possibility for feedback loops:  a change IV $\Delta x_i$ at location $l$ induces a change $\Delta y_{l'} = w_{l'l}\theta \Delta x_i$ in $y$ at neighbor $l'$, which feeds back to location $l$ as $w_{ll'} \Delta y_{l'} \rho$.  (This and the next paragraph generally follow @LeSageIntroductionSpatialEconometrics2009, §2.7.)  

Spatial econometricians have introduced the notion of *impacts* for the interpretation of regression coefficients under spatial feedbacks.  In non-spatial linear regression models without interaction, the coefficient $\beta_i$ for IV $x_i$ is identical to the partial derivative $\partial y / \partial x_i$.  $\beta_i$ can therefore be interpreted directly as the marginal effect of $x_i$ on $y$ (bracketing concerns about causal inference, etc.).  But in the spatial Durbin model, the partial derivative 
$$ \frac{\partial y}{\partial x_i} = (I_n - W\rho)^{-1}(I_n \beta_i + W \theta_i = S_i(W) $$
depends not just on the coefficients $\beta_i$ and $\theta_i$, but also the autoregression coefficient $\rho$.  And the value of this partial derivative at location $l$ depends on its connection to other locations, as encoded in $W$.  The *total impacts* for IV $x_i$ are formally defined as the mean row sum of $S_i(W)$, which corresponds to averaging $S_i(W)$ across locations.  These impacts are efficiently calculated using the `impacts()` function in the R package `spdep`, which uses a Monte Carlo method based on the traces of the powers of $W$.  The Monte Carlo method produces a sequence of estimates for total impact, which can be analyzed as an estimate of the sampling distribution for the true total impact.  

Figures \ref{fig.impacts_1}, \ref{fig.impacts_10}, \ref{fig.impacts_369} show total impact Monte Carlo estimates for IVs for CTD values of 1, 10, and 30-60-90 respectively. Total impact estimates were made for the spatial Durbin models fit on the observed dataset as well as each bootstrap resample of the IV estimates.  Impact estimates for the resampled datasets were combined into a single estimated sampling distribution for each geography-CTD-IV combination.  

*[rewrite captions]*

![Total IV impacts, CTD = 1 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_1}](11_impacts_1.png)

![Total IV impacts, CTD = 10 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_10}](11_impacts_10.png)

![Total IV impacts, CTD = 30, 60, and 90 km.  Violin plots show the distribution of Monte Carlo impact draws, with 100 draws for each of 500 block bootstrap resamples.  Horizontal lines in violin plots show 5%, 50%, and 95% quantiles.  Point ranges summarize 500 Monte Carlo impact draws for spatial Durbin models fitted on the complete dataset.  Point at 50%, top and bottom ends of ranges at 5% and 95% quantiles.  \label{fig.impacts_369}](11_impacts_369.png)

Comparing the three plots, the estimates tends to be much more uncertain for smaller CTD values.  Indeed, for CTD = 1 km, several estimates are uncertain across 10s of orders of magnitude (evaluating uncertainty by the width of the 5%-95% quantile interval).  By contrast, for CTD values of 30 or greater, uncertainty is often less than about 2 orders of magnitude, and in some cases substantially less than 1 order of magnitude.  In short, whatever physical validity different CTD values have, higher CTD values lead to more precise effects estimates.  The largest uncertainties are for the smallest subpopulations in the study area, namely, Asian, Black, Indigenous, and children proportions.  

The figures show both the bootstrap resamples (triangles and dashed lines) and observed data (circles and solid lines).  There is general agreement between both resamples and observed data for both point estimates and uncertainty.  Even when the point estimates are somewhat different, the degree of uncertainty (width of the 5%-95% quantile interval) is often very similar.  This is surprising because the resample-derived estimates include additional variation from measurement error.  

Comparing the resamples and observed data estimates, there are substantial indications of regression inflation rather than regression attenuation.  That is, the estimates inferred from the observed data are often further from 0 than the resample-derived estimates, which account for measurement error.  Given both sets of estimates, this trend means that the observed-derived estimates may be preferable in contexts, such as EJ-sensitive policymaking, where underestimation of potential harmful effects is more serious than overestimation.  This may appear to be in tension with the observation above, that the phenomenon of regression attenuation means that resampled estimates may be preferred when it is important to avoid underestimation in dEJ research.  However, the overarching goal is to *avoid underestimation* in order to *identify and redress distributive environmental injustice*.  It is appropriate to prefer the methods that best promote those goals in a particular case, even if that method is not the most effective way to promote those goals on average or in the long run.  

![Total impact estimates from the spatial Durbin models for CTD = 60 km. Total impact estimates $\hat\zeta$ are transformed as $\zeta_{trans} = 10^{\zeta/10}$.  $\zeta_{trans}$ can be interpreted as the relative change in local use when the corresponding IV increases by 10%.  For example, if $\zeta_{trans} = 1.5$, then local use is 50% greater when the corresponding IV is 10% greater.  Transformed values $>1$ therefore correspond to increases; values $<1$ correspond to decreases.  Line ranges give 90% confidence intervals. All estimates are based on the complete original datasets, not the block resamples.  \label{fig.impacts_backtrans}](11_impacts_backtrans.png)

|geography |CTD |IV                        | est. (transformed)| estimate| 95% CI|      |
|:---------|:---|:-------------------------|------------------:|--------:|------:|-----:|
|places    |60  |ag. employment            |               0.98|    -0.11|  -0.68|  0.40|
|tracts    |60  |ag. employment            |               1.22|     0.88|   0.51|  1.25|
|places    |60  |Asian                     |               1.06|     0.24|  -1.43|  1.69|
|tracts    |60  |Asian                     |               1.15|     0.60|   0.24|  0.97|
|places    |60  |Black                     |               0.63|    -2.00|  -3.51| -0.44|
|tracts    |60  |Black                     |               0.89|    -0.49|  -0.92| -0.07|
|places    |60  |children                  |               1.11|     0.45|  -1.73|  2.46|
|tracts    |60  |children                  |               0.91|    -0.41|  -1.64|  1.08|
|places    |60  |county ag. employment     |               1.42|     1.52|   0.62|  2.41|
|tracts    |60  |county ag. employment     |               1.38|     1.39|   0.85|  1.97|
|places    |60  |county pop. density (log) |               1.06|     0.24|   0.04|  0.45|
|tracts    |60  |county pop. density (log) |               1.02|     0.09|  -0.01|  0.19|
|places    |60  |Hispanic                  |               1.08|     0.35|  -0.16|  0.80|
|tracts    |60  |Hispanic                  |               1.19|     0.75|   0.50|  1.00|
|places    |60  |Indigenous                |               0.91|    -0.40|  -2.56|  1.80|
|tracts    |60  |Indigenous                |               0.64|    -1.90|  -4.74|  0.27|
|places    |60  |noncitizens               |               0.87|    -0.58|  -2.02|  0.82|
|tracts    |60  |noncitizens               |               0.64|    -1.94|  -2.60| -1.15|
|places    |60  |pop. density (log)        |               1.04|     0.17|   0.03|  0.30|
|tracts    |60  |pop. density (log)        |               1.06|     0.24|   0.20|  0.29|
|places    |60  |poverty                   |               1.12|     0.50|  -0.18|  1.08|
|tracts    |60  |poverty                   |               1.00|    -0.02|  -0.25|  0.20|
|places    |60  |women                     |               0.83|    -0.81|  -2.04|  0.28|
|tracts    |60  |women                     |               0.97|    -0.12|  -1.13|  0.81|


Table: Estimates of total impact (direct + indirect) from spatial Durbin models for CTD = 60km.  *IV*:  Independent variable.  *Est. (transformed)*: Point estimate for total impact, transformed as $10^{\zeta/10}$ to aid interpretation. \label{tab.impacts}

*[county-level plot]*

## Interpretation of Selected IVs ##

I focus interpretation on the total impact estimates for Hispanic, poverty, and agricultural employment proportion, and population density, for a CTD of 60 km, for the full data (not county-level) and without resampling.  These estimates are reported in figure \ref{fig.impacts_backtrans} and table \ref{tab.impacts}.  Transformed estimates are reported to aid interpretation:  when the estimated value of $\zeta$ is transformed as $\zeta_{trans} = 10^{\zeta/10}$, $\zeta_{trans}$ can be interpreted as the multiplicative change in local chlorpyrifos use when the corresponding IV increases by 10%.  For example, if $\zeta_{trans} = 1.5$ for some proportion variable, then a 10-percentage-point increase in this proportion is associated with a 50% increase in local chlorpyrifos use.  

A CTD of 60 km was chosen for this discussion for three reasons.  First, it yields good precision in the effects estimates — comparable to 90, better than 30, and much better than 10 or 1.  It also yields more skeptical or epistemically conservative estimates, in the sense that its estimates are closer to 0; although arguably this kind of skepticism is inappropriate in the context of environmental health [@HicksInductiveRiskRegulatory2018].  Finally, a CTD of 60 km is supported by the literature; namely, @MackayFateenvironmentlongrange2014 estimate the CTD for chlorpyrifos to be 62 km.  

The effects of agricultural employment and poverty proportions are ambiguous.  For agricultural employment, there is evidence is a positive effect of 1.22 (22% increase in chlorpyrifos use with each 10% increase in agricultural employment) for tracts; but a negligible effect of .98 for places.  For poverty, there is evidence of a negligible effect for tracts, and a positive effect (with greater uncertainty) for places.  

These inconsistencies between tracts and places may be due to the way places are constructed, namely, as a way to capture relatively dense population clusters in rural areas.  This process might exclude many agricultural workers and the rural poor, and that this exclusion might lead to biased effect estimates; though table \ref{tab.desc_stats} indicates that mean agricultural employment and poverty proportion are greater in places than in tracts.  County-level heterogeneity may also be a factor.  For example, poverty appears to have a positive effect for places in Stanislaus and Tulare counties, a negative effect in Butte county, and perhaps a negative effect in Kern county.  

There is consistent evidence of a positive effect for population density, with a (non-transformed) point estimate of .17 for places (95% CI .03-.30) and .24 for tracts (.20-.29).  For this covariate, an effect of .2 means that an order-of-magnitude increase in population density is associated with a .2 order-of-magnitude, or 58%, increase in chlorpyrifos use.  This is counterintuitive:  since chlorpyrifos is used primarily in agricultural areas, with low population density, we would expect to see a negative association.  

Local knowledge of the Central Valley suggests a potential explanation.  As indicated by figure \ref{fig.ctd}, with a CTD of 60, chlorpyrifos use is consistently much greater in the San Joaquin Valley (the southern part of the Central Valley).  The San Joaquin Valley has several small- and medium-sized cities, including Fresno (population approximately 500,000), Bakersfield (400,000), Stockton (300,000), and Modesto (200,000).  By contrast, after excluding Sacramento, the largest cities in the Sacramento Valley (the northern part of the Central Valley) are Redding and Chico (both approximately 90,000).  Population density may therefore be confounded, at least in part, with large-scale county or regional differences in chlorpyrifos use due to the kinds of crops grown or the prevalence of insect pests.  

However, the positive effect of population density remains apparent in most of the county-level regressions.  In Butte and Solano counties (Sacramento Valley) and Fresno, Kern, and Tulare counties (San Joaquin Valley), the regression models indicate that potential chlorpyrifos exposure increases with population density.  This effect is ambiguous for Stanislaus county (San Joaquin Valley), where there appears to be a larger but relatively highly uncertain effect for places but no clear effect for tracts.  And population density is negatively associated with potential chlorpyrifos exposure in San Joaquin county (San Joaquin Valley).  

Finally, there is consistent evidence that an increase in Hispanic proportion is associated with an increase in potential chlorpyrifos exposure.  This evidence is stronger with tracts than places, due in part to greater uncertainty of the latter estimates.  But this association appears across all CTD values and in 6 out of 11 county-level models.  For a CTD of 60 km, the estimated effect of a 10-point increase in Hispanic proportion is a 1.9-fold increase in potential chlorpyrifos exposure (95% CI 1.12-1.26).  A 60-point corresponding to the difference between a Hispanic-minority and Hispanic-majority tract, would be associated with a 1.97-4.00-fold increase in potential chlorpyrifos exposure.  

## Limitations ##

There are 8 notable limitations to this study.  First, there are limitations in the DPR public use data.  This dataset includes agricultural use only, and does not include industrial, commercial, state, or residential use of pesticides.  This may not be an issue for chlorpyrifos, which is banned for residential use in the US and appears to be used primarily for agricultural purposes.  However, for other pesticides that are widely used across sectors or in sectors not covered by the DPR data, the DPR data are likely to have significant gaps.  In addition, the DPR dataset tracks active ingredients, not the complex mixtures of product formulations that may make pesticides more toxic than the active ingredients alone.  

Second, and as discussed above, the methods used here model "local use" as a proxy for *potential* chlorpyrifos exposure, not actual exposure.  Also as discussed above, the methods used here effectively assume a highly simplified fate-and-transport model, in which the chemical is transported in a straight line through the air in all directions.  

Third, the spatial distribution of population data roughly models residential potential exposure to chlorpyrifos.  It does not model occupational exposure, children at school, or "take-home" occupational exposure (e.g., an agricultural worker brings home contaminated clothes that are handled by her children).  The weighted centroid method assumes populations are spread evenly over block groups; more refined distributions could be constructed using land use data, "egocentric neighborhoods" or raster methods.  These last two would likely require significantly more computing power than was available for most of this project.  

Fourth, highly correlated independent variables, such as Hispanic and non-citizen proportions, make it difficult to interpret effects.  Close correlations mean that the risk of omitted variable bias is small; but factor analysis or other approaches to composite variable construction could be used to incorporate correlated covariates.  

Fifth, the block bootstrap attempted here has not been validated, and it is not clear whether it is accurately accounting for attenuation bias.  Due to these limitations, the block bootstrap results were not incorporated into the discussion above, meaning those findings might suffer from attenuation bias.  

Six, there are signs of heteroscedasticity in the spatial Durbin models, especially for places.  There were indications throughout the study that county-level effects would address the non-Gaussian patterns in the data; but including county-level dummies created strong colinearities in the spatial regression models.  Other data sources might be incorporated to account for background baseline chlorpyrifos use rates, such as nearby crop species cultivated.  Or spatial random effects models — which allow effects to vary across space — might be used.  

Seventh, this project aims at identifying phenomena in the relationships between demographics and potential chlorpyrifos exposure.  It does not aim to developing theory, that is, offering explanations of the processes that produce these phenomena.  

Eighth and finally, this project examines distributive environmental injustice.  It does not attempt to examine procedural environmental injustice or misrecognition or disrespect.  While the data and methods used here are highly appropriate for identifying distributive environmental injustice phenomena, other approaches are necessary for understanding the social, political, and cultural processes that produce these unjust distributions of environmental hazard.  



# Discussion #

*[short para recapping findings]*

Reflecting on these results together with the long-term neurotoxic effects of chlorpyrifos suggests the possibility that Hispanic communities in the Central Valley may be subject to a process of cumulative disadvantage *[cite].  *[Water study] note that ethnicized political processes in the Central Valley have led to the political marginalization of Hispanic residents:  municipal boundaries in the Central Valley have often been drawn in ways that, sometimes intentionally *[?], exclude Hispanic communities.  

Consider the following scenario.  (1) Ethnic residential segregation leads to differential chlorpyrifos exposure, whether "by chance" or as a result of the political marginalization of Hispanic communities.  (2) Chlorpyrifos exposure, in turn, impairs the cognitive abilities of the children of these communities, compared to their peers in non-Hispanic White communities.  (3) These cognitive impairments reduce merit-based educational opportunities for Hispanic children as they grow up, which in turn (4) reduces the overall cognitive resources (e.g., number of bachelor's degree holders) and social capital of Hispanic communities.  Finally, closing the loop, (5) reduced cognitive resources and social capital may make these communities less politically effective, and specifically less able to effectively lobby for restrictions on chlorpyrifos use — at the county, state, or federal level — thereby increasing ethnic differences in chlorpyrifos exposure.  

In this scenario, even subtle differences in cognitive abilities may be substantially magnified by merit-based educational opportunities at step 3.  For example, as of this writing the University of California System requires California-resident prospective students to complete 15 college-preparatory courses with a grade point average of 3.0 or better and no greater lower than a C.  A small negative impact on cognitive ability might move a disproportionate share of Hispanic students from just above to just below these thresholds, denying them the opportunity to attend any UC campus.  

Each step of this scenario is likely to interact with other social, psychological, and physiological processes, from other chemical exposures to implicit bias.  However, it is highly plausible that these interactions will typically exacerbate rather than mitigate cumulative disadvantage.  (Interaction terms between Hispanic, noncitizen, and poverty proportions were considered for the regression models constructed in this study.  However, the interaction terms introduced extremely large variance inflation factors — in some cases multiple orders of magnitude — which made model interpretation impossible.)  

Empirically validating this complete scenario is beyond the scope of the present study.  Its plausibility does suggest important limitations, from a distributive environmental justice perspective, with two conventional aspects of human health risk assessment.  First, the scenario indicates that a plausible long-term adverse outcome of chlorpyrifos, at the population level, is increased exposure to chlorpyrifos.  This confounds the conventional distinction between hazard and exposure assessment.  

Second, multiple steps of this scenario involve social and political processes, and step 5 involves county, state, and federal political processes that cannot be reduced to individual-level processes.  The scenario also plays out over a multidecadal timescale.  These features of the scenario means that it cannot be assessed using methods designed for biophysiological processes operating within individual organisms.  These methods are appropriate for certain steps of the scenario; as are methods developed in epidemiology and public health.  But assessing steps 3-5 requires social science methods that are not conventionally used in any stage of human health risk assessment.  


# References #
