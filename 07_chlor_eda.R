library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)
library(tmap)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

## County borders -----
## https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
counties_sf = read_sf(str_c(data_dir, 'cb_2015_us_county_20m'), 
                      'cb_2015_us_county_20m') %>%
    filter(STATEFP == '06') %>%
    st_transform('+proj=utm +zone=10 +datum=NAD83')

## Location data ----
chlor_sf = read_rds(str_c(data_dir, '01_chlor_sf.Rds'))

places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number()) %>%
    ## Remove places w/ total population or total employed 0
    filter(total_popE != 0, total_employedE != 0) %>%
    ## Join w/ county names
    mutate(county_geoid = str_extract(GEOID, '[0-9]{5}')) %>%
    left_join(select(as.data.frame(counties_sf), 
                     county_geoid = GEOID, county = NAME)) %>%
    ## Population proportions
    mutate(whiteP = whiteE / total_popE, 
           whitePM = moe_prop(whiteE, total_popE, whiteM, total_popM),
           blackP = blackE / total_popE, 
           blackPM = moe_prop(blackE, total_popE, blackM, total_popM),
           indigenousP = indigenousE / total_popE, 
           indigenousPM = moe_prop(indigenousE, total_popE, indigenousM, total_popM),
           asianP = asianE / total_popE, 
           asianPM = moe_prop(asianE, total_popE, asianM, total_popM),
           hispanicP = hispanicE / total_popE, 
           hispanicPM = moe_prop(hispanicE, total_popE, hispanicM, total_popM),
           noncitizensP = noncitizensE / total_popE,
           noncitizensPM = moe_prop(noncitizensE, total_popE, 
                                    noncitizensM, total_popM),
           childrenP = childrenE / total_popE, 
           childrenPM = moe_prop(childrenE, total_popE, childrenM, total_popM),
           
           poverty_combE = povertyE + extreme_povertyE,
           poverty_combM = sqrt(povertyM^2 + extreme_povertyM^2), 
           poverty_combP = poverty_combE / total_popE, 
           poverty_combPM = moe_prop(poverty_combE, total_popE, 
                                     poverty_combM, total_popM), 
           hisp_povertyP = hisp_povertyE / hispanicE, 
           hisp_povertyPM = moe_prop(hisp_povertyE, hispanicE, hisp_povertyM, hispanicM), 
           ag_employedP = ag_employedE / total_employedE,
           ag_employedPM = moe_prop(ag_employedE, total_employedE, ag_employedM, total_employedM)
    ) %>%
    ## Population densities
    mutate_at(vars(whiteE, blackE, indigenousE, asianE, hispanicE, 
                   noncitizensE, childrenE, poverty_combE, 
                   hisp_povertyP, ag_employedP), 
              funs(D = . / units::drop_units(area))) %>%
    ## Bind w/ potential exposure
    left_join(read_rds(str_c(data_dir, '06_w_use_places.Rds')))
    

## EDA ----
## Non-Gaussian distributions
ggplot(tracts_sf, aes(w_use, color = ctd)) + 
    geom_density() +
    facet_wrap(~ ctd, scales = 'free')

## Breaking things out by counties helps a bit
tracts_sf %>%
    filter(ctd == 'ctd_30') %>%
ggplot(aes(w_use, color = ctd)) + 
    geom_density() +
    geom_rug() +
    facet_wrap(~ county, scales = 'free')
## What about a log?  
tracts_sf %>%
    filter(ctd %in% c('ctd_1')) %>%
    ggplot(aes(log10(w_use), color = ctd)) + 
    geom_density() + 
    geom_rug() +
    facet_wrap(~ county, scales = 'free')

## Pop. density
ggplot(tracts_sf, aes(log10(densityE), log10(w_use))) + 
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~ ctd, scales = 'free')

## How about ag employment?
ggplot(tracts_sf, aes(ag_employedP, w_use)) +
    geom_point() +
    geom_smooth(method = 'lm')


tm_shape(places_sf) + 
    tm_fill('ctd_60', style = 'cont', palette = '-RdBu') +
tm_shape(counties_sf) +
    tm_borders() +
tm_shape(chlor_sf) +
    tm_symbols(col = 'total_use', size = .1, alpha = .1, border.lwd = 0)

ggplot(places_sf, aes(hispanicP, w_use)) + 
    geom_point() +
    geom_smooth(method = 'lm')

library(lme4)
model = lmer(w_use ~ 0 + hispanicP + densityE + ag_employedP + (1|county), 
             data = places_sf)
summary(model)
## NB Heteroscedasticity, even w/ county-level random intercept
plot(model)

## Oh hey, spatial autocorrelation in the residuals
places_sf %>%
    mutate(.resid = residuals(model)) %>%
tm_shape() + 
    tm_fill('.resid', style = 'cont')
