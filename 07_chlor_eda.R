## This script
## 1. constructs the variables used in the regression models in the following scripts, then
## 2. conducts some further EDA for these regressions

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

prep_data = function (locations_file, use_file) {
    dataf = read_rds(locations_file) %>%
        mutate(idx = row_number()) %>%
        select(idx, everything()) %>%
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
        left_join(read_rds(use_file)) %>%
        ## Logged use
        mutate(log_w_use = log10(w_use)) %>%
        ## Split by CTD
        # plyr::dlply('ctd', identity)
        split(., .$ctd)
    return(dataf)
}

places_sfl = prep_data(str_c(data_dir, '02_places_sf.Rds'), 
                       str_c(data_dir, '06_w_use_places.Rds')) %>%
    map(rename, place_idx = idx)
tracts_sfl = prep_data(str_c(data_dir, '02_tracts_sf.Rds'), 
                       str_c(data_dir, '06_w_use_tracts.Rds')) %>%
    map(rename, tract_idx = idx)

write_rds(places_sfl, str_c(data_dir, '07_places_sfl.Rds'))
write_rds(tracts_sfl, str_c(data_dir, '07_tracts_sfl.Rds'))

## Spatial weights, for use in the regression scripts
construct_weights = function(sfl, k = 3) {
    sfl %>%
        .[[1]] %>%
        st_centroid %>%
        st_coordinates() %>%
        knearneigh(k = k) %>%
        knn2nb() %>%
        nb2listw(style = 'W')
}

k = 3
places_sfl %>%
    construct_weights(k = k) %>%
    write_rds(str_c(data_dir, '07_places_weights.Rds'))
tracts_sfl %>%
    construct_weights(k = k) %>%
    write_rds(str_c(data_dir, '07_tracts_weights.Rds'))


## EDA ----
## Non-Gaussian distributions
tracts_sfl %>%
    imap(~ ggplot(.x, aes(w_use)) +
            geom_density() +
            ggtitle(.y)) %>%
    cowplot::plot_grid(plotlist = .)

## Breaking things out by counties helps a bit
ggplot(tracts_sfl$ctd_30, aes(w_use)) + 
    geom_density() +
    geom_rug() +
    facet_wrap(~ county, scales = 'free')

## What about a log?  
tracts_sfl %>%
    imap(~ ggplot(.x, aes(w_use)) +
            geom_density() +
            scale_x_log10() +
            ggtitle(.y)) %>%
    cowplot::plot_grid(plotlist = .)

## Pop. density
imap(tracts_sfl, 
    ~ ggplot(.x, aes(log10(densityE), log10(w_use))) + 
        geom_point() +
        geom_smooth(method = 'lm') +
        ggtitle(.y)) %>%
    cowplot::plot_grid(plotlist = .)

## How about ag employment?
imap(tracts_sfl, 
    ~ ggplot(.x, aes(log10(ag_employedP), log10(w_use))) + 
        geom_point() +
        geom_smooth(method = 'lm') +
        ggtitle(.y)) %>%
    cowplot::plot_grid(plotlist = .)

## Map of use for CTD = 60
tm_shape(tracts_sfl$ctd_60) + 
    tm_fill('log_w_use', style = 'cont', palette = '-RdBu') +
    tm_shape(counties_sf) +
    tm_borders() +
    tm_shape(chlor_sf) +
    tm_symbols(col = 'total_use', size = .1, alpha = .05, border.lwd = 0)

# tracts_sfl %>%
#     map(as_tibble) %>%
#     bind_rows() %>% 
#     ggplot(aes(log_w_use)) +
#     geom_density() +
#     facet_wrap(~ ctd, scales = 'free')

## Faceted use map for all CTD
tracts_sfl %>%
    do.call(rbind, .) %>%
    mutate(ctd = case_when(ctd == 'ctd_1' ~ '1', 
                           ctd == 'ctd_10' ~ '10', 
                           ctd == 'ctd_30' ~ '30', 
                           ctd == 'ctd_60' ~ '60', 
                           ctd == 'ctd_90' ~ '90')) %>%
    tm_shape() +
    tm_fill('log_w_use', style = 'cont', 
            palette = '-RdBu', midpoint = 5, 
            title = 'Weighted\nlocal use') +
    tm_facets(by = 'ctd', free.scales = TRUE) +
    tm_shape(counties_sf) +
    tm_borders() +
    tm_scale_bar(position = c('left', 'bottom')) +
    tm_layout(legend.position = c('right', 'top'))
tmap_save(filename = '07_ctd.png', 
          width = 6, height = 7, units = 'in')
    
imap(places_sfl, 
     ~ ggplot(.x, aes(hispanicP, w_use)) + 
         geom_point() +
         geom_smooth(method = 'lm') +
         ggtitle(.y)) %>%
    cowplot::plot_grid(plotlist = .)

library(lme4)
model = lmer(w_use ~ 0 + hispanicP + densityE + ag_employedP + (1|county), 
             data = tracts_sfl$ctd_30)
summary(model)
## NB Heteroscedasticity, even w/ county-level random intercept
plot(model)

## Oh hey, spatial autocorrelation in the residuals
tracts_sfl$ctd_30 %>%
    mutate(.resid = residuals(model)) %>%
    tm_shape() + 
    tm_fill('.resid', style = 'cont')
