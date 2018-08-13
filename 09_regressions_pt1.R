## Vanilla regression and lagged X models

#+ setup, include = FALSE
library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)
library(tmap)
library(broom)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number()) %>%
    ## Remove places w/ total population or total employed 0
    filter(total_popE != 0, total_employedE != 0) %>%
    ## Log density
    mutate(log_densityE = log10(densityE)) %>%
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
    ## Bind w/ locally-weighted total use
    left_join(read_rds(str_c(data_dir, '06_w_use_places.Rds'))) %>%
    ## Logged use
    mutate(log_w_use = log10(w_use))

places10 = filter(places_sf, ctd == 'ctd_10')
places60 = filter(places_sf, ctd == 'ctd_60')

tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
    mutate(tracts_idx = row_number()) %>%
    ## Remove tracts w/ total population or total employed 0
    filter(total_popE != 0, total_employedE != 0) %>%
    ## Log density
    mutate(log_densityE = log10(densityE)) %>%
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
    ## Bind w/ locally-weighted total use
    left_join(read_rds(str_c(data_dir, '06_w_use_tracts.Rds'))) %>%
    ## Logged use
    mutate(log_w_use = log10(w_use))

tracts10 = filter(tracts_sf, ctd == 'ctd_10')
tracts60 = filter(tracts_sf, ctd == 'ctd_60')

## Common regression formula
reg_form = formula(log10(w_use) ~ hispanicP + blackP + indigenousP + asianP + childrenP + poverty_combP + ag_employedP + log10(densityE))

## KNN
k = 3
weights_pl = places_sf %>%
    filter(ctd == 'ctd_30') %>%
    st_centroid %>%
    st_coordinates() %>%
    knearneigh(k = k) %>%
    knn2nb() %>%
    nb2listw(style = 'W')
traces_pl = weights_pl %>%
    as('CsparseMatrix') %>%
    trW()
weights_tr = tracts_sf %>%
    filter(ctd == 'ctd_30') %>%
    st_centroid %>%
    st_coordinates() %>%
    knearneigh(k = k) %>%
    knn2nb() %>%
    nb2listw(style = 'W')
traces_tr = weights_tr %>%
    as('CsparseMatrix') %>%
    trW()


## Places, CTD = 10 ----
#+ places_10
## Vanilla (non-spatial) linear regression
model_lm_pl10 = lm(reg_form, data = places10)
summary(model_lm_pl10)
AIC(model_lm_pl10)
ggplot(augment(model_lm_pl10), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_lm_pl10), weights_pl, 500)

## Lagged X spatial regression
model_slx_pl10 = lmSLX(reg_form, data = places10, weights_pl)
summary(model_slx_pl10)
AIC(model_slx_pl10)
ggplot(augment(model_slx_pl10), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_slx_pl10), weights_pl, 500)

## Spatial Durbin model
model_durbin_pl10 = lagsarlm(reg_form, 
                             data = places10, 
                             weights_pl, trs = traces_pl, type = 'Durbin')
summary(model_durbin_pl10)
cor(places10$log_w_use, predict(model_durbin_pl10, pred.type = 'trend'))
AIC(model_durbin_pl10)
tibble(.fitted = as.numeric(predict(model_durbin_pl10, 
                                    pred.type = 'trend')), 
       .resid = residuals(model_durbin_pl10)) %>%
    ggplot(aes(.fitted, .resid)) + 
    geom_point()
moran.mc(residuals(model_durbin_pl10), weights_pl, 500)

## - AIC is much smaller for spatial Durbin model (-119 vs. 773 vs. 875)
## - Still some heteroscedasticity in spatial Durbin model, but much less than both vanilla and lagged X
## - Moran's I:  .75; .84; .15
## - R^2: .40; .56; .65

## Places, CTD = 60 ----
#+ places_60
## Vanilla (non-spatial) linear regression
model_lm_pl60 = lm(reg_form, data = places60)
AIC(model_lm_pl60)
summary(model_lm_pl60)
ggplot(augment(model_lm_pl60), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_lm_pl60), weights_pl, 500)

## Lagged X spatial regression
model_slx_pl60 = lmSLX(reg_form, data = places60, weights_pl)
summary(model_slx_pl60)
AIC(model_slx_pl60)
ggplot(augment(model_slx_pl60), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_slx_pl60), weights_pl, 500)

## Spatial Durbin model
model_durbin_pl60 = lagsarlm(reg_form, 
                             data = places60, 
                             weights_pl, trs = traces_pl, type = 'Durbin')
summary(model_durbin_pl60)
cor(places60$log_w_use, predict(model_durbin_pl60, pred.type = 'trend'))
AIC(model_durbin_pl60)
tibble(.fitted = as.numeric(predict(model_durbin_pl60, 
                                    pred.type = 'trend')), 
       .resid = residuals(model_durbin_pl60)) %>%
    ggplot(aes(.fitted, .resid)) + 
    geom_point()
moran.mc(residuals(model_durbin_pl60), weights_pl, 500)

## - AIC is much smaller for spatial Durbin model (-119 vs. 773 vs. 876)
## - Still some heteroscedasticity in spatial Durbin model, but much less than both vanilla and lagged X
## - Moran's I: .75, .84, .15
## - R^2: .38; .55; .62



## Tracts, CTD = 10 ----
#+ tracts_10
## Vanilla (non-spatial) linear regression
model_lm_tr10 = lm(reg_form, data = tracts10)
summary(model_lm_tr10)
AIC(model_lm_tr10)
ggplot(augment(model_lm_tr10), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_lm_tr10), weights_tr, 500)

## Lagged X spatial regression
model_slx_tr10 = lmSLX(reg_form, data = tracts10, weights_tr)
summary(model_slx_tr10)
AIC(model_slx_tr10)
ggplot(augment(model_slx_tr10), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_slx_tr10), weights_tr, 500)

## Spatial Durbin model
model_durbin_tr10 = lagsarlm(reg_form, 
                             data = tracts10, 
                             weights_tr, trs = traces_tr, type = 'Durbin')
summary(model_durbin_tr10)
cor(tracts10$log_w_use, predict(model_durbin_tr10, pred.type = 'trend'))
AIC(model_durbin_tr10)
tibble(.fitted = as.numeric(predict(model_durbin_tr10, 
                                    pred.type = 'trend')), 
       .resid = residuals(model_durbin_tr10)) %>%
    ggplot(aes(.fitted, .resid)) + 
    geom_point()
moran.mc(residuals(model_durbin_tr10), weights_tr, 500)

## - AIC is much smaller for spatial Durbin model (-998 vs. 1654 vs. 1816)
## - Still some heteroscedasticity in spatial Durbin model, but much less than both vanilla and lagged X
## - Moran's I: .82; .89; .08
## - R^2: .34; .44; .49


## tracts, CTD = 60 ----
#+ tracts_60
## Vanilla (non-spatial) linear regression
model_lm_tr60 = lm(reg_form, data = tracts60)
summary(model_lm_tr60)
AIC(model_lm_tr60)
ggplot(augment(model_lm_tr60), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_lm_tr60), weights_tr, 500)

## Lagged X spatial regression
model_slx_tr60 = lmSLX(reg_form, data = tracts60, weights_tr)
summary(model_slx_tr60)
AIC(model_slx_tr60)
ggplot(augment(model_slx_tr60), aes(.fitted, .resid)) + geom_point()
moran.mc(residuals(model_slx_tr60), weights_tr, 500)

## Spatial Durbin model
model_durbin_tr60 = lagsarlm(reg_form, 
                             data = tracts60, 
                             weights_tr, trs = traces_tr, type = 'Durbin')
summary(model_durbin_tr60)
cor(tracts60$log_w_use, predict(model_durbin_tr60, pred.type = 'trend'))
AIC(model_durbin_tr60)
tibble(.fitted = as.numeric(predict(model_durbin_tr60, 
                                    pred.type = 'trend')), 
       .resid = residuals(model_durbin_tr60)) %>%
    ggplot(aes(.fitted, .resid)) + 
    geom_point()
moran.mc(residuals(model_durbin_tr60), weights_tr, 500)

## - AIC is much smaller for spatial Durbin model (-5039 vs. -34 vs. 99)
## - Still maybe some heteroscedasticity in spatial Durbin model, but maybe mostly due to outliers
## - Moran's I: .88; .95; .07
## - R^2: .37; .45; .48

