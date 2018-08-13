## Vanilla regression, lagged X, and spatial Durbin models
## The focus of this script is on model selection issues: R^2, AIC, Moran's I, and residual plots

#+ setup, include = FALSE
library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)
library(tmap)
library(broom)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sfl = read_rds(str_c(data_dir, '07_places_sfl.Rds'))
tracts_sfl = read_rds(str_c(data_dir, '07_tracts_sfl.Rds'))

## Common regression formula
reg_form = formula(log_w_use ~ hispanicP + blackP + indigenousP + asianP + childrenP + poverty_combP + ag_employedP + density_log10)

## KNN
construct_traces = function(weights) {
    weights %>%
        as('CsparseMatrix') %>%
        trW()
}

weights_pl = read_rds(str_c(data_dir, '07_places_weights.Rds'))
weights_tr = read_rds(str_c(data_dir, '07_tracts_weights.Rds'))

traces_pl = construct_traces(weights_pl)
traces_tr = construct_traces(weights_tr)


## Places ----
#+ places
## Vanilla (non-spatial) linear regression
## Modest R^2, decreasing as CTD
## Moderate-high Moran I, increasing as CTD
models_lm_pl = map(places_sfl, 
                   ~ lm(reg_form, data = .))

car::vif(models_lm_pl[[3]])

stats_lm_pl = map_dfr(models_lm_pl, glance, .id = 'ctd') %>%
    mutate(specification = 'lm', 
           model = models_lm_pl, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_pl, 500) %>%
                               glance() %>%
                               pull(statistic)))
stats_lm_pl %>%
    select(specification, ctd, r.squared, AIC, moran) %>%
    knitr::kable()

## Lagged X spatial regression
## Higher R^2, lower AIC, higher Moran's I
models_slx_pl = map(places_sfl, 
                    ~ lmSLX(reg_form, data = ., weights_pl))
stats_slx_pl = map_dfr(models_slx_pl, glance, .id = 'ctd') %>%
    mutate(specification = 'slx',
           model = models_slx_pl, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_pl, 500) %>%
                               glance() %>%
                               pull(statistic)))

stats_slx_pl %>%
    select(specification, ctd, r.squared, AIC, moran) %>%
    knitr::kable()


## Spatial Durbin model
models_sd_pl = map(places_sfl, 
                   ~ lagsarlm(reg_form, data = ., 
                              weights_pl, trs = traces_pl, 
                              type = 'Durbin'))
## No glance() for sarlm
stats_sd_pl = models_sd_pl %>%
    map(summary, Nagelkerke = TRUE) %>% 
    map_dfr(~ tibble(r.squared = .$NK, 
                 # AIC = .$), where is this?
        ), .id = 'ctd') %>%
    mutate(specification = 'sd',
           model = models_sd_pl, 
           AIC = map_dbl(model, AIC),
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_pl, 500) %>%
                               glance() %>%
                               pull(statistic)))

## Combine and plot the model test statistics
## Across all CTDs, spatial Durbin models consistently have lowest AIC and Moran's I, and highest R^2
bind_rows(stats_lm_pl, 
          stats_slx_pl, 
          stats_sd_pl) %>%
    select(specification, ctd,
           r.squared, AIC, moran) %>%
    gather(key = statistic, value, r.squared:moran) %>%
    ggplot(aes(ctd, value, color = specification)) +
    geom_point() +
    facet_wrap(~ statistic, scales = 'free')

## Residual plots
## Spatial Durbin models consistently have less heteroscedasticity than other two; 
## though heteroscedasticity might still be an issue
bind_rows(stats_lm_pl, 
          stats_slx_pl, 
          stats_sd_pl) %>%
    mutate(specification = fct_inorder(specification)) %>%
    select(specification, ctd, .fit = fitted, .resid = residuals) %>%
    unnest() %>%
    ggplot(aes(.fit, .resid)) + 
    geom_point() +
    geom_smooth() +
    facet_wrap(specification ~ ctd, scales = 'free', nrow = 3)



## Tracts ----
## Vanilla (non-spatial) linear regression
models_lm_tr = map(tracts_sfl, 
                   ~ lm(reg_form, data = .))

car::vif(models_lm_tr[[3]])

stats_lm_tr = map_dfr(models_lm_tr, glance, .id = 'ctd') %>%
    mutate(specification = 'lm', 
           model = models_lm_tr, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_tr, 500) %>%
                               glance() %>%
                               pull(statistic)))
stats_lm_tr %>%
    select(specification, ctd, r.squared, AIC, moran) %>%
    knitr::kable()

## Lagged X spatial regression
## Higher R^2, lower AIC, higher Moran's I
models_slx_tr = map(tracts_sfl, 
                    ~ lmSLX(reg_form, data = ., weights_tr))
stats_slx_tr = map_dfr(models_slx_tr, glance, .id = 'ctd') %>%
    mutate(specification = 'slx',
           model = models_slx_tr, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_tr, 500) %>%
                               glance() %>%
                               pull(statistic)))

stats_slx_tr %>%
    select(specification, ctd, r.squared, AIC, moran) %>%
    knitr::kable()


## Spatial Durbin model
models_sd_tr = map(tracts_sfl, 
                   ~ lagsarlm(reg_form, data = ., 
                              weights_tr, trs = traces_tr, 
                              type = 'Durbin'))
## No glance() for sarlm
stats_sd_tr = models_sd_tr %>%
    map(summary, Nagelkerke = TRUE) %>% 
    map_dfr(~ tibble(r.squared = .$NK, 
                     # AIC = .$), where is this?
    ), .id = 'ctd') %>%
    mutate(specification = 'sd',
           model = models_sd_tr, 
           AIC = map_dbl(model, AIC),
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, 
                           ~ moran.mc(., weights_tr, 500) %>%
                               glance() %>%
                               pull(statistic)))

## Combine and plot the model test statistics
## Across all CTDs, spatial Durbin models consistently have lowest AIC and Moran's I, and highest R^2
bind_rows(stats_lm_tr, 
          stats_slx_tr, 
          stats_sd_tr) %>%
    select(specification, ctd,
           r.squared, AIC, moran) %>%
    gather(key = statistic, value, r.squared:moran) %>%
    ggplot(aes(ctd, value, color = specification)) +
    geom_point() +
    facet_wrap(~ statistic, scales = 'free')

## Residual plots
## Spatial Durbin models consistently have less heteroscedasticity than other two; 
## though heteroscedasticity might still be an issue
bind_rows(stats_lm_tr, 
          stats_slx_tr, 
          stats_sd_tr) %>%
    mutate(specification = fct_inorder(specification)) %>%
    select(specification, ctd, .fit = fitted, .resid = residuals) %>%
    unnest() %>%
    ggplot(aes(.fit, .resid)) + 
    geom_point() +
    geom_smooth() +
    facet_wrap(specification ~ ctd, scales = 'free', nrow = 3)
