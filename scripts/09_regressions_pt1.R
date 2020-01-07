## Vanilla regression, lagged X, and spatial Durbin models
## The focus of this script is on model selection issues: R^2, AIC, Moran's I, and residual plots

#+ setup, include = FALSE
library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)
library(tmap)
library(broom)

library(assertthat)
library(car)

## A very crude augment() function for spatial Durbin models
## Returns a selected covariate and observed and predicted DV values
augment.sarlm = function(model, covar) {
    covar = enquo(covar)
    
    tibble(!!covar := model$X[,quo_name(covar)], 
           y_obs = model$y, 
           y_hat = as.numeric(predict(model)))
}

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sfl = read_rds(str_c(data_dir, '07_places_sfl.Rds'))
tracts_sfl = read_rds(str_c(data_dir, '07_tracts_sfl.Rds'))

## Common regression formula
reg_form = formula(log_w_use ~ hispanicP + noncitizensP +
                       blackP + indigenousP + asianP + womenP + 
                       childrenP + poverty_combP + 
                       ag_employedP + density_log10 + ag_employed_cP + density_log10_c)
covars = reg_form %>%
    as.character() %>%
    .[3] %>%
    str_split(fixed(' + ')) %>%
    unlist()

## Threshold for acceptable VIF
vif_threshold = 10

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

## County stats table ----
#+ tables
county_var_names = c('log.w.use' = 'weighted local use (log)', 
                     'ag.employed.cP' = 'ag. employment', 
                     'density.log10.c' = 'pop. density (log)')

county_stats = function (sfl, geography) {
    sfl %>%
        .[[1]] %>%
        as.data.frame() %>%
        select(county, ag_employed_cP, density_log10_c) %>%
        filter(!duplicated(.)) %>%
        rename(`ag. employment` = ag_employed_cP, 
               `pop. density (log)` = density_log10_c) %>%
        mutate(geography = geography)
}
county_stats_table = bind_rows(county_stats(tracts_sfl, 'tracts'), 
                               county_stats(places_sfl, 'places')) %>%
    select(geography, everything()) %>%
    arrange(geography, county)
knitr::kable(county_stats_table, format = 'markdown', 
             digits = 2)
knitr::kable(county_stats_table, format = 'markdown', 
             digits = 2) %>%
    write_lines('09_county_stats_table.txt')

logsumexp = function (x) log10(sum(10^x))
county_use = function(sfl, geography) {
    sfl %>%
        bind_rows() %>%
        as.data.frame() %>%
        select(ctd, county, log_w_use) %>%
        group_by(ctd, county) %>%
        summarize(mean = mean(log_w_use), 
                  sd = sd(log_w_use),
                  min = min(log_w_use), 
                  max = max(log_w_use), 
                  total = logsumexp(log_w_use)) %>%
        mutate(geography = geography) %>%
        ungroup()
}
county_use_table = bind_rows(county_use(tracts_sfl, 'tracts'),
                             county_use(places_sfl, 'places')) %>%
    select(ctd, geography, everything()) %>%
    mutate(ctd = as.integer(str_extract(ctd, '[0-9]+'))) %>%
    rename(CTD = ctd) %>%
    arrange(geography, county, CTD)
# knitr::kable(county_use_table, format = 'markdown', 
#              digits = 1)
# knitr::kable(county_use_table, format = 'markdown', 
#              digits = 1) %>%
#     write_lines('09_county_use_table.txt')


## IV table ----
## Nicely-printing variable names
var_names = c('ag.employedP' = 'ag. employment', 
              'asianP' = 'Asian', 
              'blackP' = 'Black', 
              'childrenP' = 'children', 
              'hispanicP' = 'Hispanic', 
              'noncitizensP' = 'noncitizens',
              'indigenousP' = 'Indigenous', 
              'density.log10' = 'pop. density (log)', 
              'poverty.combP' = 'poverty', 
              'womenP' = 'women')

moran = function(vec, weights) {
    moran.test(vec, weights) %>%
        .$estimate %>%
        .['Moran I statistic'] %>%
        `names<-`(NULL)
}
desc_stats = function(sfl, geography, weights) {
    sfl %>%
        .[[1]] %>%
        as.data.frame() %>% 
        select(one_of(covars), -contains('_c')) %>%
        rename_all(~ str_replace_all(., '_', '.')) %>%
        # rename(poverty.combP = poverty_combP,
        #        ag.employedP = ag_employedP, 
        #        density.log10 = density_log10) %>%
        summarize_all(funs(mean, sd,
                           min, max,
                           moran = moran(., weights))) %>%
        gather(key = var_stat, value) %>%
        separate(var_stat, into = c('var', 'stat'), sep = '_') %>%
        spread(key = stat, value) %>%
        mutate(geography = geography) %>%
        mutate(var = var_names[var]) %>%
        select(geography, variable = var, 
               mean, sd, min, max, moran)
}
desc_stats_table = bind_rows(
    desc_stats(places_sfl, 'places', weights_pl),
    desc_stats(tracts_sfl, 'tracts', weights_tr))
knitr::kable(desc_stats_table, format = 'markdown', 
             digits = 2)
knitr::kable(desc_stats_table, format = 'markdown', 
             digits = 2) %>%
    write_lines('09_desc_stats_table.txt')

dv_stats = function(sfl, geography, weights) {
    sfl %>% 
        as.data.frame() %>% 
        summarize_at(vars(log_w_use), 
                     funs(mean, sd,
                          min, max,
                          moran = moran(., weights))) %>% 
        mutate(geography = geography)
}

dv_stats_table = bind_rows(
    map_dfr(places_sfl, dv_stats, 
            'places', weights_pl, .id = 'ctd'),
    map_dfr(tracts_sfl, dv_stats, 
            'tracts', weights_tr, .id = 'ctd')
) %>% 
    separate(ctd, into = c('trash', 'CTD'), 
             sep = '_', convert = TRUE) %>% 
    select(CTD, geography, mean:moran)
knitr::kable(dv_stats_table, format = 'markdown', 
             digits = 2)
knitr::kable(dv_stats_table, format = 'markdown', 
             digits = 2) %>% 
    write_lines('09_dv_stats_table.txt')


## Places ----
#+ places
## Vanilla (non-spatial) linear regression
## Modest R^2, decreasing as CTD
## Moderate-high Moran I, increasing as CTD
models_lm_pl = map(places_sfl, 
                   ~ lm(reg_form, data = .))

## VIF below threshold?  
assert_that(!any(vif(models_lm_pl[[3]]) >= vif_threshold),
            msg = 'VIF at or above threshold')
## Any observations dropped? 
are_equal(nrow(models_lm_pl$ctd_30$model), nrow(places_sfl[[3]]), 
          msg = 'Observations dropped in model')

stats_lm_pl = map_dfr(models_lm_pl, glance, .id = 'ctd') %>%
    mutate(specification = 'lm', 
           model = models_lm_pl, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, moran, weights_pl))
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
           moran = map_dbl(residuals, moran, weights_pl))

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
    map_dfr(~ tibble(specification = 'sd', 
                     model = list(.)), 
            .id = 'ctd') %>%
    mutate(y = map(model, ~ .$y), 
           fitted = map(model, fitted),
           r.squared = map2_dbl(y, fitted, cor)^2, 
           AIC = map_dbl(model, AIC),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, moran, weights_pl))

## Combine and plot the model test statistics
## Across all CTDs, spatial Durbin models consistently have lowest AIC and Moran's I, and highest R^2
bind_rows(stats_lm_pl, 
          stats_slx_pl, 
          stats_sd_pl) %>%
    mutate(specification = fct_inorder(specification)) %>%
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

## Observed and predicted response values vs. hispanicP
augment_sd_pl = map_dfr(models_sd_pl, augment, hispanicP, .id = 'ctd')

ggplot(augment_sd_pl, aes(hispanicP)) +
    geom_point(aes(y = y_obs, color = 'observed')) +
    geom_point(aes(y = y_hat, color = 'predicted'), 
               alpha = .5) +
    geom_linerange(aes(ymin = y_obs, ymax = y_hat)) +
    scale_color_brewer(palette = 'Set1') +
    facet_wrap(vars(ctd), scales = 'free')


## Tracts ----
## Vanilla (non-spatial) linear regression
models_lm_tr = map(tracts_sfl, 
                   ~ lm(reg_form, data = .))

## VIF below threshold?  
assert_that(!any(vif(models_lm_tr[[3]]) >= vif_threshold),
            msg = 'VIF at or above threshold')
## Any observations dropped? 
are_equal(nrow(models_lm_tr$ctd_30$model), nrow(tracts_sfl[[3]]), 
          msg = 'Observations dropped in model')

stats_lm_tr = map_dfr(models_lm_tr, glance, .id = 'ctd') %>%
    mutate(specification = 'lm', 
           model = models_lm_tr, 
           fitted = map(model, fitted),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, moran, weights_tr))
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
           moran = map_dbl(residuals, moran, weights_tr))

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
    map_dfr(~ tibble(specification = 'sd', 
                     model = list(.)), 
            .id = 'ctd') %>%
    mutate(y = map(model, ~ .$y), 
           fitted = map(model, fitted),
           r.squared = map2_dbl(y, fitted, cor)^2, 
           AIC = map_dbl(model, AIC),
           residuals = map(model, residuals), 
           moran = map_dbl(residuals, moran, weights_tr))

## Combine and plot the model test statistics
## Across all CTDs, spatial Durbin models consistently have lowest AIC and Moran's I, and highest R^2
bind_rows(stats_lm_tr, 
          stats_slx_tr, 
          stats_sd_tr) %>%
    mutate(specification = fct_inorder(specification)) %>%
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

## Predicted and observed y vs. hispanicP
augment_sd_tr = map_dfr(models_sd_tr, augment, hispanicP, .id = 'ctd')

ggplot(augment_sd_tr, aes(hispanicP)) +
    geom_point(aes(y = y_obs, color = 'observed')) +
    geom_point(aes(y = y_hat, color = 'predicted'), 
               alpha = .5) +
    geom_linerange(aes(ymin = y_obs, ymax = y_hat)) +
    scale_color_brewer(palette = 'Set1') +
    facet_wrap(vars(ctd), scales = 'free')


## Peeking at the impacts
summary(models_sd_tr[[3]])
impacts(models_sd_tr[[3]], tr = traces_tr, R = 500) %>%
    .$sres %>%
    .$total %>%
    as_tibble() %>%
    gather(key = term, value = estimate) %>%
    group_by(term) %>%
    summarize(mean = mean(estimate), se = sd(estimate), 
              conf_low = quantile(estimate, probs = .025), conf_high = quantile(estimate, probs = .975))
