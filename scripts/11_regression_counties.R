library(tidyverse)
library(sf)
library(spdep)
library(doSNOW)
library(tictoc)

## Load data ----
# data_dir = '../data/'
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sfl = read_rds(str_c(data_dir, '07_places_sfl.Rds'))
tracts_sfl = read_rds(str_c(data_dir, '07_tracts_sfl.Rds'))

# weights_pl = read_rds(str_c(data_dir, '07_places_weights.Rds'))
# weights_tr = read_rds(str_c(data_dir, '07_tracts_weights.Rds'))

tracts = tracts_sfl$ctd_60 %>%
    split(.$county) %>%
    keep(~ nrow(.) > 50)

## Weights and traces ----
construct_weights = function(sf, k = 3) {
    sf %>%
        st_centroid() %>%
        st_coordinates() %>%
        knearneigh(k = k) %>%
        knn2nb() %>%
        nb2listw(style = 'W')
}

weights_tr = map(tracts, construct_weights)
traces_tr = map(weights_tr, ~ trW(as(., 'CsparseMatrix')))

datal = list(places = tracts, 
             weights = weights_tr, 
             traces = traces_tr) %>%
    transpose()


## Fit models ----
reg_form = formula(log_w_use ~ hispanicP + noncitizensP +
                       blackP + indigenousP + asianP + womenP + 
                       childrenP + poverty_combP + 
                       ag_employedP + density_log10)

models = map(datal, 
             ~ lagsarlm(reg_form,
                        data = .$places,
                        listw = .$weights,
                        trs = .$traces,
                        type = 'Durbin',
                        zero.policy = NULL))


## Impacts ----
impacts = map2(models, traces_tr, ~ impacts(.x, tr = .y, R = 100))

impacts_draws = impacts %>%
    modify(~ .$sres) %>% 
    modify(~ .$total) %>%
    map(as.data.frame) %>%
    map(mutate, R = row_number()) %>%
    bind_rows(.id = 'county') %>%
    as_tibble()

impacts_long = gather(impacts_draws, key = variable, value = estimate, 
       hispanicP:density_log10)

quants = c(.05, .5, .95)
ggplot(impacts_long, aes(fct_rev(county), estimate)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # geom_violin(draw_quantiles = .5) +
    stat_summary(fun.ymin = function (x) quantile(x, 
                                                  probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, 
                                                  probs = quants[3])) +
    facet_wrap(~ variable, scales = 'free_x') +
    coord_flip() +
    theme_minimal()
