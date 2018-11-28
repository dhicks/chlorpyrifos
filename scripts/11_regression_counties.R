## This script constructs Spatial Durbin models for each county separately.  
## This approach was a late addition to the project.  Given the computational needs for the resampled models, I chose to calculate these models separately, rather than integrating them into script 10.  Note that this means changes to the model specifications used in script 10 are not necessarily reflected here.  

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

dataf = list(places = places_sfl, 
           tracts = tracts_sfl) %>%
    modify_depth(2, ~ split(., .$county)) %>%
    ## Recombine into a single second-order dataframe
    modify_depth(3, ~ tibble(data = list(.))) %>%
    modify_depth(2, bind_rows, .id = 'county') %>%
    modify(bind_rows, .id = 'ctd') %>%
    bind_rows(.id = 'geography') %>%
    ## Filter
    ## 30 obs for places, 50 for tracts
    mutate(nrow = map_int(data, nrow)) %>%
    filter(ifelse(geography == 'places', 
                  nrow >= 30, 
                  nrow >= 50))

# tracts = tracts_sfl$ctd_60 %>%
#     split(.$county) %>%
#     keep(~ nrow(.) > 50)
# 
# places = places_sfl$ctd_60 %>%
#     split(.$county) %>%
#     keep(~ nrow(.) >= 30)


## Weights and traces ----
construct_weights = function(sf, k = 3) {
    sf %>%
        st_centroid() %>%
        st_coordinates() %>%
        knearneigh(k = k) %>%
        knn2nb() %>%
        nb2listw(style = 'W')
}

dataf = dataf %>%
    mutate(weights = map(data, construct_weights), 
           traces = map(weights, ~ trW(as(., 'CsparseMatrix'))))


## Fit models ----
## This should be the same as in 10, only w/o the county-level variables
reg_form = formula(log_w_use ~ hispanicP + noncitizensP +
                       blackP + indigenousP + asianP + womenP + 
                       childrenP + poverty_combP + 
                       ag_employedP + density_log10)

# models = map(datal, 
#              ~ lagsarlm(reg_form,
#                         data = .$places,
#                         listw = .$weights,
#                         trs = .$traces,
#                         type = 'Durbin',
#                         zero.policy = NULL))

## ~15 sec
tic()
dataf = dataf %>% 
    # slice(1:5) %>%
    rowwise() %>%
    mutate(model = list(lagsarlm(reg_form, 
                            data = data, 
                            listw = weights, 
                            trs = traces, 
                            type = 'Durbin', 
                            zero.policy = NULL)), 
           impacts = list(impacts(model, 
                                  tr = traces, 
                                  R = 400))) %>%
    ungroup()
toc()


## Impacts ----
# impacts = map2(models, traces_pl, ~ impacts(.x, tr = .y, R = 100))

impacts_draws = dataf %>%
    select(geography, ctd, county, impacts) %>%
    rowwise() %>%
    mutate(draws = list(as.tibble(impacts$sres$total))) %>%
    unnest(draws, .id = 'co_model_idx')

impacts_long = gather(impacts_draws, key = variable, value = estimate, 
       hispanicP:density_log10)

# quants = c(.05, .5, .95)
# impacts_long %>%
#     # filter(ctd == 'ctd_60') %>%
#     ggplot(aes(fct_rev(county), estimate, color = geography)) +
#     geom_hline(yintercept = 0, linetype = 'dashed') +
#     # geom_violin(draw_quantiles = .5) +
#     stat_summary(fun.ymin = function (x) quantile(x, 
#                                                   probs = quants[1]),
#                  fun.y = median, 
#                  fun.ymax = function (x) quantile(x, 
#                                                   probs = quants[3])) +
#     facet_wrap(ctd ~ variable, scales = 'free_x', ncol = 10) +
#     coord_flip() +
#     theme_minimal()


## Outputs ----
write_rds(dataf, str_c(data_dir, '11_county_models.Rds'))
write_rds(impacts_draws, str_c(data_dir, '11_county_impacts.Rds'))
write_rds(impacts_long, str_c(data_dir, '11_county_impacts_long.Rds'))
