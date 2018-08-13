## Spatial Durbin models in Stan
## Setup ----
library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)

library(rstan)
options(mc.cores = parallel::detectCores())

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

## Spatial weights as sparse matrices
k = 3
weights_pl = places_sf %>%
    filter(ctd == 'ctd_30') %>%
    st_centroid %>%
    st_coordinates() %>%
    knearneigh(k = k) %>%
    knn2nb() %>%
    nb2listw(style = 'W') %>%
    as('RsparseMatrix')


## Stan ----
response = 'log_w_use'
vars = c('hispanicP', 'blackP', 'indigenousP',
         'asianP', 'childrenP', 'poverty_combP',
         'ag_employedP', 'log_densityE')
# vars = c('hispanicP', 'ag_employedP', 'log_densityE')
dataf = places_sf %>%
    filter(ctd == 'ctd_60') %>%
    as_tibble() %>%
    select(one_of(response, vars))


parsed = stanc(file = '12_sdm.stan', verbose = TRUE)
compiled = stan_model(stanc_ret = parsed, verbose = FALSE)
pars = c('alpha', 'beta', 'rho', 'gamma')

tictoc::tic()
samples = sampling(compiled, data = list(N = nrow(dataf), 
                                         p = ncol(dataf) - 1,
                                         y = dataf[[1]], 
                                         X = dataf[-1],
                                         # W = as.matrix(W_3nn)
                                         Wn = length(weights_pl@j), 
                                         Wrow = length(weights_pl@p),
                                         Wv = weights_pl@j+1,
                                         Wu = weights_pl@p+1,
                                         Ww = weights_pl@x),
                   chains = 2, iter = 2000, 
                   control = list(adapt_delta = 0.9, 
                                  max_treedepth = 15))
tictoc::toc()

# summary(samples, pars = pars)$summary
# pairs(samples, pars = pars)
