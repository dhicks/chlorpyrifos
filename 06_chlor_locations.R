#' This script calculates distances and relative potential exposures between chlorpyrifos use locations and tracts and places
library(tidyverse)
library(tidycensus)
library(sf)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

ctd_df = tibble(ctd = c(1, 10, 30, 60, 90), 
                decay_rate = (1-.63)^(1/ctd)) %>%
    mutate(ctd = str_c('ctd_', ctd), 
           ctd = fct_inorder(ctd))

chlor_sf = read_rds(str_c(data_dir, '01_chlor_sf.Rds')) %>%
    mutate(log_total_use = log10(total_use), 
           chlor_idx = row_number())
tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
    mutate(tracts_idx = row_number())
places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number())


## Places ----
## ~65s
tictoc::tic()
## Calculate distances between centroids and uses
distances_pl = places_sf %>%
    st_centroid() %>%
    st_distance(chlor_sf) %>%
    units::drop_units() %>%
    as_tibble() %>%
    mutate(places_idx = row_number()) %>%
    gather(key = chlor_idx, value = distance, -places_idx) %>%
    mutate(chlor_idx = as.integer(str_extract(chlor_idx, '[0-9]+')), 
           ## Distances in km
           distance = distance / 1000)
tictoc::toc()

## Identify uses within regions, and make the whole thing long
distances_pl = st_within(chlor_sf, places_sf, sparse = FALSE) %>%
    as_tibble() %>%
    mutate(chlor_idx = row_number()) %>%
    gather(key = places_idx, value = within, -chlor_idx) %>%
    mutate(places_idx = as.integer(str_extract(places_idx, '[0-9]+'))) %>%
    right_join(distances_pl)

## Calculate total potential exposure
## ~35s
tictoc::tic()
pot_exp_pl = distances_pl %>%
    inner_join(chlor_sf, by = 'chlor_idx') %>%
    select(places_idx, chlor_idx, year,
           distance, within, total_use) %>%
    crossing(ctd_df) %>%
    mutate(decay_coef = ifelse(within, 1, decay_rate^distance)) %>%
    group_by(ctd, places_idx) %>%
    summarize(pot_exp = sum(total_use * decay_coef, 
                            na.rm = TRUE)) %>%
    group_by(ctd) %>%
    mutate(pot_exp_z = (pot_exp - mean(pot_exp, na.rm = TRUE)) / 
               sd(pot_exp, na.rm = TRUE)) %>%
    ungroup() %>%
    select(-pot_exp) %>%
    spread(key = ctd, value = pot_exp_z, fill = 0)
tictoc::toc()

## Save output
write_rds(distances_pl, str_c(data_dir, '05_dist_places.Rds'))
write_rds(pot_exp_pl, str_c(data_dir, '05_pot_exp_places.Rds'))


## Tracts ----
## ~111s
tictoc::tic()
distances_tr = tracts_sf %>%
    st_centroid() %>%
    st_distance(chlor_sf) %>%
    units::drop_units() %>%
    as_tibble() %>%
    mutate(tracts_idx = row_number()) %>%
    gather(key = chlor_idx, value = distance, -tracts_idx) %>%
    mutate(chlor_idx = as.integer(str_extract(chlor_idx, '[0-9]+')), 
           distance = distance / 1000)
tictoc::toc()

distances_tr = st_within(chlor_sf, tracts_sf, sparse = FALSE) %>%
    as_tibble() %>%
    mutate(chlor_idx = row_number()) %>%
    gather(key = tracts_idx, value = within, -chlor_idx) %>%
    mutate(tracts_idx = as.integer(str_extract(tracts_idx, '[0-9]+'))) %>%
    right_join(distances_tr)

## Calculate total potential exposure
## ~146s
## NB ~14GB memory
tictoc::tic()
pot_exp_tr = distances_tr %>%
    inner_join(chlor_sf, by = 'chlor_idx') %>%
    select(tracts_idx, chlor_idx, year,
           distance, within, total_use) %>%
    crossing(ctd_df) %>%
    mutate(decay_coef = ifelse(within, 1, decay_rate^distance)) %>%
    group_by(ctd, tracts_idx) %>%
    summarize(pot_exp = sum(total_use * decay_coef, 
                            na.rm = TRUE)) %>%
    group_by(ctd) %>%
    mutate(pot_exp_z = (pot_exp - mean(pot_exp, na.rm = TRUE)) / 
               sd(pot_exp, na.rm = TRUE)) %>%
    ungroup() %>%
    select(-pot_exp) %>%
    spread(key = ctd, value = pot_exp_z, fill = 0)
tictoc::toc()

write_rds(distances_tr, str_c(data_dir, '05_dist_tracts.Rds'))
write_rds(pot_exp_tr, str_c(data_dir, '05_pot_exp_tracts.Rds'))
