#' This script calculates distances and relative potential exposures between chlorpyrifos use locations and tracts and places
library(tidyverse)
library(tidycensus)
library(sf)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

ctd_df = tibble(ctd = c(1, 10, 30, 60, 90), 
                decay_rate = (1-.63)^(1/ctd)) %>%
    mutate(ctd = str_c('ctd_', ctd))

chlor_sf = read_rds(str_c(data_dir, '01_chlor_sf.Rds')) %>%
    mutate(log_total_use = log10(total_use), 
           chlor_idx = as.character(row_number()))
tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
    mutate(tracts_idx = row_number())
places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number())

tracts_centroids = read_rds(str_c(data_dir, '05_w_cntr_tracts.Rds'))
places_centroids = read_rds(str_c(data_dir, '05_w_cntr_places.Rds'))


## Places ----
## ~70s
tictoc::tic()
distances_pl = places_centroids %>% 
    ## Calculate distances between centroids and uses
    st_distance(chlor_sf) %>% 
    units::drop_units() %>% 
    ## Add indices
    `rownames<-`(places_centroids$GEOID) %>%
    `colnames<-`(chlor_sf$chlor_idx) %>%
    as_tibble(rownames = 'GEOID') %>%
    gather(key = 'chlor_idx', value = 'distance', -GEOID) %>% 
    mutate(## Distances in km
           distance = distance / 1000)
tictoc::toc()

## Identify uses within regions
distances_pl = st_within(chlor_sf, places_sf, sparse = FALSE) %>% 
    `rownames<-`(chlor_sf$chlor_idx) %>%
    `colnames<-`(places_sf$GEOID) %>%
    as_tibble(rownames = 'chlor_idx') %>%
    ## Make the whole thing long
    gather(key = 'GEOID', value = 'within', -chlor_idx) %>%
    inner_join(distances_pl)

## Calculate total weighted local use
## ~35s
tictoc::tic()
w_use_pl = distances_pl %>%
    inner_join(chlor_sf, by = 'chlor_idx') %>%
    select(GEOID, chlor_idx, year,
           distance, within, total_use) %>%
    crossing(ctd_df) %>%
    mutate(decay_coef = ifelse(within, 1, decay_rate^distance)) %>% 
    group_by(ctd, GEOID) %>%
    summarize(w_use = sum(total_use * decay_coef, 
                            na.rm = TRUE))
tictoc::toc()

# left_join(places_sf, w_use_pl) %>%
#     ggplot(aes(w_use)) + 
#     geom_density() + 
#     geom_rug() + 
#     facet_wrap(~ ctd + county, scales = 'free')
# library(tmap)
# left_join(places_sf, w_use_pl) %>%
#     tm_shape() +
#     tm_polygons(col = 'w_use', style = 'cont', 
#                 popup.vars = c('w_use', 'GEOID'))

## Save output
write_rds(distances_pl, str_c(data_dir, '06_dist_places.Rds'))
write_rds(w_use_pl, str_c(data_dir, '06_w_use_places.Rds'))


## Tracts ----
## ~111s
tictoc::tic()
distances_tr = tracts_centroids %>%
    ## Distances
    st_distance(chlor_sf) %>%
    units::drop_units() %>%
    ## Add indices
    `rownames<-`(tracts_centroids$tract) %>%
    `colnames<-`(chlor_sf$chlor_idx) %>%
    as_tibble(rownames = 'GEOID') %>%
    gather(key = 'chlor_idx', value = 'distance', -GEOID) %>% 
    mutate(## Distances in km
        distance = distance / 1000)
tictoc::toc()

distances_tr = st_within(chlor_sf, tracts_sf, sparse = FALSE) %>%
    `rownames<-`(chlor_sf$chlor_idx) %>%
    `colnames<-`(tracts_sf$GEOID) %>%
    as_tibble(rownames = 'chlor_idx') %>%
    ## Make the whole thing long
    gather(key = 'GEOID', value = 'within', -chlor_idx) %>%
    inner_join(distances_tr)

## Calculate total potential exposure
## ~146s
## NB >14GB memory
tictoc::tic()
w_use_tr = distances_tr %>%
    inner_join(chlor_sf, by = 'chlor_idx') %>%
    select(GEOID, chlor_idx, year,
           distance, within, total_use) %>%
    crossing(ctd_df) %>%
    mutate(decay_coef = ifelse(within, 1, decay_rate^distance)) %>%
    group_by(ctd, GEOID) %>%
    summarize(w_use = sum(total_use * decay_coef, 
                          na.rm = TRUE))
tictoc::toc()

# left_join(tracts_sf, w_use_tr) %>%
#     ggplot(aes(w_use)) +
#     geom_density() +
#     geom_rug() +
#     facet_wrap(~ ctd, scales = 'free')
# library(tmap)
# w_use_tr %>%
#     filter(ctd == 'ctd_1') %>%
#     left_join(tracts_sf, .) %>%
#     tm_shape() +
#     tm_fill(col = 'w_use', style = 'cont',
#                 popup.vars = c('w_use', 'GEOID'))

write_rds(distances_tr, str_c(data_dir, '06_dist_tracts.Rds'))
write_rds(w_use_tr, str_c(data_dir, '06_w_use_tracts.Rds'))
