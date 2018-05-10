#' This script aggregates chlorpyrifos use around tracts and places
library(tidyverse)
library(tidycensus)
library(sf)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

chlor_sf = read_rds(str_c(data_dir, '01_chlor_sf.Rds')) %>%
    mutate(log_total_use = log10(total_use), 
           chlor_idx = row_number())
tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
    mutate(tracts_idx = row_number())
places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number())

## Tracts ----
## ~68s
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

write_rds(distances_tr, str_c(data_dir, '04_dist_tracts.Rds'))

# distances_tr %>%
#     inner_join(tracts_sf) %>%
#     inner_join(chlor_sf, by = 'chlor_idx') %>%
#     select(tracts_idx, chlor_idx, year, 
#            distance, total_use#, 
#            # geometry = geometry.x
#     ) %>%
#     # st_as_sf() %>%
#     mutate(decay_coef_60 = distance^-.24, 
#            weighted_use_60 = total_use * decay_coef_60, 
#            decay_coef_30 = distance^-.29, 
#            weighted_use_30 = total_use * decay_coef_30) %>%
#     group_by(tracts_idx, year) %>%
#     summarize(weighted_use_60 = sum(weighted_use_60), 
#               weighted_use_30 = sum(weighted_use_30)) %>%
#     gather(key = calculation, value = weighted_use, 
#            -tracts_idx, -year) %>%
#     ggplot(aes(weighted_use, color = year)) + 
#     # geom_violin(draw_quantiles = .5) + 
#     geom_density() + #scale_x_log10() +
#     facet_wrap(~ calculation)

## Places ----
## ~37s
tictoc::tic()
distances_pl = places_sf %>%
    st_centroid() %>%
    st_distance(chlor_sf) %>%
    units::drop_units() %>%
    as_tibble() %>%
    mutate(places_idx = row_number()) %>%
    gather(key = chlor_idx, value = distance, -places_idx) %>%
    mutate(chlor_idx = as.integer(str_extract(chlor_idx, '[0-9]+')), 
           distance = distance / 1000)
tictoc::toc()

distances_pl = st_within(chlor_sf, places_sf, sparse = FALSE) %>%
    as_tibble() %>%
    mutate(chlor_idx = row_number()) %>%
    gather(key = places_idx, value = within, -chlor_idx) %>%
    mutate(places_idx = as.integer(str_extract(places_idx, '[0-9]+'))) %>%
    right_join(distances_pl)

write_rds(distances_pl, str_c(data_dir, '04_dist_places.Rds'))

# distances_pl %>%
#     inner_join(places_sf) %>%
#     inner_join(chlor_sf, by = 'chlor_idx') %>%
#     select(places_idx, chlor_idx, year,
#            distance, within, total_use#,
#            # geometry = geometry.x
#            ) %>%
#     # st_as_sf() %>%
#     mutate(decay_coef_60 = ifelse(within, 1, .0062^distance),
#            weighted_use_60 = total_use * decay_coef_60, 
#            decay_coef_30 = ifelse(within, 1, .012^distance),
#            weighted_use_30 = total_use * decay_coef_30,
#            decay_coef_1 = ifelse(within, 1, .37^distance),
#            weighted_use_1 = total_use * decay_coef_1) %>%
#     group_by(places_idx) %>%
#     summarize_at(vars(contains('weighted_use')),
#                  funs(sum)) %>%
#     gather(key = calculation, value = weighted_use,
#            -places_idx) %>% 
#     mutate(weighted_use_trimmed = ifelse(weighted_use > 1, 
#                                        weighted_use, 
#                                        0)) %>%
#     ggplot(aes(weighted_use_trimmed)) +
#     geom_density() +
#     # scale_x_log10() +
#     facet_wrap(~ calculation, scales = 'free_y')
