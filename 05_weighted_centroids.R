library(tidyverse)
library(tidycensus)
library(sf)

options(tigris_use_cache = TRUE)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

county_df = str_c(data_dir, '01_counties.Rda') %>%
    read_rds() %>%
    filter(study_area)

## Retrieve census blocks ----
## NB not fast
blocks_sf = get_decennial('block', 
                      variables = 'P0010001',
                      year = 2010, 
                      state = 'CA', 
                      county = county_df$county,
                      # county = 'Kern',
                      cache = TRUE,
                      geometry = TRUE) %>%
    st_transform('+proj=utm +zone=10 +datum=NAD83') %>%
    mutate(tract = str_extract(GEOID, '[0-9]{11}'), 
           centroid = st_centroid(geometry))

## Weighted centroids for tracts ----
## NB very fast
weighted_centroids_tracts = blocks_sf %>%
    ## Extract block centroid coordinates
    st_set_geometry('centroid') %>%
    bind_cols({st_coordinates(.) %>%
            as_tibble() %>%
            set_names(c('x', 'y'))}) %>%
    ## Drop geometry columns when we summarize()
    as.data.frame() %>%
    ## Weighted mean coordinates within each tract
    group_by(tract) %>%
    summarize(x = weighted.mean(x, value),
              y = weighted.mean(y, value)) %>%
    # summarize(x = mean(x), 
    #           y = mean(y)) %>%
    ## Results as sf
    st_as_sf(coords = c('x', 'y'), crs = st_crs(blocks_sf)) %>%
    arrange(tract)

tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds'))
tracts_sf %>%
    st_centroid() %>% 
    arrange(GEOID) %>%
    st_distance(weighted_centroids_tracts, by_element = TRUE) %>%
    units::drop_units() %>%
    {. / 1000} %>%
    tibble(centroid_shift = .) %>%
    ggplot(aes(centroid_shift)) + 
    geom_density() + 
    geom_rug()

write_rds(weighted_centroids_tracts, 
          str_c(data_dir, '05_w_cntr_tracts.Rds'))

## Weighted centroids for places ----
places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number())

## ~20s
blocks_by_places_sf = places_sf %>%
    select(GEOID) %>%
    st_join(blocks_sf, left = FALSE)

weighted_centroids_places = blocks_by_places_sf %>%
    ## Extract block centroid coordinates
    st_set_geometry('centroid') %>%
    bind_cols({st_coordinates(.) %>%
            as_tibble() %>%
            set_names(c('x', 'y'))}) %>%
    ## Drop geometry columns when we summarize()
    as.data.frame() %>%
    ## Weighted mean coordinates within each place
    group_by(GEOID = GEOID.x) %>%
    summarize(x = weighted.mean(x, value), 
              y = weighted.mean(y, value)) %>% 
    ## Drop 1 place (Silver City) whose blocks all have 0 population
    filter(!is.na(x)) %>%
    ## Results as sf
    st_as_sf(coords = c('x', 'y'), crs = st_crs(blocks_sf)) %>%
    arrange(GEOID)

places_sf %>%
    st_centroid() %>% 
    arrange(GEOID) %>%
    filter(GEOID %in% weighted_centroids_places$GEOID) %>%
    st_distance(weighted_centroids_places, by_element = TRUE) %>%
    units::drop_units() %>%
    {. / 1000} %>%
    tibble(centroid_shift = .) %>%
    ggplot(aes(centroid_shift)) + 
    geom_density() + 
    geom_rug()

write_rds(weighted_centroids_places, 
          str_c(data_dir, '05_w_cntr_places.Rds'))
