library(tidyverse)
library(sf)
library(tmap)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

## County borders -----
## https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
counties_sf = read_sf(str_c(data_dir, 'cb_2015_us_county_20m'), 
                       'cb_2015_us_county_20m') %>%
    filter(STATEFP == '06') %>%
    st_transform('+proj=utm +zone=10 +datum=NAD83')

## Chlorpyrifos use ----
chlor_sf = read_rds(str_c(data_dir, '01_chlor_sf.Rds'))

## Tracts and places ----
tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds'))
places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds'))
    

## Use map -----
tm_shape(counties_sf) + 
    tm_borders(col = 'black') +
tm_shape(tracts_sf) +
    tm_polygons(col = 'blue', alpha = .25) +
tm_shape(places_sf) +
    tm_polygons(col = 'yellow', alpha = .75) +
tm_shape(subset(chlor_sf, year == 2015)) +
    tm_dots(col = 'total_use', palette = 'Reds') +
tm_scale_bar(position = c('left', 'bottom'))


## Distribution of use ----
## I'd say this is close enough to lognormal
## Annual values
ggplot(chlor_sf, aes(log10(total_use), color = year)) + 
    geom_density()

## Aggregated
chlor_sf %>%
    group_by(CO_MTRS) %>%
    summarize(total_use = sum(total_use)) %>%
    ggplot(aes(log10(total_use))) + 
    geom_density()
