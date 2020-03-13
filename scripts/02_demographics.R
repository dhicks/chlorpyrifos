## tidycensus uses the Census API. Obtain a key at <https://api.census.gov/data/key_signup.html>
library(tidyverse)
library(tidycensus)
library(sf)

# data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'
data_dir = '../data/'

county_df = str_c(data_dir, '01_counties.Rda') %>%
    read_rds() %>%
    filter(study_area)

## Variables of interest -----
## List of ACS variables
# vars_df = load_variables(2015, "acs5", cache = TRUE)

vars_of_interest = c(
         'total_pop' = 'B01003_001', # total population
         ## Gender
         'women' = 'B01001_026', # women/girls, all ages
         'men' = 'B01001_002', # men/boys, all ages
         ## Race and ethnicity
         'white' = 'B01001H_001', # white-alone non-hispanic
         'black' = 'B01001B_001', # black-alone
         'indigenous' = 'B01001C_001', # indigenous-alone
         'asian' = 'B01001D_001', # asian-alone
         'hispanic' = 'B01001I_001', # hispanic or latino
         ## Birth x citizenship
         # 'total_citizenship' = 'B05002_001', # total for birth x citizenship
         'noncitizens' = 'B05002_021', # foreign-born noncitizens
         ## Children
         # 'total_age' = 'B06001_001', # total for age
         'children' = 'B06001_002', # children under 5 yo
         ## Poverty
         # 'total_poverty' = 'C17002_001', # total for poverty
         'extreme_poverty' = 'C17002_002', # IP ratio < .5
         'poverty' = 'C17002_003', # IP ratio .5-.99
         ## Hispanic poverty
         # 'total_hisp_poverty' = 'B17001I_001', # Hispanic total
         'hisp_poverty' = 'B17001I_002', # Hispanic poverty
         ## Migration
         ## Not working? 
         # 'total_migration' = 'B99072_001', # total for migration
         # 'migrant' = 'B99072_002' # different house
         ## Employed in agriculture
         'total_employed' = 'C24050_001', # civilian, employed, >= 16 yo
         'ag_employed' = 'C24050_002' # employed in ag
         )

## Tracts -----
tracts_sf = get_acs(geography = 'tract', 
                     variables = vars_of_interest,
                     year = 2015,
                     state = 'CA', 
                     county = county_df$county, 
                     cache = TRUE, output = 'wide',
                     geometry = TRUE) %>% 
    st_transform('+proj=utm +zone=10 +datum=NAD83') %>%
    mutate(area = st_area(.),
           area = units::set_units(area, 'km^2'),
           densityE = total_popE / units::drop_units(area),
           densityM = total_popM / units::drop_units(area),
           density_log10 = log10(densityE))

write_rds(tracts_sf, str_c(data_dir, '02_tracts_sf.Rds'))

## Places -----
options(tigris_use_cache = TRUE)
counties_unfltd_shp = tigris::counties(state = 'CA')
counties_shp = subset(counties_unfltd_shp, NAME %in% county_df$county)

places_shp_unfltd = tigris::places('CA')
places_sf = raster::intersect(places_shp_unfltd, counties_shp) %>%
    st_as_sf() %>%
    st_transform('+proj=utm +zone=10 +datum=NAD83')

places_df = get_acs(geography = 'place', 
                    variables = vars_of_interest, 
                    year = 2015, 
                    state = 'CA', 
                    # county = county_df$county, 
                    cache = TRUE, output = 'wide', 
                    geometry = FALSE)

places_sf = places_sf %>% 
    inner_join(places_df, by = c('GEOID.1' = 'GEOID')) %>%
    select(GEOID = GEOID.1, 
           name = NAME, 
           county = NAME.2,
           total_popE:ag_employedM) %>%
    mutate(area = st_area(.),
           area = units::set_units(area, 'km^2'),
           densityE = total_popE / units::drop_units(area), 
           densityM = total_popM / units::drop_units(area),
           density_log10 = log10(densityE))

write_rds(places_sf, str_c(data_dir, '02_places_sf.Rds'))

# library(tmap)
# tm_shape(tracts_sf) + 
#     tm_fill(col = 'MAP_COLORS') +
#     tm_shape(places_sf) + 
#     tm_borders(col = 'black') +
#     tm_text('NAME', size = .6)
# 
# ## Places cover 85% of the population
# sum(places_sf$total_popE) / sum(tracts_sf$total_popE)
# ## But <8% of the area
# sum(places_sf$area) / sum(tracts_sf$area)
# 
