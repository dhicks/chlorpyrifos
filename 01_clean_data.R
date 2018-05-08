## TODO: multiple years
##  + check for error 75 (outliers)

#' This script cleans the PUR data for chlorpyrifos use in the Central Valley.  Output files are in `shp` (Esri Shapefile) and `Rds` (R-specific serialization) formats. 
#' Important notes:
#' 
#' - Output contains only uses in the counties listed in `cv_counties`
#' - Output contains only total use of chemical (in pounds) at each section.
#' - Uses are mapped to section centroids rather than sections or fields. 
#' - Chlorpyrifos use is identified by matching the string 'chlorpyrifos' in the chemical list.  This matches chlorpyrifos as well as chlorpyrifos-methyl and chlorpyifos oxon.  In 2015, only 2/11192 uses were of chlorpyrifos-methyl or chlorpyrifos oxon.  
#' - Shapefiles are in UTM 10 projection. 

library(tidyverse)
library(sf)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

counties_file = str_c(data_dir, '01_counties.Rda')
if (!file.exists(counties_file)) {
  counties_of_interest = c('Shasta', 'Tehama', 'Glenn', 'Butte', 
                           'Colusa', 'San Joaquin', 'Stanislaus', 
                           'Merced', 'Madera', 'Fresno', 'Kings', 
                           'Tulare', 'Kern')
  pur_county_file = str_c(data_dir, 'pur2015/county.txt')
  county_df = pur_county_file %>%
    read_csv() %>%
    rename(county = couty_name) %>% 
    filter(county %in% toupper(counties_of_interest)) %>%
    mutate(county = str_to_title(county))
  
  write_rds(county_df, counties_file)
} else {
  county_df = read_rds(counties_file)
}

## Pesticide use data -----
pur_years = 11:15

pur_files = county_df %>%
    crossing(pur_years) %>%
    mutate(dir = str_c('pur20', pur_years), 
           file = str_c('udc', pur_years, '_', county_cd, '.txt'), 
           path = str_c(dir, '/', file), 
           fullpath = str_c(data_dir, path), 
           pur_years = str_c('20', pur_years)) %>%
    pull(fullpath)
names(pur_files) <- rep(str_c('20', pur_years), 
                        length.out = length(pur_files))

pur_data_unfltd = pur_files %>%
    map(read_csv) %>%
    map(mutate_at, vars(range, township, county_cd, 
                        grower_id, license_no),
        'as.character') %>%
    bind_rows(.id = 'year') %>%
    full_join(county_df)

## Chemical data -----
chem_file = str_c(data_dir, 'pur2015/chemical.txt')
chem_df = read_csv(chem_file)

## These appear to be the same across years, at least for chlorpyrifos
# chem_df2 = read_csv(str_c(data_dir, 'pur2011/chemical.txt'))

chlor_data = chem_df %>%
    filter(str_detect(chemname, 'CHLORPYRIFOS')) %>%
    inner_join(pur_data_unfltd) %>%
    ## 1 entry w/ NA lbs_chm_used
    filter(!is.na(lbs_chm_used))

## Error checking -----
chlor_errors = read_csv(str_c(data_dir, 'pur2015/errors2015.txt')) %>%
    inner_join(chlor_data)
## Table of types of errors
count(chlor_errors, error_code, error_description)

## Changes were only made for inconsistent values (codes 43, 47), usually identifiers
## So *not* when the record might be a duplicate (80), the application year might be wrong (21), or the units of measure might be wrong (38)
str_c(data_dir, 'pur2015/changes2015.txt') %>%
    read_csv() %>%
    inner_join(chlor_errors)

## When the record might be a duplicate (80), duplicate_set doesn't actually identify the duplicates? 
chlor_errors %>%
    filter(year == 2015) %>%
    count(duplicate_set) %>%
    arrange(desc(n))

## Grower # 50155005245 has a potential duplicate (use # 181070, row 5 below)
## But there doesn't seem to be another entry with the same lbs_chm_used
chlor_data %>%
    filter(grower_id == '50155005245')
## So I guess duplicates have already been deleted?  

## This grower seems to have some issues with their paperwork
chlor_errors %>%
    filter(grower_id == '50155005245')
## Uses 94056 and 94061 seem to be identical, but have different site_loc_id values

## Based on this, I'm going to infer that duplicates have already been removed and I don't need to do any further cleaning

## Section shapefiles -----
## ~10 seconds
sections_unfltd = read_sf(str_c(data_dir, '/plsnet'), 'plsnet')

sections_sf = sections_unfltd %>%
    inner_join(county_df, by = c('COUNTY_CD' = 'county_cd')) %>%
    st_transform('+proj=utm +zone=10 +datum=NAD83') %>%
    st_centroid()

## Combine section shapefiles w/ chlorpyrifos use data ----
chlor_sf = chlor_data %>%
    group_by(year, comtrs) %>%
    summarize(total_use = sum(lbs_chm_used)) %>%
    inner_join(sections_sf, ., 
                      by = c('CO_MTRS' = 'comtrs'))

write_rds(chlor_sf, str_c(data_dir, '01_chlor_sf.Rds'))
st_write(chlor_sf, str_c(path.expand(data_dir), 
                         'chlorpyrifos_2015/', 
                         'chlorpyrifos_2015.shp'), 
         delete_layer = TRUE)


