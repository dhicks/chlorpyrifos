#' This script cleans the PUR data for chlorpyrifos use in the Central Valley.  Output files are in `shp` (Esri Shapefile) and `Rds` (R-specific serialization) formats. 
#' Important notes:
#' 
#' - Output contains chlorpyrifos use for all California counties.  Counties to be analyzed downstream (Central Valley counties other than Sacramento) are designed with the variable `study_area`. 
#' - Output contains only total annual use of chemical (in pounds) at each section for 2011-2015. 
#' - Uses are mapped to section centroids rather than sections or fields. 
#' - Chlorpyrifos use is identified by matching the string 'chlorpyrifos' in the chemical list.  This matches chlorpyrifos as well as chlorpyrifos-methyl and chlorpyifos oxon.  In 2015, only 2/11k uses were of chlorpyrifos-methyl or chlorpyrifos oxon.  
#' - Shapefiles are in UTM 10 projection. 

library(tidyverse)
library(sf)
library(tictoc)
library(vroom)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

counties_file = str_c(data_dir, '01_counties.Rda')
if (!file.exists(counties_file)) {
    study_area = c('Shasta', 'Tehama', 'Glenn', 'Butte', 
                   'Colusa', 'Yuba', 'Sutter',
                   'Yolo', 'Solano', 'San Joaquin', 'Stanislaus', 
                   'Merced', 'Madera', 'Fresno', 'Kings', 
                   'Tulare', 'Kern')
    pur_county_file = str_c(data_dir, 'pur2015/county.txt')
    county_df = pur_county_file %>%
        read_csv() %>%
        rename(county = couty_name) %>% 
        mutate(county = str_to_title(county), 
               study_area = county %in% study_area)
    
    write_rds(county_df, counties_file)
} else {
    county_df = read_rds(counties_file)
}

## Pesticide use data -----
# pur_years = 11:15
pur_years = 11

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

# ~172 s
tic()
pur_data_unfltd = pur_files %>%
    map(vroom, 
        delim = ',',
        col_types = cols(range = 'c', 
                         township = 'c', section = 'c', county_cd = 'c', 
                         grower_id = 'c', license_no = 'c', 
                         applic_time = 'c', section = 'c', 
                         site_loc_id = 'c',
                         acre_planted = 'n', acre_treated = 'n')) %>%
    bind_rows(.id = 'year') %>%
    full_join(county_df)
toc()

## Chemical data -----
chem_file = str_c(data_dir, 'pur2015/chemical.txt')
chem_df = read_csv(chem_file)

## These appear to be the same across years, at least for chlorpyrifos
# chem_df2 = read_csv(str_c(data_dir, 'pur2011/chemical.txt'))

chlor_data_uncleaned = chem_df %>%
    filter(str_detect(chemname, 'CHLORPYRIFOS')) %>%
    inner_join(pur_data_unfltd) %>%
    ## 1 entry w/ NA lbs_chm_used
    filter(!is.na(lbs_chm_used))

## TODO: move this below error checking
## Site/commodity code ----
site_file = str_c(data_dir, 'pur2015/site.txt')
site_df = read_csv(site_file, 
                   col_types = cols(site_code = 'c'))
chlor_data_uncleaned %>% 
    left_join(site_df, by = 'site_code') %>% 
    group_by(site_name) %>% 
    summarize(lbs_total_use = sum(lbs_chm_used), 
              median_rate = median(lbs_chm_used / acre_treated)) %>% 
    arrange(desc(lbs_total_use))

## TODO: move below error checking
## Month used ----
chlor_data_uncleaned %>% 
    mutate(date = lubridate::mdy(applic_dt), 
           month = lubridate::month(date)) %>% 
    group_by(month) %>% 
    summarize(lbs_total_use = sum(lbs_chm_used)) %>% 
    ggplot(aes(as.factor(month), lbs_total_use, 
               group = 1L)) +
    geom_point() +
    geom_line()


## Error checking -----
chlor_errors = pur_years %>%
    str_c('pur20', ., '/errors20', ., '.txt') %>%
    str_c(data_dir, .) %>%
    set_names(str_c('20', pur_years)) %>%
    map_dfr(read_csv, .id = 'year') %>%
    inner_join(chlor_data_uncleaned)
## Table of types of errors
count(chlor_errors, error_code, error_description)

## Changes were only made for inconsistent values (codes 43, 47), usually identifiers
## So *not* when the record might be a duplicate (80), the application year might be wrong (21), or the units of measure might be wrong (38)
str_c(data_dir, 'pur2015/changes2015.txt') %>%
    read_csv() %>%
    inner_join(chlor_errors)

## Outlier correction (code 75)
## It looks like these adjustments have already been made in the datafile
filter(chlor_errors, error_code == 75) %>%
    select(year, comments, lbs_prd_used, acre_treated)
str_c(data_dir, 'pur2011/changes2011.txt') %>%
    read_csv() %>%
    inner_join(chlor_errors) %>%
    filter(error_code == 75)

## When the record might be a duplicate (80), duplicate_set doesn't actually identify the duplicates? 
## Or does in some cases, and in other cases they've already been removed? 
chlor_errors %>%
    filter(error_code == 80) %>%
    count(duplicate_set, year) %>%
    arrange(desc(n)) #%>% View

## In 2015, grower # 50155005245 has a potential duplicate (use # 181070, row 5 below)
## But there doesn't seem to be another entry with the same lbs_chm_used
chlor_data_uncleaned %>%
    filter(grower_id == '50155005245') #%>% View
## So I guess duplicates have already been deleted?  

## This grower seems to have some issues with their paperwork
chlor_errors %>%
    filter(grower_id == '50155005245')
## Uses 94056 and 94061 seem to be identical, but have different site_loc_id values

## OTOH, consider these two from 20122
chlor_data_uncleaned %>%
    filter(use_no %in% c(2468421, 5093585))# %>% View
## These are identical, except that the application time is 0440 vs 1640 and some of the record identifiers are different.  

## So we'll remove the second record when duplicate_set appears twice
chlor_data = chlor_errors %>%
    filter(error_code == 80) %>%
    count(duplicate_set, year) %>%
    arrange(desc(n)) %>%
    filter(n > 1) %>%
    ## Rejoin w/ chlor_errors to get the use numbers
    inner_join(chlor_errors) %>%
    select(duplicate_set, year, use_no) %>%
    ## Get the second record
    group_by(duplicate_set, year) %>%
    summarize(use_no = last(use_no)) %>%
    ## Anti-join to chlor_data_uncleaned
    anti_join(chlor_data_uncleaned, .)

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
    summarize(total_use = sum(lbs_chm_used) * 0.45359237) %>%
    inner_join(sections_sf, ., 
               by = c('CO_MTRS' = 'comtrs'))

write_rds(chlor_sf, str_c(data_dir, '01_chlor_sf.Rds'))
st_write(chlor_sf, str_c(path.expand(data_dir), 
                         'chlorpyrifos/', 
                         'chlorpyrifos.shp'), 
         delete_layer = TRUE)


