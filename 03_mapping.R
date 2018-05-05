## County shapefiles -----
## https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
counties_shp = readOGR(str_c(data_dir, 'cb_2015_us_county_20m'), 
                       'cb_2015_us_county_20m')
# proj4string(counties_shp)
counties_shp = spTransform(counties_shp, 
                           CRS("+proj=utm +zone=10 +datum=NAD83"))
counties_shp = subset(counties_shp, 
                      STATEFP == '06' & NAME %in% cv_counties)

## County-level use map -----
tm_shape(counties_shp) +
    tm_borders() +
    tm_fill('lightgrey') +
    tm_shape(chlor_sf) +
    tm_dots('total_use')



counties_total_use = counties_shp
counties_total_use@data = chlor_data %>%
    group_by(county) %>%
    summarize(total_lbs = sum(lbs_chm_used)) %>%
    right_join(counties_total_use@data, by = c('county'='NAME'))

tm_shape(counties_total_use) + 
    tm_polygons('total_lbs', style = 'cont', title = 'total use (lbs)') + 
    tm_text('county') + 
    tm_layout(main.title = 'Chlorpyrifos use, 2015', 
              scale = .5, main.title.size = 1)
save_tmap(tm = last_map(), filename = 'chlor_use.png', width = 3, height = 4, units = 'in')

