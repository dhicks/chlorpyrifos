#' This document conducts an EDA on census data only, first at the tract level and then at the places level.  
library(tidyverse)
library(sf)

library(tidycensus)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

#' # Tracts #
tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
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
                                     poverty_combM, total_popM)) %>%
    ## Population densities
    mutate_at(vars(whiteE, blackE, indigenousE, asianE, hispanicE, 
                   noncitizensE, childrenE, poverty_combE), 
              funs(D = . / units::drop_units(area)))
## Gives a warning about NaNs; 
## but there aren't any in the output
# as.data.frame() %>%
# select(-geometry) %>%
# transmute_all(funs(is.nan)) %>%
# summarize_all(sum)

glimpse(tracts_sf)

## Distributions of proportions ----
tracts_sf %>%
    as.data.frame() %>%
    select(ends_with('P')) %>%
    gather(key = variable, value) %>%
    ggplot(aes(value, color = variable, fill = variable)) +
    geom_density() +
    geom_rug() +
    facet_wrap(~ variable, scales = 'free') +
    scale_x_continuous(labels = scales::percent_format())

#' There are a few tracts with a modest proportion of Asian and Black residents ($> 20%$); but only a few.  Almost no tracts have more than 5% Indigenous residents, and none have more than about 20%.  Children are typically ~5-12% of the population.  The poverty rate varies dramatically, with a median somewhere around 30% and some values above 50%.  Hispanic and White proportions are the most diverse.  

tracts_sf %>%
    mutate(w_plus_h = whiteP + hispanicP) %>%
    ggplot(aes(w_plus_h)) + 
    geom_density() + 
    geom_rug()

#' In most tracts, a supermajority of people are either White or Hispanic.  

ggplot(tracts_sf, aes(hispanicP, poverty_combP)) + 
    geom_point() +
    geom_smooth(method = 'lm')

#' A greater Hispanic proportion is associated with a greater poverty rate.  

## Correlations ----
tracts_sf %>%
    as_tibble() %>%
    select(densityE, whiteP, blackP, childrenP, hispanicP, indigenousP, noncitizensP, poverty_combP, whiteP) %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'var1') %>%
    as_tibble() %>%
    gather(key = 'var2', value = 'cor', -var1) %>%
    mutate(cor.print = round(cor, digits = 1)) %>%
    ggplot(aes(var1, var2, fill = cor, label = cor.print)) + 
    geom_tile() +
    geom_text() +
    scale_fill_gradient2()
    
#' White proportion has moderate to very strong negative corelations with every other variable (except Indigenous).  Very strong correlation between Hispanic and noncitizen proportion.  Moderate correlations between each pair of poverty, children, noncitizens, and Hispanic.  

ggplot(tracts_sf, aes((densityE), hispanicP)) + 
    geom_point() +
    geom_smooth()

#' No indication of a relationship between population density and Hispanic, linear or nonlinear.  


## White/Hispanic segregation ----
## Evenness/dissimilarity
tracts_sf %>%
    as.data.frame() %>%
    mutate(w_h_dissim = abs(hispanicE/sum(hispanicE) - whiteE/sum(whiteE))) %>%
    summarize(w_h_dissim = .5 * sum(w_h_dissim))
#' Evenness is moderate, at 47%

ggplot(tracts_sf, aes(abs(hispanicP - whiteP))) + 
    geom_density() +
    geom_rug()

## Exposure/interaction
tracts_sf %>%
    as.data.frame() %>%
    mutate(w_h_exposure = abs(whiteE/sum(whiteE) * hispanicE / total_popE), 
           h_w_exposure = abs(hispanicE/sum(hispanicE) * whiteE/total_popE)) %>%
    summarize(w_h_exposure = sum(w_h_exposure), 
              h_w_exposure = sum(h_w_exposure))
#' Exposure is moderate-low, at about 30% in both directions

## Correlation
ggplot(tracts_sf, aes(hispanicP, whiteP)) + 
    geom_point()
with(tracts_sf, cor(hispanicP, whiteP))

#' Very strong negative correlation between the two proportions. But I guess, in this kind of two-group context, we would get a very strong negative correlation even if dissimilarity were low and exposure were high.  

tracts_sf %>%
    as.data.frame() %>%
    as_tibble() %>%
    summarize_at(vars(whiteE, hispanicE), funs(sum)) %>%
    mutate(hw_ratio = hispanicE / whiteE)
#' About 20% more Hispanics than Whites

## Spatial weights ----
## We'll construct KNN for K=3-8
library(spdep)

weights = 3:8 %>%
    set_names() %>%
    map(~ {tracts_sf %>%
            st_centroid() %>%
            st_coordinates() %>%
            knearneigh(k = .x) %>%
            knn2nb() %>%
            nb2listw(style = 'W')})

## Moran's I ----
moran.i = function(vec, weights, ...) {
    moran.test(vec, 
               weights, ...)$estimate['Moran I statistic']
}

## Moran's I for overall density
weights %>%
    map(~moran.i(log10(tracts_sf$densityE), .)) %>%
    unlist()
#' Moderate population clustering, .41-.45

weights %>%
tibble(weights = ., k = names(.)) %>%
    crossing(tibble(variable = c('whiteE_D', 
                                 'blackE_D', 
                                 'indigenousE_D', 
                                 'asianE_D', 
                                 'hispanicE_D', 
                                 'noncitizensE_D', 
                                 'childrenE_D', 
                           'poverty_combE_D'))) %>%
    rowwise() %>%
    mutate(var_value = {tracts_sf %>% 
            as.data.frame() %>%
            pull(variable) %>%
            {. + 1} %>% log10() %>%
            list()}, 
           moran_i = moran.i(var_value, weights)) %>%
    select(k, variable, moran_i) %>%
    arrange(desc(moran_i)) %>%
    mutate(variable = fct_inorder(variable)) %>%
    ggplot(aes(variable, moran_i, color = k, group = k)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = c(.41, .45), linetype = 'dashed') +
    coord_flip()

#' The 6 neighborings all give similar values of Moran's $I$, with slightly lower values as $K$ increases.  The dashed line corresponds to the range of $I$ for total population density, calculated above.  
#' 
#' Asian and black residents have moderate-high clustering.  White, Hispanic, and noncitizen residents have moderate clustering.  Children and impoverished residents seem to have clustering values the same as or just above the overall population.  Indigenous people have weak positive clustering.  


#' # Places #

places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    ## There are 9 places with total population 0
    filter(total_popE != 0) %>%
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
                                     poverty_combM, total_popM)) %>%
    ## Population densities
    mutate_at(vars(whiteE, blackE, indigenousE, asianE, hispanicE, 
                   noncitizensE, childrenE, poverty_combE), 
              funs(D = . / units::drop_units(area)))

glimpse(places_sf)

## Distributions of proportions ----
places_sf %>%
    as.data.frame() %>%
    select(ends_with('P')) %>%
    gather(key = variable, value) %>%
    ggplot(aes(value, color = variable, fill = variable)) +
    geom_density() +
    geom_rug() +
    facet_wrap(~ variable, scales = 'free') +
    scale_x_continuous(labels = scales::percent_format())

#' Slightly higher proportions across the board.  But no dramatic differences.  

places_sf %>%
    mutate(w_plus_h = whiteP + hispanicP) %>%
    ggplot(aes(w_plus_h)) + 
    geom_density() + 
    geom_rug()

#' Again, white+Hispanic supermajority

ggplot(places_sf, aes(hispanicP, poverty_combP)) + 
    geom_point() +
    geom_smooth(method = 'lm')

#' Again, greater Hispanic proportion is associated with a greater poverty rate.  


## Correlations ----
places_sf %>%
    as_tibble() %>%
    select(densityE, whiteP, blackP, childrenP, hispanicP, indigenousP, noncitizensP, poverty_combP, whiteP) %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'var1') %>%
    as_tibble() %>%
    gather(key = 'var2', value = 'cor', -var1) %>%
    mutate(cor.print = round(cor, digits = 1)) %>%
    ggplot(aes(var1, var2, fill = cor, label = cor.print)) + 
    geom_tile() +
    geom_text() +
    scale_fill_gradient2()

#' White is still anticorrelated with everything except Indigenous.  Strong or moderate correlations between noncitizens, poverty, or Hispanic. (Not children.)  Moderate correlations between Hispanic and density, and moderate-weak between Hispanic and children and density and noncitizens.  

ggplot(places_sf, aes((densityE), hispanicP)) + 
    geom_point() +
    geom_smooth()

#' Monotonic nonlinear relationship between density and Hispanic. 


## White/Hispanic segregation ----
## Evenness/dissimilarity
places_sf %>%
    as.data.frame() %>%
    mutate(w_h_dissim = abs(hispanicE/sum(hispanicE) - whiteE/sum(whiteE))) %>%
    summarize(w_h_dissim = .5 * sum(w_h_dissim))
#' 36% dissimilarity, lower than with tracts

ggplot(places_sf, aes(abs(hispanicP - whiteP))) + 
    geom_density() +
    geom_rug()

## Exposure/interaction
places_sf %>%
    as.data.frame() %>%
    mutate(w_h_exposure = abs(whiteE/sum(whiteE) * hispanicE / total_popE), 
           h_w_exposure = abs(hispanicE/sum(hispanicE) * whiteE/total_popE)) %>%
    summarize(w_h_exposure = sum(w_h_exposure), 
              h_w_exposure = sum(h_w_exposure))
#' Slightly higher White-Hispanic exposure, but still moderate-low

## Correlation
ggplot(places_sf, aes(hispanicP, whiteP)) + 
    geom_point()
with(tracts_sf, cor(hispanicP, whiteP))

#' Very strong negative correlation between the two proportions. But I guess, in this kind of two-group context, we would get a very strong negative correlation even if dissimilarity were low and exposure were high.  

places_sf %>%
    as.data.frame() %>%
    as_tibble() %>%
    summarize_at(vars(whiteE, hispanicE), funs(sum)) %>%
    mutate(hw_ratio = hispanicE / whiteE)
#' Hispanic-white ratio slightly higher, at 27%

## Spatial weights ----
## We'll construct KNN for K=3-8
library(spdep)

weights_places = 3:8 %>%
    set_names() %>%
    map(~ {places_sf %>%
            st_centroid() %>%
            st_coordinates() %>%
            knearneigh(k = .x) %>%
            knn2nb() %>%
            nb2listw(style = 'W')})

coords = places_sf %>%
    st_centroid() %>%
    st_coordinates()
plot(places_sf, max.plot = 1)
plot(weights_places$`8`, coords = coords, add = TRUE, col = 'blue')
plot(weights_places$`3`, coords = coords, add = TRUE, col = 'red')

## Both systems of neighbors produce an archipelago of tight clusters and longer connections.  Neither seems to produce ridiculously extended "neighbor" connections.  

## Moran's I ----
moran.i = function(vec, weights, ...) {
    moran.test(vec, 
               weights, ...)$estimate['Moran I statistic']
}

## Moran's I for overall density
weights_places %>%
    map(~moran.i(log10(places_sf$densityE), .)) %>%
    unlist()
#' Higher moderate population clustering, .45-.52

weights_places %>%
    tibble(weights = ., k = names(.)) %>%
    crossing(tibble(variable = c('whiteE_D', 
                                 'blackE_D', 
                                 'indigenousE_D', 
                                 'asianE_D', 
                                 'hispanicE_D', 
                                 'noncitizensE_D', 
                                 'childrenE_D', 
                                 'poverty_combE_D'))) %>%
    rowwise() %>%
    mutate(var_value = {places_sf %>% 
            as.data.frame() %>%
            pull(variable) %>%
            {. + 1} %>% log10() %>%
            list()}, 
           moran_i = moran.i(var_value, weights)) %>%
    select(k, variable, moran_i) %>%
    arrange(desc(moran_i)) %>%
    mutate(variable = fct_inorder(variable)) %>%
    ggplot(aes(variable, moran_i, color = k, group = k)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = c(.45, .52), linetype = 'dashed') +
    coord_flip()

#' With places, most groups have low and below-average clustering. Impoverished people, noncitizens, and Hispanics have moderate clustering, and Hispanics and noncitizens are above the overall average.  
