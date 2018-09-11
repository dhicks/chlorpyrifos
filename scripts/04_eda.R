#' This document conducts an EDA on census data only, first at the tract level and then at the places level.  
library(tidyverse)
library(sf)
library(spdep)

library(tidycensus)

data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

load_sf = function(rds_file) {
    rds_file %>%
        str_c(data_dir, .) %>%
        read_rds() %>%
        ## Remove locations w/ total population or total employed 0
        filter(total_popE != 0, total_employedE != 0) %>%
        ## Population proportions
        mutate(womenP = womenE / total_popE,
               womenPM = moe_prop(womenE, total_popE, womenM, total_popM),
               whiteP = whiteE / total_popE, 
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
                                         poverty_combM, total_popM), 
               hisp_povertyP = hisp_povertyE / hispanicE, 
               hisp_povertyPM = moe_prop(hisp_povertyE, hispanicE, hisp_povertyM, hispanicM), 
               ag_employedP = ag_employedE / total_employedE,
               ag_employedPM = moe_prop(ag_employedE, total_employedE, ag_employedM, total_employedM)
        ) %>%
        ## Population densities
        mutate_at(vars(womenE, whiteE, blackE, indigenousE, asianE, hispanicE, 
                       noncitizensE, childrenE, poverty_combE, 
                       hisp_povertyP, ag_employedP), 
                  funs(D = . / units::drop_units(area)))
}

#' # Tracts #
tracts_sf = load_sf('02_tracts_sf.Rds')
    
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

#' There are a few tracts with a modest proportion of Asian and Black residents ($> 20\%$); but only a few.  Almost no tracts have more than 5% Indigenous residents, and none have more than about 20%.  Children are typically ~5-12% of the population.  The poverty rate varies dramatically, with a median somewhere around 30% and some values above 50%.  Hispanic and White proportions are the most diverse.  Very little variation in proportion of women, though there are a few extreme tracks with values < 25% or > 80%

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
    select(densityE, ends_with('P')) %>%
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

tracts_sf %>%
    as_tibble() %>%
    select(densityE, ends_with('P'), -whiteP) %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'var1') %>%
    as_tibble() %>%
    gather(key = 'var2', value = 'cor', -var1) %>%
    filter(abs(cor) > .4, var1 < var2) %>%
    arrange(desc(abs(cor)))
    
#' White proportion has moderate to very strong negative corelations with every other variable (except Indigenous).  Very strong correlation between Hispanic and noncitizen proportion and between Hispanic and general poverty.  Strong correlations between agricultural employment and noncitizens and Hispanic.  Moderate correlations between each pair of poverty, children, noncitizens, and Hispanic, and between agricultural employment and poverty.  

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
coords_tracts = tracts_sf %>%
    st_centroid() %>%
    st_coordinates()
## Contiguity
weights_tracts_contig = tracts_sf %>%
    pull(geometry) %>%
    as_Spatial() %>%
    poly2nb() %>%
    nb2listw(style = 'W')
## KNN
weights_tracts_knn = 3:8 %>%
    set_names() %>%
    map(~ {knearneigh(coords_tracts, k = .x) %>%
            knn2nb() %>%
            nb2listw(style = 'W')})
## Inverse spatial weights w/in 50 km
dnn_tracts = dnearneigh(coords_tracts, d1 = 0, d2 = 50 * 1000)
weights_tracts_d = nbdists(dnn_tracts, coords = coords_tracts) %>%
    map( ~ 1/.) %>%
    nb2listw(dnn_tracts, glist = ., style = 'W', 
             zero.policy = TRUE)

weights_tracts = c(weights_tracts_knn, 
                   contiguity = list(weights_tracts_contig),
                   distance = list(weights_tracts_d))

plot(tracts_sf, max.plot = 1)
plot(weights_tracts$contiguity, coords = coords_tracts, 
     add = TRUE, col = 'blue')

## Moran's I ----
moran.i = function(vec, weights, ...) {
    moran.test(vec, 
               weights, ...)$estimate['Moran I statistic']
}

## Moran's I for overall density
moran_i_tracts = weights_tracts %>%
    map_dfr(~moran.i(log10(tracts_sf$densityE), .)) %>%
    gather(key = 'k', value = 'I')
moran_i_tracts
#' Moderate population clustering, ~.40-.45 for KNN

weights_tracts %>%
    tibble(weights = ., k = names(.)) %>%
    crossing(tibble(variable = c('womenE_D',
                                 'whiteE_D', 
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
    geom_hline(data = moran_i_tracts, 
               aes(yintercept = I, color = k), 
               linetype = 'dashed') +
    coord_flip()

#' The 6 KNN neighborings all give similar values of Moran's $I$, with slightly lower values as $K$ increases.  The dashed lines correspond to the values of $I$ for total population density, calculated above.  Distance-based weights have consistently lower values of Moran's $I$, but order the groups in basically the same way.  Continuity weights have consistently higher values of I, with almost no difference between different groups.  
#' 
#' Asian and black residents have moderate-high clustering.  White, Hispanic, and noncitizen residents have moderate clustering.  Children and impoverished residents seem to have clustering values the same as or just above the overall population.  Indigenous people have weak positive clustering.  
#' 


#' # Places #

places_sf = load_sf('02_places_sf.Rds')

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
    select(densityE, womenP, whiteP, blackP, childrenP, hispanicP, indigenousP, noncitizensP, poverty_combP, whiteP) %>%
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

places_sf %>%
    as_tibble() %>%
    select(densityE, ends_with('P'), -whiteP) %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'var1') %>%
    as_tibble() %>%
    gather(key = 'var2', value = 'cor', -var1) %>%
    filter(abs(cor) > .4, var1 < var2) %>%
    arrange(desc(abs(cor)))

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
library(spdep)

coords_places = places_sf %>%
    st_centroid() %>%
    st_coordinates()
## Contiguity
weights_places_contig = places_sf %>%
    pull(geometry) %>%
    as_Spatial() %>%
    poly2nb() %>%
    nb2listw(style = 'W', zero.policy = TRUE)
## KNN
weights_places_knn = 3:8 %>%
    set_names() %>%
    map(~ {knearneigh(coords_places, k = .x) %>%
            knn2nb() %>%
            nb2listw(style = 'W')})
## Inverse spatial weights w/in 50 km
dnn_places = dnearneigh(coords_places, d1 = 0, d2 = 50 * 1000)
weights_places_d = nbdists(dnn_places, coords = coords_places) %>%
    map( ~ 1/.) %>%
    nb2listw(dnn_places, glist = ., style = 'W', zero.policy = TRUE)

weights_places = c(weights_places_knn, 
                   contiguity = list(weights_places_contig),
                   distance = list(weights_places_d))

plot(places_sf, max.plot = 1)
plot(weights_places$contiguity, coords = coords_places, 
     add = TRUE, col = 'blue')

## All systems of neighbors produce an archipelago of tight clusters and longer connections.  Neither seems to produce ridiculously extended "neighbor" connections.  
## 
## Contiguity weights produce large numbers of islands:  246/397 (62%) have no neighbors

weights_places$contiguity$neighbours


## Moran's I ----
moran.i = function(vec, weights, ...) {
    moran.test(vec, 
               weights, ...)$estimate['Moran I statistic']
}

## Moran's I for overall density
moran_i_places = weights_places %>%
    map_dfr(~moran.i(log10(places_sf$densityE), ., 
                     zero.policy = TRUE)) %>%
    gather(key = 'k', value = 'I')
moran_i_places
#' Higher moderate population clustering, .45-.53. Distance weights are more consistent w/ KNN here.  Contiguity weights are much lower.  

weights_places %>%
    tibble(weights = ., k = names(.)) %>%
    crossing(tibble(variable = c('womenE_D',
                                 'whiteE_D', 
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
           moran_i = moran.i(var_value, weights, 
                             zero.policy = TRUE)) %>%
    select(k, variable, moran_i) %>%
    arrange(desc(moran_i)) %>%
    mutate(variable = fct_inorder(variable)) %>%
    ggplot(aes(variable, moran_i, color = k, group = k)) + 
    geom_point() +
    geom_line() +
    geom_hline(data = moran_i_places, 
               aes(yintercept = I, color = k), 
               linetype = 'dashed') +
    coord_flip()

#' With places, most groups have low and below-average clustering. Impoverished people, noncitizens, and Hispanics have moderate clustering, and Hispanics and noncitizens are above the overall average.  Distance values are generally similar to but a bit lower than the KNN.  Contiguity values are generally quite a bit lower.  
