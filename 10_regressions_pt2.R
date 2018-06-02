## Spatial Durbin models
library(tidyverse)
library(sf)
library(tidycensus)
library(spdep)
library(doSNOW)
library(tictoc)

registerDoSNOW(makeCluster(2))

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sf = read_rds(str_c(data_dir, '02_places_sf.Rds')) %>%
    mutate(places_idx = row_number()) %>%
    ## Remove places w/ total population or total employed 0
    filter(total_popE != 0, total_employedE != 0) %>%
    ## Log density
    mutate(log_densityE = log10(densityE)) %>%
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
                                     poverty_combM, total_popM), 
           hisp_povertyP = hisp_povertyE / hispanicE, 
           hisp_povertyPM = moe_prop(hisp_povertyE, hispanicE, hisp_povertyM, hispanicM), 
           ag_employedP = ag_employedE / total_employedE,
           ag_employedPM = moe_prop(ag_employedE, total_employedE, ag_employedM, total_employedM)
    ) %>%
    ## Population densities
    mutate_at(vars(whiteE, blackE, indigenousE, asianE, hispanicE, 
                   noncitizensE, childrenE, poverty_combE, 
                   hisp_povertyP, ag_employedP), 
              funs(D = . / units::drop_units(area))) %>%
    ## Bind w/ locally-weighted total use
    left_join(read_rds(str_c(data_dir, '06_w_use_places.Rds'))) %>%
    ## Logged use
    mutate(log_w_use = log10(w_use))

tracts_sf = read_rds(str_c(data_dir, '02_tracts_sf.Rds')) %>%
    mutate(tracts_idx = row_number()) %>%
    ## Remove tracts w/ total population or total employed 0
    filter(total_popE != 0, total_employedE != 0) %>%
    ## Log density
    mutate(log_densityE = log10(densityE)) %>%
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
                                     poverty_combM, total_popM), 
           hisp_povertyP = hisp_povertyE / hispanicE, 
           hisp_povertyPM = moe_prop(hisp_povertyE, hispanicE, hisp_povertyM, hispanicM), 
           ag_employedP = ag_employedE / total_employedE,
           ag_employedPM = moe_prop(ag_employedE, total_employedE, ag_employedM, total_employedM)
    ) %>%
    ## Population densities
    mutate_at(vars(whiteE, blackE, indigenousE, asianE, hispanicE, 
                   noncitizensE, childrenE, poverty_combE, 
                   hisp_povertyP, ag_employedP), 
              funs(D = . / units::drop_units(area))) %>%
    ## Bind w/ locally-weighted total use
    left_join(read_rds(str_c(data_dir, '06_w_use_tracts.Rds'))) %>%
    ## Logged use
    mutate(log_w_use = log10(w_use))


## Functions for fitting and resampling ----
perturb_ivs = function(data, formula) {
    dv = formula %>%
        as.character() %>%
        .[2]
    
    vars = formula %>%
        as.character() %>%
        .[3] %>%
        str_split(fixed(' + ')) %>%
        unlist()
    # return(str_c(vars, collapse = '|'))
    
    data_wide = data %>%
        as.data.frame() %>%
        select(matches(str_c(vars, collapse = '|'))) %>%
        mutate(row_idx = row_number())
    data_E = data_wide %>%
        select(-ends_with('M'), row_idx) %>%
        gather(key = 'variable', value = 'estimate', -row_idx) %>%
        mutate(prefix = str_extract(variable, '[a-z]+'))
    data_M = data_wide %>%
        select(row_idx, ends_with('M')) %>%
        gather(key = 'variable', value = 'moe', -row_idx) %>%
        mutate(prefix = str_extract(variable, '[a-z]+'),
               se = moe / qnorm(.95))
    
    data_long = inner_join(data_E, data_M, by = c('row_idx', 'prefix'))
    
    data_long = mutate(data_long, 
                       perturbed_value = rnorm(nrow(data_long), 
                                               mean = estimate, 
                                               sd = se))
    design_wide = data_long %>%
        select(row_idx, variable = variable.x, perturbed_value) %>%
        spread(key = variable, value = perturbed_value) %>%
        arrange(row_idx) %>%
        select(-row_idx)
    missed_vars = setdiff(vars, names(design_wide))
    design_wide = data %>%
        as.data.frame() %>%
        select(one_of(dv), one_of(missed_vars)) %>%
        bind_cols(design_wide)
    
    return(design_wide)
}

{
    ## Probably for annoying scoping reasons, the parallel setup in resample_and_model() can't find fit_model() unless it's defined inside the former.  
    # fit_model = function(data,
    #                      weights, # spatial weights
    #                      regression_formula,
    #                      zero.policy = NULL,
    #                      return_model = FALSE # return the complete fitted model?
    # ) {
    #     ## Calculate traces
    #     traces = trW(as(weights, 'CsparseMatrix'))
    # 
    #     ## Fit model
    #     model = lagsarlm(regression_formula,
    #                      data = data,
    #                      listw = weights,
    #                      trs = traces,
    #                      type = 'Durbin',
    #                      zero.policy = zero.policy)
    #     ## Calculate impacts
    #     impacts = impacts(model, tr = traces, R = 100)
    # 
    #     ## Moran's I of residuals
    #     # moran = moran.mc(residuals(model),
    #     #                  weights, nsim = 500,
    #     #                  zero.policy = zero.policy)
    #     moran = moran.test(residuals(model), weights,
    #                        zero.policy = zero.policy)
    # 
    #     ## Return results
    #     results = list(rho = data.frame(rho = model$rho,
    #                                     rho.se = model$rho.se),
    #                    aic = AIC(model),
    #                    impacts = impacts,
    #                    moran = moran,
    #                    formula = regression_formula)
    #     if (return_model) {
    #         results$model = model
    #     }
    #     return(results)
    # }
}

resample_and_model = function(data, 
                              regression_formula, 
                              filter_condition = NULL, # string to filter data
                              k = 3, # for KNN weights
                              zero.policy = NULL,
                              do_bootstrap = FALSE, # Do resamples? 
                              n_resamples = 1, # Num. resample datasets
                              seed = NULL, # seed for RNG
                              return_model = FALSE # return the complete fitted models? 
) {
    fit_model = function(data,
                         weights, # spatial weights
                         regression_formula,
                         zero.policy = NULL,
                         return_model = FALSE # return the complete fitted model?
    ) {
        ## Calculate traces
        traces = trW(as(weights, 'CsparseMatrix'))
        
        ## Fit model
        model = lagsarlm(regression_formula,
                         data = data,
                         listw = weights,
                         trs = traces,
                         type = 'Durbin',
                         zero.policy = zero.policy)
        ## Calculate impacts
        impacts = impacts(model, tr = traces, R = 100)
        
        ## Moran's I of residuals
        # moran = moran.mc(residuals(model),
        #                  weights, nsim = 500,
        #                  zero.policy = zero.policy)
        moran = moran.test(residuals(model), weights,
                           zero.policy = zero.policy)
        
        ## Return results
        results = list(rho = data.frame(rho = model$rho,
                                        rho.se = model$rho.se),
                       aic = AIC(model),
                       impacts = impacts,
                       moran = moran,
                       formula = regression_formula)
        if (return_model) {
            results$model = model
        }
        return(results)
    }
    
    ## Filter data using filter_condition
    if (!is.null(filter_condition)) {
        data = filter_(data, filter_condition)
    }
    
    ## Construct spatial weights
    weights_knn = data %>%
        st_centroid %>%
        st_coordinates() %>%
        knearneigh(k = k) %>%
        knn2nb() %>%
        nb2listw(style = 'W')
    
    ## Resampling
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (!do_bootstrap) {
        ## If we're not doing bootstrap, just use the unmodified data
        resample_datasets = list(list(data = data, weights = weights_knn))
    } else {
        ## Sample blocks
        n_blocks = floor(nrow(data) / (k+1))
        ## Blocks are defined by their "centers"
        centers = sample(1:nrow(data), n_resamples*n_blocks, replace = TRUE)
        
        ## Go from index of "centers" to index of all locations
        resample_locations = weights_knn$neighbours[centers] %>%
            ## Concatenate the centers, each center in a singleton list
            c({centers %>% list() %>% transpose() %>% flatten()}) %>%
            tibble(location = .) %>%
            mutate(resample = rep_len(1:n_resamples, 
                                      length.out = 2*n_resamples*n_blocks)) %>%
            unnest() %>%
            split(.$resample) %>%
            map(~.$location)
        
        ## Subset the spatial weights matrix and coerce back to listw
        weights_matrix = as(weights_knn, 'CsparseMatrix')
        resample_weights = map(resample_locations, ~ weights_matrix[., .]) %>% 
            ## Row standardize
            map(~ ./rowSums(.)[row(.)]) %>% 
            map(function (x) {x[is.na(x)] = 0; return(x)}) %>% 
            map(quietly(mat2listw), style = 'M') %>%
            transpose() %>%
            .$result %>%
            ## Pass a lagsarlm() check for row standardization
            map(function (x) {x$style = 'W'; return(x)})
        
        ## Subset data
        resample_data = map(resample_locations, ~ data[.,]) %>%
            map(perturb_ivs, regression_formula)
        
        ## Put everything together
        resample_datasets = list(data = resample_data, 
                                 weights = resample_weights) %>%
            transpose()
    }
    
    ## Fit model
    fitted_models = resample_datasets %>%
        ## Something in the interaction w/ doSNOW causes warnings
        ## cf <https://github.com/hadley/plyr/issues/203>
        plyr::llply(function (x) 
            fit_model(x$data, 
                      weights = x$weights, 
                      regression_formula,
                      zero.policy = zero.policy, 
                      return_model = return_model), 
            .parallel = TRUE, 
            .paropts = list(#.export = 'regression_formula', 
                .packages = c('spdep'))) %>%
        map(function (x) {x$k = k; return(x)})
    
    return(fitted_models)
}


## Run resamples ----

reg_form = formula(log_w_use ~ hispanicP + blackP + indigenousP + asianP + childrenP + poverty_combP + ag_employedP + log_densityE)
# tictoc::tic()
# returned = resample_and_model(tracts_sf,
#                               reg_form,
#                               filter_condition = 'ctd == "ctd_60"', 
#                               k = 3, 
#                               zero.policy = TRUE, 
#                               return_model = FALSE,
#                               do_bootstrap = TRUE, 
#                               n_resamples = 2)
# tictoc::toc()

models_meta_df = cross_df(list(dataset = c('places_sf', 'tracts_sf'), 
                               k = as.integer(seq(3, 9, by = 2)), 
                               ctd = unique(places_sf$ctd)
)) %>%
    mutate(model_idx = as.character(row_number()))

tic()
resamples = foreach(row = iter(models_meta_df, by = 'row'), 
                .packages = c('tidyverse', 'sf', 'spdep'), 
                .export = c('places_sf', 'tracts_sf')
                ) %do% {
    resample_and_model(eval(parse(text = row$dataset)), 
                       reg_form,
                       filter_condition = str_c('ctd == \"', 
                                                row$ctd, '\"'), 
                       k = row$k, 
                       zero.policy = TRUE, 
                       do_bootstrap = TRUE, n_resamples = 10)
                }
toc()

write_rds(models_meta_df, str_c(data_dir, '10_models_meta.Rds'))
write_rds(resamples, str_c(data_dir, '10_resamples.Rds'))
