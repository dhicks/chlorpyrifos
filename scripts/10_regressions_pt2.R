## Spatial Durbin models
## NB On a standard-performance setup (eg, laptop) this script can take 10+ hours to run
## The parallelization below is configured for a high-performance compute server
## Make sure the configuration is appropriate for your setup
## Hard-coded bootstrap settings:  
## 500 block resamples for each primary model
## 100 simulations for the impacts on each resample model
library(tidyverse)
# library(sf)
# library(tidycensus)
library(spdep)
library(doSNOW)
library(tictoc)

## parallel::detectCores()
registerDoSNOW(makeCluster(25))

## Load data ----
data_dir = '../data/'
# data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

places_sfl = read_rds(str_c(data_dir, '07_places_sfl.Rds'))
tracts_sfl = read_rds(str_c(data_dir, '07_tracts_sfl.Rds'))

weights_pl = read_rds(str_c(data_dir, '07_places_weights.Rds'))
weights_tr = read_rds(str_c(data_dir, '07_tracts_weights.Rds'))

## Functions for fitting and resampling ----
perturb_ivs = function(data, formula) {
    ## This function uses the tidycensus MOEs to perturb the independent variables / draw samples from the measurement distributions for each variable
    
    ## Extract dependent variable from formula
    dv = formula %>%
        as.character() %>%
        .[2]
    
    ## Extract IVs
    vars = formula %>%
        as.character() %>%
        .[3] %>%
        str_split(fixed(' + ')) %>%
        unlist()
    # return(str_c(vars, collapse = '|'))
    
    ## Wide format, w/ IVs + their MOEs in columns
    data_wide = data %>%
        as.data.frame() %>%
        select(matches(str_c(vars, collapse = '|'))) %>%
        mutate(row_idx = row_number())
    ## Arrange the IVs in long format
    data_E = data_wide %>%
        select(-ends_with('M'), row_idx) %>%
        gather(key = 'variable', value = 'estimate', -row_idx) %>%
        mutate(prefix = str_extract(variable, '[a-z\\_]+'))
    ## Arrange the MOEs in long format
    data_M = data_wide %>%
        select(row_idx, ends_with('M')) %>%
        gather(key = 'variable', value = 'moe', -row_idx) %>%
        mutate(prefix = str_extract(variable, '[a-z\\_]+'),
               se = moe / qnorm(.95))
    
    ## Join them together, and draw from the measurement distribution for each variable
    data_long = inner_join(data_E, data_M, by = c('row_idx', 'prefix'))
    data_long = mutate(data_long, 
                       perturbed_value = rnorm(nrow(data_long), 
                                               mean = estimate, 
                                               sd = se))
    
    ## Arrange the draws back into a wide design matrix
    design_wide = data_long %>%
        select(row_idx, variable = variable.x, perturbed_value) %>%
        spread(key = variable, value = perturbed_value) %>%
        arrange(row_idx) %>%
        select(-row_idx)
    
    ## Recombine with the DV and any IVs that didn't have MOEs
    missed_vars = setdiff(vars, names(design_wide))
    design_wide = data %>%
        as.data.frame() %>%
        select(one_of(dv), one_of(missed_vars)) %>%
        bind_cols(design_wide)
    
    return(design_wide)
}

## fit_model() does the actual fitting
fit_model = function(data,
                     weights, # spatial weights
                     regression_formula,
                     zero.policy = NULL,
                     R.obs = 500, # impact draws
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
    impacts = impacts(model, tr = traces, R = R)
    
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


## This function constructs a certain number of resamples using perturb_ivs(), 
## then fits a SDM on each using fit_model().  
## If initialized, a parallel backend will be used for the fitting stage.  
resample_and_model = function(data, weights,
                              regression_formula, 
                              filter_condition = NULL, # string to filter data
                              k = 3, # n for knn; used in block resampling
                              R = 100, # num MCMC impact draws 
			      zero.policy = NULL,
                              perturb_ivs = FALSE, # Perturb independent values? 
                              n_resamples = 1, # Num. resample datasets
                              seed = NULL, # seed for RNG
                              return_model = FALSE # return the complete fitted models? 
) {
    ## Filter data using filter_condition
    if (!is.null(filter_condition)) {
        data = filter_(data, filter_condition)
    }
    
    ## Resampling
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (!perturb_ivs) {
        ## If we're not doing bootstrap, just use the unmodified data
        resample_datasets = list(data = data)
    } else {
        # resample_datasets = rep(list(data), n_resamples)
        resample_datasets = replicate(n_resamples, perturb_ivs(data, regression_formula), simplify = FALSE)
    }
    
    ## Fit model
    pb = txtProgressBar(max = n_resamples, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    fitted_models = foreach(x = resample_datasets, 
                            .packages = 'spdep', 
                            .export = 'fit_model',
                            .options.snow = list(progress = progress)
    ) %dopar% {
        fit_model(x, 
                  weights = weights, 
                  regression_formula = regression_formula, 
                  zero.policy = zero.policy, 
                  R = R,
                  return_model = return_model)
    }
    
    return(fitted_models)
}

## 8 sec for 1 resample of tract CTD -> 8 x 500 resamples x 5 CTD = 20k sec for tracts = ~5.6 hours
# tic()
# thing = resample_and_model(tracts_sfl[[1]], weights_tr, reg_form, n_resamples = 1, perturb_ivs = TRUE, R = 100)
# toc()
# str(thing, max.level = 2)

## Run resamples ----

reg_form = formula(log_w_use ~ hispanicP + noncitizensP +
                       blackP + indigenousP + asianP + womenP + 
                       childrenP + poverty_combP + 
                       ag_employedP + density_log10 + ag_employed_cP + density_log10_c)

print('Constructing model metadataframe')
models_meta_df = tibble(geography = c('places', 'tracts'), 
                        ctd = list(names(places_sfl), names(tracts_sfl)),
                        data = list(places_sfl, tracts_sfl), 
                        weights = list(weights_pl, weights_tr)) %>%
    unnest(ctd, data, .drop = FALSE)

write_rds(models_meta_df, str_c(data_dir, '10_models_meta.Rds'))

print('Fitting full spatial Durbin models')
tic()
durbin = foreach(row = iter(models_meta_df, by = 'row')
                 # .packages = c('tidyverse', 'sf', 'spdep'),
                 # .verbose = TRUE
) %do% {
    resample_and_model(row$data[[1]], row$weights[[1]],
                       regression_formula = reg_form,
                       seed = 78910,
                       zero.policy = TRUE,
		       R = 500,
                       perturb_ivs = FALSE)
}
toc()

write_rds(durbin, str_c(data_dir, '10_durbin_models.Rds'))

print('Fitting resamples')
tic()
resamples = foreach(row = iter(models_meta_df, by = 'row'),
                    .packages = c('tidyverse', 'spdep'),
                    .verbose = FALSE
) %do% {
    resample_and_model(row$data[[1]], row$weights[[1]],
                       regression_formula = reg_form,
                       seed = 1369,
                       zero.policy = TRUE, R = 100,
                       perturb_ivs = TRUE, n_resamples = 500)
}
toc()

write_rds(resamples, str_c(data_dir, '10_resamples.Rds'))
