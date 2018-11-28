library(tidyverse)
library(coda)

library(tictoc)
library(assertthat)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

models_meta = read_rds(str_c(data_dir, '10_models_meta.Rds')) %>%
    select(-weights, -data) %>%
    mutate(CTD = fct_inorder(str_extract(ctd, '[0-9]+')), 
           model_idx = as.character(row_number()))
resamples = read_rds(str_c(data_dir, '10_resamples.Rds'))
durbin = read_rds(str_c(data_dir, '10_durbin_models.Rds'))

## Impacts from the county-level models
impacts_co = read_rds(str_c(data_dir, '11_county_impacts_long.Rds'))


quants = c(.05, .5, .95)

## Rho ----
rho = durbin %>%
    modify_depth(2, 'rho') %>%
    modify_depth(1, bind_rows) %>%
    bind_rows(.id = 'model_idx') %>%
    inner_join(models_meta)
rho_rs = resamples %>%
    modify_depth(2, 'rho') %>%
    modify_depth(1, bind_rows, .id = 'resample') %>%
    bind_rows(.id = 'model_idx') %>%
    inner_join(models_meta)
ggplot(rho_rs, aes(CTD, rho, color = geography)) + 
    geom_violin(draw_quantiles = quants) +
    geom_point(data = rho, position = position_dodge(width = .9)) +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    ylab(expression(rho)) +
    theme_minimal()
ggsave('12_rho.png', width = 6, height = 2)

## Moran's I of residuals ----
moran_i = durbin %>%
    modify_depth(2, 'moran') %>%
    modify_depth(2, 'estimate') %>%
    modify_depth(2, 'Moran I statistic') %>%
    modify_depth(1, unlist) %>%
    map(~ tibble(moran_i = .)) %>%
    bind_rows(.id = 'model_idx') %>%
    inner_join(models_meta)
moran_i_rs = resamples %>%
    modify_depth(2, 'moran') %>%
    modify_depth(2, 'estimate') %>%
    modify_depth(2, 'Moran I statistic') %>%
    modify_depth(1, unlist) %>%
    map(~ tibble(moran_i = .)) %>%
    bind_rows(.id = 'model_idx') %>%
    group_by(model_idx) %>%
    mutate(resample = as.character(row_number())) %>%
    ungroup() %>%
    inner_join(models_meta)
ggplot(moran_i_rs, aes(CTD, moran_i, color = geography)) +
    geom_violin(draw_quantiles = quants) +
    geom_point(data = moran_i, position = position_dodge(width = .9)) +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    ylab("Moran's I of residuals") +
    theme_minimal()
ggsave('12_moran.png', width = 6, height = 2)

## Relationship btwn Moran's I and rho
inner_join(rho_rs, moran_i_rs) %>% 
    ggplot(aes(rho, moran_i, color = geography)) + 
    geom_point(alpha = .1) +
    geom_smooth(method = 'lm') +
    facet_wrap(~ ctd, scales = 'free') +
    theme_minimal()

## Impacts ----
impacts = durbin %>%
    modify_depth(2, 'impacts') %>%
    modify_depth(2, ~ .$sres) %>% 
    modify_depth(2, ~ .$total) %>% 
    modify_depth(2, as.data.frame) %>% 
    modify_depth(1, bind_rows, .id = 'resample') %>% 
    bind_rows(.id = 'model_idx') %>%
    inner_join(models_meta)
impacts_rs = resamples %>%
    modify_depth(2, 'impacts') %>%
    modify_depth(2, ~ .$sres) %>% 
    modify_depth(2, ~ .$total) %>% 
    modify_depth(2, as.data.frame) %>% 
    modify_depth(1, bind_rows, .id = 'resample') %>% 
    bind_rows(.id = 'model_idx') %>%
    inner_join(models_meta)

## Nicely-printing variable names
var_names = c('ag_employedP' = 'ag. employment', 
              'asianP' = 'Asian', 
              'blackP' = 'Black', 
              'childrenP' = 'children', 
              'hispanicP' = 'Hispanic', 
              'noncitizensP' = 'noncitizens',
              'indigenousP' = 'Indigenous', 
              'density_log10' = 'pop. density (log)', 
              'poverty_combP' = 'poverty', 
              'womenP' = 'women', 
              'ag_employed_cP' = 'county ag. employment', 
              'density_log10_c' = 'county pop. density (log)')

## Long-format impacts for ggplot
impacts_long = gather(impacts, 
                      key = 'variable', value = 'value', 
                      hispanicP:density_log10_c) %>%
    mutate(var_print = var_names[variable])
## ~20 sec
## Grouped version is way faster for some reason
tic()
impacts_rs_long = impacts_rs %>%
    group_by(model_idx, resample) %>%
    # sample_frac(1) %>%
    gather(key = 'variable', value = 'value', 
           hispanicP:density_log10_c) %>%
    mutate(var_print = var_names[variable]) %>%
    ungroup()
toc()
assert_that(are_equal(length(var_names)*nrow(impacts_rs), nrow(impacts_rs_long)))
impacts_co = impacts_co %>%
    rename(value = estimate) %>%
    mutate(var_print = var_names[variable], 
           CTD = {ctd %>% 
                   str_extract('[0-9]+') %>%
                   fct_inorder()})


impacts_plot = impacts_rs_long %>%
    filter(!(ctd %in% c('ctd_1', 'ctd_10'))) %>%
    ggplot(aes(CTD, value, color = geography)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # geom_violin(draw_quantiles = quants, 
    #             position = position_dodge(width = 1.25)) +
    stat_summary(fun.ymin = function (x) quantile(x, 
                                                  probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, 
                                                  probs = quants[3]), 
                 position = position_dodge(width = 1), 
                 aes(linetype = 'resamples')) +
    facet_wrap(~ var_print, scales = 'free') +
    scale_y_continuous(name = 'Total impacts' 
                       # labels = function (x) round(10^(x/10), 1)
                       # labels = scales::math_format(10^.x)
    ) +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    theme_minimal()
impacts_plot +
    stat_summary(data = subset(impacts_long, 
                               !(ctd %in% c('ctd_1', 'ctd_10'))),
                 # size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = quants[3]), 
                 aes(linetype = 'observed'))
ggsave('12_impacts_369.png', width = 6, height = 4.5)

impacts_rs_long %>%
    filter(ctd == 'ctd_1') %>%
    {impacts_plot %+% .} +
    stat_summary(data = subset(impacts_long, 
                               ctd == 'ctd_1'),
                 size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = quants[3]))
ggsave('12_impacts_1.png', width = 6, height = 4.5)

impacts_rs_long %>%
    filter(ctd == 'ctd_10') %>%
    {impacts_plot %+% .} +
    stat_summary(data = subset(impacts_long, 
                               ctd == 'ctd_10'),
                 size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = quants[3]))
ggsave('12_impacts_10.png', width = 6, height = 4.5)


## 30-60-90 w/ county-level estimates
impacts_long %>%
    filter(!(ctd %in% c('ctd_1', 'ctd_10')), 
           !str_detect(var_print, 'county')) %>%
    mutate(county = 'Full Data') %>%
    ggplot(aes(county, 
               value, color = geography)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # geom_violin(draw_quantiles = quants, 
    #             position = position_dodge(width = 1.25)) +
    stat_summary(position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = quants[3]), 
                 size = .1,
                 aes(linetype = 'full data', 
                     shape = 'full data')) +
    stat_summary(data = filter(impacts_co, !(ctd %in% c('ctd_1', 'ctd_10'))),
                 fun.ymin = function (x) quantile(x,
                                                  probs = quants[1]),
                 fun.y = median,
                 fun.ymax = function (x) quantile(x,
                                                  probs = quants[3]),
                 position = position_dodge(width = 1),
                 size = .1,
                 aes(linetype = 'county-level', 
                     shape = 'county-level')) +
    facet_wrap(~ CTD + var_print, scales = 'free_x', ncol = 5) +
    scale_x_discrete(limits = fct_inorder(
        c('Full Data', 
          rev(c('Butte', 'Fresno', 
          'Kern', 'San Joaquin', 
          'Solano', 'Stanislaus', 
          'Tulare'))))) +
    scale_y_continuous(name = 'Total impacts') +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    scale_shape(name = 'model') +
    scale_linetype(name = 'model') +
    coord_flip() +
    theme_minimal()
ggsave('12_impacts_co.png', width = 10, height = 10)


## Backtransformed impacts ----
impacts_long %>%
    filter(ctd %in% c('ctd_90', 'ctd_30', 'ctd_60')) %>% 
    mutate(value.backtrans = 10^(value/10)) %>%
    ggplot(aes(CTD, value.backtrans, color = geography)) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    stat_summary(position = position_dodge(width = .75), 
                 fun.ymin = function (x) quantile(x, probs = quants[1]),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = quants[3])) +
    facet_wrap(~ var_print, scales = 'free') +
    scale_y_continuous(name = 'Total impacts (fold)'
                       # labels = scales::percent_format()
    ) +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    theme_minimal()
ggsave('12_impacts_backtrans.png', width = 6, height = 4.5)

impacts_long %>%
    filter(ctd %in% c('ctd_60')) %>% 
    mutate(value.backtrans = 10^(value/10)) %>%
    group_by(geography, CTD, variable = var_print) %>%
    summarize(estimate = median(value), 
              ci_low = quantile(value, probs = quants[1]), 
              ci_high = quantile(value, probs = quants[3])) %>%
    ungroup() %>%
    mutate(estimate_trans = 10^(estimate/10), 
           geography = str_extract(geography, '[a-z]+')) %>%
    arrange(variable, geography) %>%
    select(geography, CTD, variable, estimate_trans, everything()) %>%
    knitr::kable(digits = 2, 
                 col.names = c('geography', 'CTD', 
                               'IV', 
                               'est. (transformed)', 
                               'estimate', '95% CI', '')) %>%
    write_lines('12_impacts_table.txt')

