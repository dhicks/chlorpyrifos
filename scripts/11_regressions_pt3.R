library(tidyverse)
library(coda)

## Load data ----
data_dir = '~/Google Drive/Coding/EJ datasets/CA pesticide/'

models_meta = read_rds(str_c(data_dir, '10_models_meta.Rds')) %>%
    select(-weights, -data) %>%
    mutate(CTD = fct_inorder(str_extract(ctd, '[0-9]+')), 
           model_idx = as.character(row_number()))
resamples = read_rds(str_c(data_dir, '10_resamples.Rds'))
durbin = read_rds(str_c(data_dir, '10_durbin_models.Rds'))

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
ggsave('11_rho.png', width = 6, height = 2)

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
ggsave('11_moran.png', width = 6, height = 2)

##
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
              'indigenousP' = 'Indigenous', 
              'density_log10' = 'pop. density (log)', 
              'poverty_combP' = 'poverty')

## Long-format impacts for ggplot
impacts_long = gather(impacts, 
                      key = 'variable', value = 'value', 
                      hispanicP:density_log10) %>%
    mutate(var_print = var_names[variable])
impacts_rs_long = gather(impacts_rs, 
                         key = 'variable', value = 'value', 
                         hispanicP:density_log10) %>%
    mutate(var_print = var_names[variable])

impacts_plot = impacts_rs_long %>%
    filter(!(ctd %in% c('ctd_1', 'ctd_10'))) %>%
    ggplot(aes(CTD, value, color = geography)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_violin(draw_quantiles = quants, 
                position = position_dodge(width = 1.25)) +
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
                 size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = .05),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = .95))
ggsave('11_impacts_369.png', width = 6, height = 4.5)

impacts_rs_long %>%
    filter(ctd == 'ctd_1') %>%
    {impacts_plot %+% .} +
    stat_summary(data = subset(impacts_long, 
                               ctd == 'ctd_1'),
                 size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = .05),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = .95))
ggsave('11_impacts_1.png', width = 6, height = 4.5)

impacts_rs_long %>%
    filter(ctd == 'ctd_10') %>%
    {impacts_plot %+% .} +
    stat_summary(data = subset(impacts_long, 
                               ctd == 'ctd_10'),
                 size = .25,
                 position = position_dodge(width = .25), 
                 fun.ymin = function (x) quantile(x, probs = .05),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = .95))
ggsave('11_impacts_10.png', width = 6, height = 4.5)

impacts_long %>%
    filter(ctd %in% c('ctd_90', 'ctd_30', 'ctd_60')) %>% 
    mutate(value.backtrans = 10^(value/10)) %>%
    ggplot(aes(CTD, value.backtrans, color = geography)) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    stat_summary(position = position_dodge(width = .75), 
                 fun.ymin = function (x) quantile(x, probs = .05),
                 fun.y = median, 
                 fun.ymax = function (x) quantile(x, probs = .95)) +
    facet_wrap(~ var_print, scales = 'free') +
    scale_y_continuous(name = 'Total impacts (trans.)' 
                       # labels = function (x) round(10^(x/10), 1)
                       # labels = scales::math_format(10^.x)
    ) +
    scale_color_brewer(palette = 'Set1', 
                       labels = c('places', 'tracts')) +
    theme_minimal()
ggsave('11_impacts_backtrans.png', width = 6, height = 4.5)

impacts_long %>%
    filter(ctd %in% c('ctd_60')) %>% 
    mutate(value.backtrans = 10^(value/10)) %>%
    group_by(geography, CTD, variable = var_print) %>%
    summarize(estimate = median(value), 
              ci_low = quantile(value, probs = .05), 
              ci_high = quantile(value, probs = .95)) %>%
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
    write_lines('11_impacts_table.txt')

