#' """ Workshop 3: Ecosystem stability
#'     author: BSC6926-B53B
#'     date: 3/21/25"""

library(tidyverse)
library(hypervolume)

## Stability of seagrass ecosystems in Florida Bay
## data
# The data used for this example comes from the [South Florida Fisheries Habitat Assessment Program](https://myfwc.com/research/habitat/seagrasses/fhap/) which monitors seagrass habitats annually using quadrat samples. The [benthic cover data](https://github.com/SeascapeEcologyLab-workshops/BSC6926-B53B_Spring2025/blob/main/data/FLbay_SAV.csv) consists of data from 4 basins and measures 6 metrics of the SAV community. Data is averaged across stations for each metric. 
# 
# - BASIN = Basin sampled
# - YEAR = year of monitoring
# - STATION = monitoring station 
# - TT = *Thalassia testudium* percent cover
# - HW = *Halodule wrightii* percent cover
# - SF = *Syringodium filiforme* percent cover
# - TMA = total macroalgae percent cover
# - TDR = total drift algae percent cover
# - sg_rich = seagrass species richness\

# load sav monitoring data 
df = read_csv('data/FLbay_SAV.csv') 
head(df)


# ## Prepare data
# Because hypervolumes can be generated with any continuous data as an axes, many of the times the units are not combatible. Blonder et al. [2014](https://doi-org.ezproxy.fiu.edu/10.1111/geb.12146) & [2018](https://doi-org.ezproxy.fiu.edu/10.1111/2041-210X.12865) to convert all of the axes into the same units. This can be done by taking the z-score of the values to convert units into standard deviations. Z-scoring data can be done with the formula:
#   $$ z = \frac{x_{i}-\overline{x}}{sd} $$ Where $x_{i}$ is a value, $\overline{x}$ is the mean, and $sd$ is the standard deviation. By z-scoring each axis, 0 is the mean of that axis, a value of 1 means that the value is 1 standard deviation above the global mean of that axis, and a value of -1 is 1 standard deviation below the global mean of the axis. In R this can be done manually or with the `scale()` function. 
# 
# Hypervolumes cannot be made when all values for a single axis are the same (e.g. all values 0 for a species cover in a basin for a year), so we can add a tiny bit of variation in order to make the hypervolume.
# 
# We then can `nest()` the data to take advantage of the `purr` package and `map()`.

# z-score and nest data to make hypervolume
set.seed(14)
df = df |> 
  # z score data across all sites and years
  mutate(across(c(TT:sg_rich), scale), 
         # add tiny amount so when all values the same can make hv       
         across(c(TT:sg_rich), 
                ~map_dbl(., ~. + rnorm(1, mean = 0, sd = 0.0001)))) |> 
  # remove station from dataset
  select(-STATION) |> 
  # nest data by basin and year
  group_by(BASIN, YEAR) |> 
  nest() 
head(df)


## Generate hypervolumes
# Hypervolumes are a multidimensional tool that is based on Hutchinson's *n*-dimensional niche concept and we can build them with the `hypervolume` package. \
# 
# With a nested dataset of our columns that we want to build hypervolumes for we can use `mutate()` and `map()` to generate the hypervolume. 
# 
# We can also use `map()` and `get_centroid()` to extract centroid values of each hypervolume. 

# generate hypervolumes
df = df |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(BASIN,YEAR,sep = '_'),
                                                     samples.per.point = 100,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         centroid = map(hv, \(hv) get_centroid(hv)),
         size = map_dbl(hv, \(hv) get_volume(hv)))

saveRDS(df, 'data/SAV_hvs.rds') 


head(df)

# If wanting to save you can save output as `.rds`

df = readRDS('data/SAV_hvs.rds')

# Comparison metrics
# We can use the overlap to understand the similarity between years, centroid distance to compare mean conditions between years, and the log of the size ratio between years to understand the stability. This can be done by creating a data frame with all of the possible year combinations, and merging dataframes together to easily join. We can then use `map()` to calculate the metrics between.

# comparison of across each year
df_y= tibble(y1 = unique(df$YEAR),
             y2 = unique(df$YEAR)) |> 
  expand(y1,y2)

# make all unique year comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df1 = df |> 
  select(BASIN, y1 = YEAR, hv1 = hv, hv1_size = size, cent1 = centroid)

df2 = df |> 
  select(BASIN, y2 = YEAR, hv2 = hv, hv2_size = size, cent2 = centroid)


# create data frame of all data and make yearly comparisons
df_ov = tibble(BASIN = rep(unique(df$BASIN),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
               y2 = rep(df_y$y2, times = length(unique(df$BASIN)))) |> 
  inner_join(df1, by = c('BASIN', 'y1')) |> 
  inner_join(df2, by = c('BASIN', 'y2')) |> 
  mutate(ychange = y2-y1,
  # calculate the differnces in size 
         lsr = log(hv2_size/hv1_size),
  # join hypervolumees in a set for cverlap
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
  # calculate overlap
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
  # calculate centroid distance 
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  #unnest centroid differences
  unnest_wider(ov) |> 
  # select only metrics of interest
  select(BASIN, y1, y2, ychange,lsr,
         dist_cent, jaccard, sorensen,
         uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2)

# save output
write_csv(df_ov, 'data/SAV_centDist.csv')

df_ov = read_csv('data/SAV_centDist.csv')

df_ov

# Plot Overlap, centroid distance, and log size ratio of all comparisons.
df_ov = read_csv('data/SAV_centDist.csv') |> 
mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP','BLK')))
#overlap
ggplot(df_ov, aes(BASIN, sorensen, fill = BASIN))+
  geom_point(aes(color = BASIN), size = 1, 
             position=position_jitterdodge(dodge.width = 0.75, jitter.width = 1))+
  # geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  geom_boxplot(alpha = 0.6, outliers = F)+
  labs(x = 'Basin', y = 'Overlap')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# centroid distance
ggplot(df_ov, aes(BASIN, dist_cent, fill = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(aes(color = BASIN), size = 1, 
             position=position_jitterdodge(dodge.width = 0.75, jitter.width = 1))+
  # geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  geom_boxplot(alpha = 0.6, outliers = F)+
  labs(x = 'Basin', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


# log size ratio
ggplot(df_ov, aes(BASIN, lsr, fill = BASIN))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 1)+
  geom_point(aes(color = BASIN), size = 1, 
             position=position_jitterdodge(dodge.width = 0.75, jitter.width = 1))+
  # geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  geom_boxplot(alpha = 0.6, outliers = F)+
  labs(x = 'Basin', y = 'log(y2 size/y1 size)')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

## Year to Year comparisons
# A common way to look at stability is to plot change over time by plotting the year to year comparisons of each metric. 

# filter only single year comparisons
d = df_ov |> 
  filter(ychange == 1)

# overlap
ggplot(d, aes(y2, sorensen, color = BASIN))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Overlap')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# centroid distance 
ggplot(d, aes(y2, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# log size ratio
ggplot(d, aes(y2, lsr, color = BASIN))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'log(y2 size/y1 size)')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


## Compare to baseline
# Another way to look at stability is to plot change relative to a baseline. Here we plot them relative to the first year of the dataset (2007).

# filter only single year comparisons
d = df_ov |> 
  filter(y1 == 2007)

# overlap
ggplot(d, aes(y2, sorensen, color = BASIN))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Overlap')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# centroid distance 
ggplot(d, aes(y2, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


# log size ratio
ggplot(d, aes(y2, lsr, color = BASIN))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'log(y2 size/y1 size)')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

## Trend in stability 
# We can look across the number of years to understand the trend in stability. By fitting three possible models we can determine the trend of time For each metric. When intercept model is the best, we can determine that the trend is static and not changing with the number of years between comparison. If linear, the centroid distance can indicate a shift in state overtime, and a quadratic with a maximum at middle values can indicate a disturbance with recovery in state. 
library(MuMIn)
# overlap 
df_o = df_ov |> 
  group_by(BASIN) |>
  nest() |> 
  # fit intercept, linear, and quadratic model
  mutate(m_int = map(data, \(df)lm(sorensen~1, data = df)),
         m_lin = map(data, \(df)lm(sorensen~ychange, data = df)),
         m_quad = map(data, \(df)lm(sorensen~ychange + I(ychange^2), data = df)),
         AICc_int = map_dbl(m_int, \(x) AICc(x)),
         AICc_lin = map_dbl(m_lin, \(x) AICc(x)),
         AICc_quad = map_dbl(m_quad, \(x) AICc(x)),
         model = case_when(
           AICc_int - min(c(AICc_int,AICc_lin,AICc_quad)) <= 4 ~ 'Intercept',
           AICc_lin - AICc_quad <= 4 ~ 'Linear',
           T ~ 'Quadratic'))

# unnest data 
d = df_o |> 
  select(BASIN, data, model) |> 
  unnest(cols = c(data)) |>  
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP','BLK')))

ggplot(d, aes(ychange, sorensen, color = BASIN))+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~BASIN,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'Overlap')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# centroid distance
df_cd = df_ov |> 
  group_by(BASIN) |>
  nest() |> 
  # fit intercept, linear, and quadratic model
  mutate(m_int = map(data, \(df)lm(dist_cent~1, data = df)),
         m_lin = map(data, \(df)lm(dist_cent~ychange, data = df)),
         m_quad = map(data, \(df)lm(dist_cent~ychange + I(ychange^2), data = df)),
         AICc_int = map_dbl(m_int, \(x) AICc(x)),
         AICc_lin = map_dbl(m_lin, \(x) AICc(x)),
         AICc_quad = map_dbl(m_quad, \(x) AICc(x)),
         model = case_when(
           AICc_int - min(c(AICc_int,AICc_lin,AICc_quad)) <= 4 ~ 'Intercept',
           AICc_lin - AICc_quad <= 4 ~ 'Linear',
           T ~ 'Quadratic'))

# unnest data 
d = df_cd |> 
  select(BASIN, data, model) |> 
  unnest(cols = c(data)) |>  
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP','BLK')))

ggplot(d, aes(ychange, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~BASIN,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'Centroid distance')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# log size ratio
df_lsr = df_ov |> 
  group_by(BASIN) |>
  nest() |> 
  # fit intercept, linear, and quadratic model
  mutate(m_int = map(data, \(df)lm(lsr~1, data = df)),
         m_lin = map(data, \(df)lm(lsr~ychange, data = df)),
         m_quad = map(data, \(df)lm(lsr~ychange + I(ychange^2), data = df)),
         AICc_int = map_dbl(m_int, \(x) AICc(x)),
         AICc_lin = map_dbl(m_lin, \(x) AICc(x)),
         AICc_quad = map_dbl(m_quad, \(x) AICc(x)),
         model = case_when(
           AICc_int - min(c(AICc_int,AICc_lin,AICc_quad)) <= 4 ~ 'Intercept',
           AICc_lin - AICc_quad <= 4 ~ 'Linear',
           T ~ 'Quadratic'))

# unnest data 
d = df_lsr |> 
  select(BASIN, data, model) |> 
  unnest(cols = c(data)) |>  
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP','BLK')))

ggplot(d, aes(ychange, lsr, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~BASIN,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'log(y2 size/y1 size)')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))



## Variable importance 
# Because hypervolumes are multivariate, each axis has the potential to influence the overall change. We can determine the influence of each axis on the overall change by removing an axis and recalculating each metric. We can then compare the hypervolume without that axis to the metrics with all axes. We can determine the importance of a variable using the following formula:
# $$ imp_x = 1 - r_i$$
# where $imp_x$ is the importance of axis $x$, and $r_i$ is pearson correlation coefficient between the metric $i$ (overlap, centroid distance, or log size ratio) of the hypervolume with all axes to the hypervolume calculate without axis $x$. High importance values indicate that removing the axis had big changes on the values of that metric. 

# get data used for hvs
df = readRDS('data/SAV_hvs.rds') |> 
  select(BASIN, YEAR, data) |> 
  unnest(data)

#axes 
ax = c("TT", "HW", "SF", "TMA", "TDR", "sg_rich")

for (i in 1:length(ax)){
  d = df |> 
    select(-ax[i])|> 
    group_by(BASIN, YEAR) |>
    nest() |> 
    mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(BASIN,YEAR,sep = '_'),
                                                       samples.per.point = 1000,
                                                       kde.bandwidth = estimate_bandwidth(data), 
                                                       sd.count = 3, 
                                                       quantile.requested = 0.95, 
                                                       quantile.requested.type = "probability", 
                                                       chunk.size = 1000, 
                                                       verbose = F)),
           hv_size = map_dbl(hv, \(hv) get_volume(hv)),
           axis = ax[i])
  
  if (i == 1){
    df_tot = d
  }else{
    df_tot = bind_rows(df_tot, d)
  }
}

# make all comparisons for each axis
# comparison of across each year
df_y= tibble(y1 = unique(df$YEAR),
             y2 = unique(df$YEAR)) |> 
  expand(y1,y2)

# make all unique year comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df1 = df_tot |> 
  select(BASIN, axis, y1 = YEAR, hv1 = hv, hv1_size = hv_size)

df2 = df_tot |> 
  select(BASIN, axis, y2 = YEAR, hv2 = hv, hv2_size = hv_size)


# create data frame of all data and make yearly comparisons
df_ax = tibble(BASIN = rep(unique(df$BASIN),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
               y2 = rep(df_y$y2, times = length(unique(df$BASIN)))) |>
  slice(rep(1:n(), each=length(ax)))|> 
  mutate(axis = rep(ax, times=nrow(df_y)*length(unique(df$BASIN)))) |> 
  inner_join(df1, by = c('BASIN', 'y1', 'axis')) |> 
  inner_join(df2, by = c('BASIN', 'y2', 'axis')) |> 
  mutate(ychange = y2-y1,
  # calculate the differences in size 
         lsr_ax = log(hv2_size/hv1_size),
  # join hypervolumes in a set for overlap
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
  # calculate overlap
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
  # calculate centroid distance 
         dist_cent_ax = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  #unnest overlap
  unnest_wider(ov) |> 
  # select only metrics of interest
  select(BASIN, y1, y2, ychange, axis,
         lsr_ax, dist_cent_ax, sorensen_ax = sorensen)

# save output
write_csv(df_ax, 'data/SAV_metrics_ax.csv')

# Calculate correlation and plot
# load and combine with whole data
df_ov = read_csv('data/SAV_centDist.csv') |> 
  select(BASIN, y1, y2, lsr, dist_cent, sorensen)

df_all = read_csv('data/SAV_metrics_ax.csv')|> 
  left_join(df_ov) |> 
  group_by(BASIN, axis) |> 
  nest() |> 
  mutate(cor_ov = map(data, \(df) cor.test(df$sorensen, df$sorensen_ax)),
         imp_ov = map_dbl(cor_ov, \(x) 1 - x$estimate),
         cor_cd = map(data, \(df) cor.test(df$dist_cent, df$dist_cent_ax)),
         imp_cd = map_dbl(cor_cd, \(x) 1 - x$estimate),
         cor_lsr = map(data, \(df) cor.test(df$lsr, df$lsr_ax)),
         imp_lsr = map_dbl(cor_lsr, \(x) 1 - x$estimate),
         BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP','BLK')),
         axis = factor(axis, levels = c('TDR', 'TMA', 'sg_rich',
                                        'SF', 'HW', 'TT')))

# overlap
ggplot(df_all, aes(axis, imp_ov, fill = BASIN))+
  geom_col()+
  labs(x = 'Variable', y = 'Overlap variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~BASIN,  nrow = 2)+
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# centroid distance
ggplot(df_all, aes(axis, imp_cd, fill = BASIN))+
  geom_col()+
  labs(x = 'Variable', y = 'Centroid distance variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~BASIN,  nrow = 2)+
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# log size ratio
ggplot(df_all, aes(axis, imp_lsr, fill = BASIN))+
  geom_col()+
  labs(x = 'Variable', y = 'log(y2 size/y1 size) variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~BASIN,  nrow = 2)+
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


plot(dfs_hvs$hv[[1]])
plot(dfs_hvs$hv[[2]])
plot(dfs_hvs$hv[[3]])


hvj = hypervolume_join(dfs_hvs$hv[[1]], dfs_hvs$hv[[2]], dfs_hvs$hv[[3]])
plot(hvj)

# comparison of across each year
df_y= tibble(y1 = unique(df$YEAR),
             y2 = unique(df$YEAR)) |> 
  expand(y1,y2)

# make all unique year comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df1 = df |> 
  select(BASIN, y1 = YEAR, hv1 = hv, hv1_size = size, cent1 = centroid)

df2 = df |> 
  select(BASIN, y2 = YEAR, hv2 = hv, hv2_size = size, cent2 = centroid)


# create data frame of all data and make yearly comparisons
df_ov = tibble(BASIN = rep(unique(df$BASIN),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
               y2 = rep(df_y$y2, times = length(unique(df$BASIN)))) |> 
  inner_join(df1, by = c('BASIN', 'y1')) |> 
  inner_join(df2, by = c('BASIN', 'y2')) |> 
  mutate(ychange = y2-y1,
  # calculate the differnces in size 
         lsr = log(hv2_size/hv1_size),
  # join hypervolumees in a set for cverlap
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
  # calculate overlap
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
  # calculate centroid distance 
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  #unnest centroid differences
  unnest_wider(ov) |> 
  # select only metrics of interest
  select(BASIN, y1, y2, ychange,lsr,
         dist_cent, jaccard, sorensen,
         uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2)

# save output
write_csv(df_ov, 'data/SAV_centDist.csv')
