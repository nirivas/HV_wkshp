

# Load Libraries ----------------------------------------------------------

library(tidyverse)

library(hypervolume)



# Load Data ---------------------------------------------------------------

df = read.csv("data/LDWFseinedata.csv") |> 
  filter(basin == "Calcasieu") |> 
  rename_with(tolower) |> 
  mutate(across(where(is.character), tolower),
         year = substr(date, 1,4),
         month = substr(date, 5,6),
         day = substr(date, 7,8)) |>  
  rename(turbidity = tubidity) |> 
  select(-date)



# Create HVs --------------------------------------------------------------

set.seed(1)

## Scaling ####
df2 <- df |> 
  group_by(lat, lon, year, month, day) |> 
  slice(1) |> 
  ungroup() |>
  mutate(across(salinity:airtemp, scale)) |> 
  select(year,salinity, watertemp, turbidity,airtemp) |> 
  group_by(year) |> 
  nest()

head(df2)

## HVs Creation ####

df2 = df2 |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(year),
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         centroid = map(hv, \(hv) get_centroid(hv)),
         size = map_dbl(hv, \(hv) get_volume(hv)))

beepr::beep(3)

head(df2)

#saveRDS(df2, 'data/first_hvs.rds') 


# Yearly Comparisons ------------------------------------------------------

## Load Data ####
df2 = readRDS('data/first_hvs.rds') |> 
  mutate(year=as.numeric(year))


# comparison of across each year
df_y= tibble(y1 = unique(df2$year),
             y2 = unique(df2$year)) |> 
  expand(y1,y2)

# make all unique year comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df3 = df2 |> 
  select( y1 = year, hv1 = hv, hv1_size = size, cent1 = centroid)

df4 = df2 |> 
  select( y2 = year, hv2 = hv, hv2_size = size, cent2 = centroid)

# create data frame of all data and make yearly comparisons
df_ov = tibble(year = rep(unique(df2$year),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df2$year))),
               y2 = rep(df_y$y2, times = length(unique(df2$year)))) |> 
  inner_join(df3, by = c('y1')) |> 
  inner_join(df4, by = c('y2')) |> 
  mutate(ychange = y2-y1,
         # calculate the differnces in size 
         lsr = log(hv2_size/hv1_size),
         # join hypervolumees in a set for cverlap
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         # calculate overlap
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         # calculate centroid distance 
         dist_cent = 
           map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  #unnest centroid differences
  unnest_wider(ov) |> 
  # select only metrics of interest
  select(year, y1, y2, ychange,lsr,
         dist_cent, jaccard, sorensen,
         uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2)

beepr::beep(2)


# save output
#write_csv(df_ov, 'data/first_centDist.csv')



## Plotting ####



### Year to Year comparisons ####
# A common way to look at stability is to plot change over time by plotting the year to year comparisons of each metric. 

df_ov=read_csv("data/first_centDist.csv") 

# filter only single year comparisons
d = df_ov |> 
  filter(ychange == 1)

# overlap
ggplot(d, aes(y2, sorensen,color="blue"))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Overlap')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  #facet_wrap(~BASIN, nrow = 2)+
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
ggplot(d, aes(y2, dist_cent))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  #facet_wrap(~BASIN, nrow = 2)+
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
ggplot(d, aes(y2, lsr))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'log(y2 size/y1 size)')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  #facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

### Compare to baseline ####
# Another way to look at stability is to plot change relative to a baseline. Here we plot them relative to the first year of the dataset (1986).

# filter only single year comparisons
d = df_ov |> 
  filter(y1 == 1986)

# overlap
ggplot(d, aes(y2, sorensen))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Overlap')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  #facet_wrap(~BASIN, nrow = 2)+
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
ggplot(d, aes(y2, dist_cent))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
 # facet_wrap(~BASIN, nrow = 2)+
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
ggplot(d, aes(y2, lsr))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 1)+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'log(y2 size/y1 size)')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  #facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))



### Trend in stability ####
# We can look across the number of years to understand the trend in stability. By fitting three possible models we can determine the trend of time For each metric. When intercept model is the best, we can determine that the trend is static and not changing with the number of years between comparison. If linear, the centroid distance can indicate a shift in state overtime, and a quadratic with a maximum at middle values can indicate a disturbance with recovery in state. 

library(MuMIn)
# overlap 
df_o = df_ov |> 
 # group_by(basin) |>
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
  select(data, model) |> 
  unnest(cols = c(data)) 

ggplot(d, aes(ychange, sorensen))+
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
  select(data, model) |> 
  unnest(cols = c(data))


ggplot(d, aes(ychange, dist_cent))+
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
  select(data, model) |> 
  unnest(cols = c(data)) 

ggplot(d, aes(ychange, lsr))+
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



# Species -----------------------------------------------------------------

# Using white shrimp, blue crap, brown shrimp

# 1. White Shrimp (Litopenaeus setiferus)
# Estuarine-dependent species
#
# Spawns offshore, juveniles migrate into estuaries
#
# Highly sensitive to salinity and temperature
#
# 🔶 2. Brown Shrimp (Farfantepenaeus aztecus) — Closely related competitor
# Similar habitat use to white shrimp
#
# Slightly different seasonal and salinity preferences
#
# Useful for testing niche differentiation or competitive overlap
#
# 🔷 3. Blue Crab (Callinectes sapidus) — Functionally different predator
# Also estuarine-dependent
#
# Omnivorous predator, unlike shrimp
#
# Occupies broader salinity/temperature range
#
# Useful for testing niche separation between trophic levels
#
# 🧠 This trio lets you ask:
#   "To what extent do co-occurring estuarine species with differing ecologies (two shrimp and one predator) share environmental niche space in Calcasieu estuary?"
#
# Or alternatively:
#
#   "Do closely related penaeid shrimp partition environmental niche space, and how does that compare to a more distantly related estuarine predator?"





#Filter out species

dfs_hvs <- df |>
  filter(species %in% c("white shrimp", "blue crab", "brown shrimp")) |>
  mutate(across(c(salinity, watertemp, turbidity, airtemp), scale)) |>
  select(species, salinity, watertemp, turbidity, airtemp, num) |>
  group_by(species) |>
  nest(weight = num, data = salinity:airtemp) |> 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, name = paste(species),
                                                     samples.per.point = 1000,
                                                     weight = weight$num,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         centroid = map(hv, \(hv) get_centroid(hv)),
         size = map_dbl(hv, \(hv) get_volume(hv)))

beepr::beep(3)

#saveRDS(dfs_hvs, 'data/first_hvs_species.rds')

dfs_hvs=readRDS('data/first_hvs_species.rds')

head(dfs_hvs)


plot(dfs_hvs$hv[[1]])
plot(dfs_hvs$hv[[2]])
plot(dfs_hvs$hv[[3]])


hvj = hypervolume_join(dfs_hvs$hv[[1]], dfs_hvs$hv[[2]], dfs_hvs$hv[[3]])
plot(hvj)


# comparison of across each species
df_y= tibble(y1 = unique(dfs_hvs$species),
             y2 = unique(dfs_hvs$species)) |> 
  expand(y1,y2)

# make all unique species comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df1 = dfs_hvs |> 
  select(y1 = species, hv1 = hv, hv1_size = size, cent1 = centroid)

df2 = dfs_hvs |> 
  select(y2 = species, hv2 = hv, hv2_size = size, cent2 = centroid)


# create data frame of all data and make species comparisons
df_ov = tibble(species = rep(unique(dfs_hvs$species),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(dfs_hvs$species))),
               y2 = rep(df_y$y2, times = length(unique(dfs_hvs$species)))) |> 
  inner_join(df1, by = c( 'y1')) |> 
  inner_join(df2, by = c( 'y2')) |> 
  mutate(
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
  select(y1, y2,lsr,
         dist_cent, jaccard, sorensen,
         uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2)

# save output
#write_csv(df_ov, 'data/species_centDist.csv')

head(df_ov)
# Species comparisons

df_ov=read_csv('data/species_centDist.csv')






dfs_hvs |>
  select(species, size) |>
  ggplot(aes(x = species, y = size, fill = species)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  labs(
    title = "Overall Niche Size of Estuarine Species",
    x = "Species",
    y = "Hypervolume Size"
  ) +
  theme_minimal(base_size = 14)+
  geom_text(aes(label = round(size, 3)), vjust = -0.3, size = 5)


dfs_hvs |>
  select(species, size) |>
  ggplot(aes(x = species, y = 1, size = size, fill = species)) +
  geom_point(shape = 21, color = "black") +
  geom_text(aes(label = round(size, 2)), vjust = -3, size = 9) +  # nicely formatted labels
  scale_size_continuous(range = c(10, 40)) +
  labs(
    title = "Niche Size of Estuarine Species",
    x = "Species",
    y = NULL
  ) +
  guides(size = "none", fill = "none") +  # remove legends
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )


df_ov_dedup <- df_ov |> 
  group_by(y1, y2) |> 
  summarise(across(where(is.numeric), mean), .groups = "drop")

ggplot(df_ov_dedup, aes(x = y1, y = y2, fill = sorensen)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(sorensen, 2)), size = 4) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1)) +
  labs(
    title = "Pairwise Niche Overlap (Sorensen Similarity)",
    x = "Species 1",
    y = "Species 2",
    fill = "Sorensen"
  ) +
  theme_classic(base_size = 14)

ggsave("sorensen_heatmap.png", plot = last_plot(), width = 10, height = 8, dpi = 300)



# Q 3 ---------------------------------------------------------------------



hv_species = dfs_hvs |> 
  select(species, hv)


# Keep only environmental data
env_scaled = df |>
  filter(species %in% c("white shrimp", "blue crab", "brown shrimp")) |>
  select(year, salinity, watertemp, turbidity, airtemp) |>
  drop_na() |>
  mutate(across(everything(), as.numeric)) |>
  mutate(across(salinity:airtemp, scale)) |> 
  group_by(year) |> 
  nest()



# Expand grid of species × years
suitability_df <- crossing(
  species = c("white shrimp", "blue crab", "brown shrimp"),
  year = unique(env_scaled$year)
)

# Join HVs and data
suitability_df <- suitability_df |>
  left_join(hv_species, by = "species") |>
  left_join(env_scaled, by = "year")

# Calculate % of suitable environments
suitability_df <- suitability_df |>
  mutate(prop_suitable = map2_dbl(hv, data, ~{
    points_inside <- hypervolume::hypervolume_inclusion_test(.x, .y)
    mean(points_inside)  # proportion of sites inside the niche
  }))

# Caused by warning in `hypervolume::hypervolume_inclusion_test()`:
#   ! Results may have a high error rate. Consider setting fast.or.accurate='accurate'.

#saveRDS(suitability_df, 'data/suitability_df.rds')

suitability_df = readRDS('data/suitability_df.rds')

ggplot(suitability_df, aes(x = year, y = prop_suitable, color = species)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "Year", y = "Proportion of Suitable Sites",
       title = "Habitat Suitability Over Time")


suitability_df2 <- suitability_df |>
  mutate(
    n_points = map_int(data, nrow),
    prop_suitable = ifelse(n_points < 20, NA, 
                           map2_dbl(hv, data, ~{
                             points_inside <- hypervolume::hypervolume_inclusion_test(.x, .y)
                             mean(points_inside)
                           }))
  )



species_counts <- df |>
  filter(species %in% c("white shrimp", "brown shrimp", "blue crab"),
         year %in% c("2015", "2016", "2017", "2018", "2019")) |>
  group_by(year, species) |>
  summarise(n_caught = sum(num), .groups = "drop") |>
  arrange(year, species)

print(species_counts)

target_years <- as.character(1986:2019)
target_species <- c("white shrimp", "brown shrimp", "blue crab")

# Fill in missing combinations with 0
species_counts <- df |>
  filter(species %in% target_species,
         year %in% target_years) |>
  group_by(year, species) |>
  summarise(n_caught = sum(num), .groups = "drop") |>
  complete(year = target_years, species = target_species, fill = list(n_caught = 0)) |>
  arrange(year, species)

ggplot(species_counts, aes(x = as.numeric(year), y = n_caught, color = species)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Abundance of Estuarine Species (1986–2019)",
    x = "Year",
    y = "Number Caught",
    color = "Species"
  ) +
  theme_classic(base_size = 14)


species_counts_zoom <- species_counts |>
  filter(as.numeric(year) >= 2017 & as.numeric(year) <= 2019)

ggplot(species_counts_zoom, aes(x = as.numeric(year), y = n_caught, color = species)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 2017:2019, labels = 2017:2019) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Abundance of Estuarine Species (2017–2019)",
    x = "Year",
    y = "Number Caught",
    color = "Species"
  ) +
  theme_classic(base_size = 14)




white_hv <- dfs_hvs |>
  filter(species == "white shrimp") |>
  pull(hv)

white_hv <- white_hv[[1]]

env_2019 <- env_scaled |>
  filter(year == "2019") |>
  pull(data)

env_2019 <- env_2019[[1]]

points_inside <- hypervolume::hypervolume_inclusion_test(white_hv, env_2019)
prop_suitable_2019 <- mean(points_inside)

points_inside <- hypervolume::hypervolume_inclusion_test(
  white_hv,
  env_2019,
  fast.or.accurate = "accurate"
)

prop_suitable_2019 <- mean(points_inside)


white_centroid <- dfs_hvs |>
  filter(species == "white shrimp") |>
  pull(centroid)

white_centroid <- white_centroid[[1]]

# Get means of scaled 2019 environmental variables
env_2019_mean <- colMeans(env_2019)

# Difference from centroid
env_2019_mean - white_centroid


env_2019 <- env_scaled |>
  filter(year == "2019") |>
  pull(data)

env_2019_data <- env_2019[[1]]

suitability_2019_summary <- dfs_hvs |>
  select(species, hv, centroid) |>
  mutate(
    env_2019 = list(env_2019_data),
    
    prop_suitable = map2_dbl(hv, env_2019, ~ {
      if (nrow(.y) < 5) return(NA_real_)  # safeguard
      points_inside <- hypervolume::hypervolume_inclusion_test(.x, .y, fast.or.accurate = "accurate")
      mean(points_inside)
    }),
    
    env_2019_mean = map(env_2019, colMeans),
    
    diff_from_centroid = map2(env_2019_mean, centroid, ~ .x - .y)
  ) |>
  select(species, prop_suitable, diff_from_centroid) |>
  unnest_wider(diff_from_centroid)


write_csv(suitability_2019_summary, "data/suitability_2019_summary.csv")

