

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

head(df)



df_before = df |> 
  group_by(year,month,day) |> 
  summarise(salinity = mean(salinity, na.rm = T),
            watertemp = mean(watertemp, na.rm = T),
            turbidity = mean(turbidity, na.rm = T),
            airtemp = mean(airtemp, na.rm = T)) |> 
  ungroup() |>
  filter(year == "1986") |> 
  mutate(across(salinity:airtemp, scale)) |> 
  select(salinity, watertemp, turbidity,airtemp) 


hv_before = hypervolume_gaussian(df_before, name = 'Before',
                                 samples.per.point = 1000,
                                 kde.bandwidth = estimate_bandwidth(df_before), 
                                 sd.count = 3, 
                                 quantile.requested = 0.95, 
                                 quantile.requested.type = "probability", 
                                 chunk.size = 1000, 
                                 verbose = F) 




hv_before
plot(hv_before)
plot(hv_before, show.3d = TRUE)


df2 = df |> 
  group_by(year,month,day) |> 
  summarise(salinity = mean(salinity, na.rm = T),
            watertemp = mean(watertemp, na.rm = T),
            turbidity = mean(turbidity, na.rm = T),
            airtemp = mean(airtemp, na.rm = T)) |> 
  ungroup() |> 
  select(year,salinity, watertemp, turbidity,airtemp) |> 
  mutate(across(salinity:airtemp, scale))  |> 
  group_by(year) |> 
  nest()

df2

df2 = df2 |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = year,
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)))

head(df2)

plot(df2$hv[[5]], show.3d = F, 
     xlim = c(-3,3), ylim = c(-3,3), zlim = c(-3,3),
     xlab = "Salinity", ylab = "Water Temp", zlab = "Turbidity",
     main = "1986" )



# Lets do this ------------------------------------------------------------
set.seed(1)

df2 <- df |> 
  group_by(lat, lon, year, month, day) |> 
  slice(1) |> 
  ungroup() |>
  mutate(across(salinity:airtemp, scale),
         across(c(salinity:airtemp), 
                ~map_dbl(., ~. + rnorm(1, mean = 0, sd = 0.0001)))) |> 
  select(year,salinity, watertemp, turbidity,airtemp) |> 
  group_by(year) |> 
  nest()

head(df2)

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

saveRDS(df2, 'data/first_hvs.rds') 

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
write_csv(df_ov, 'data/first_centDist.csv')

df_ov = read_csv('data/first_centDist.csv')

df_ov
