# title: "Workshop 1: Introduction to hypervolumes"
# author: "BSC6926-B53B"
# date: "3/14/25"


## Getting to know the basics 

## Hypervolumes 
# Hypervolumes are a multidimensional tool that is based on Hutchinson's *n*-dimensional niche concept and we can build them with the `hypervolume` package. 
# 
# ### Preparing the data
# Typically we have a dataset that has the abundance of the community and another that has trait data. Therefore we need to combine them. Also, we can't have `NAs`, so we have to filter any missing data one we combine the datasets. 

library(tidyverse)
# abundance data
ab = read_csv('data/abundHermine.csv')

# trait data
tr = read_csv('data/fishTraits.csv')

# combine
df = left_join(ab, tr, by = 'species') |> 
  drop_na()

# Now we can make the hypervolume by selecting the data to be included and z-scoring

df_before = df |> 
  filter(period == 'Before') |> 
  mutate(across(trophic_level:generation_time, scale)) |> 
  select(trophic_level, temp_preference, generation_time)

df_before

### Building hypervolumes
# With a nested dataset of our columns that we want to build hypervolumes for we can use `mutate()` and `map()` to generate the hypervolume. 

library(hypervolume)
#library(rgl) #for plotting in 3D

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

# This can also be done to multiple at the same time with nest().

df = df |> 
  select(period, trophic_level, temp_preference, generation_time) |> 
  mutate(across(trophic_level:generation_time, scale)) |> 
  group_by(period) |> 
  nest()

df

df = df |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = period,
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)))

head(df)
df$hv

### plotting hypervolumes 
# We can plot multiple hypervolumes by joining them together 

hvj = hypervolume_join(df$hv[[1]], df$hv[[2]])

#plot 3d
plot(hvj, pairplot = T, colors=c('red','black'),
     show.3d=T,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=T,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

#plot biplot
plot(hvj, pairplot = T, colors=c('dodgerblue3','cyan2'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
### hypervolume metrics
# The geometry of hypervolumes are useful when characterizing and comparing hypervolumes. Hypervolume size represents the variation of the data, centroid distance compares the euclidian distance between two hypervolume centroids (mean conditions), and overlap measures the simularity of hypervolumes. 

# size 
df = df |> 
  mutate(hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         hv_cent = map(hv, \(hv) get_centroid(hv)))

head(df)

# centroid distance 
hypervolume_distance(df$hv[[1]], df$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset = hypervolume_set(df$hv[[1]], df$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset)


### Weight hypervolume input
# The above hypervolume is just based on the traits using presence of species, but we can weight the points to shape the hypervolume based on abundance 

#prep data
df_w = left_join(ab, tr, by = 'species') |> 
  drop_na() |> 
  select(period, abund, trophic_level, temp_preference, generation_time) |> 
  mutate(across(trophic_level:generation_time, scale)) |> 
  group_by(period) |> 
  nest(weight = abund, data = trophic_level:generation_time) 
df_w

# make hypervolumes
df_w = df_w |> 
  mutate(hv = map2(data,weight, \(data,weight) 
                   hypervolume_gaussian(data, 
                                        name = paste(period,'weighted',sep = '_'),
                                                                    weight = weight$abund,
                                                                    samples.per.point = 1000,
                                                                    kde.bandwidth = estimate_bandwidth(data), 
                                                                    sd.count = 3, 
                                                                    quantile.requested = 0.95, 
                                                                    quantile.requested.type = "probability", 
                                                                    chunk.size = 1000, 
                                                                    verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)))

df_w


hvj_w = hypervolume_join(df_w$hv[[1]], df_w$hv[[2]])

plot(hvj_w, pairplot = T, colors=c('dodgerblue3','cyan2'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=T,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

# centroid distance 
hypervolume_distance(df_w$hv[[1]], df_w$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset_w = hypervolume_set(df_w$hv[[1]], df_w$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset_w)

### From mean and sd
# Sometimes, we do not have enough points to meet assumptions to make a hypervolume. Therefore, we can simulate random points based on mean and sd of our axes. We can then simulate the information needed and make our hypervolumes.

# mean and sd
df_m = left_join(ab, tr, by = 'species') |> 
  drop_na() |> 
  pivot_longer(trophic_level:generation_time, names_to = 'trait', values_to = 'value') |> 
  group_by(period,trait) |> 
  summarize(mean = mean(value),
            sd = sd(value))

df_m

#generate points 
# number of points 
n = 50 

df_tot = df_m |> slice(rep(1:n(), each=n))|>
  mutate(point = map2_dbl(mean,sd, \(mean,sd) rnorm(1,mean =mean,sd =sd))) |> 
  group_by(period, trait) |> 
  mutate(num = row_number()) |>
  ungroup() |> 
  select(-mean, -sd)|>
  pivot_wider(names_from = trait, values_from = point)|> 
  select(-num) |> 
  mutate(across(generation_time:trophic_level,scale)) |> 
  group_by(period) |> 
  nest()


# generate hypervolumes
df_tot = df_tot |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = period,
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)))


#plot
hvj_tot = hypervolume_join(df_tot$hv[[1]], df_tot$hv[[2]])

plot(hvj_tot, pairplot = T, colors=c('dodgerblue3','cyan2'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=T,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

# centroid distance 
hypervolume_distance(df_tot$hv[[1]], df_tot$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset_tot = hypervolume_set(df_tot$hv[[1]], df_tot$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset_tot)

