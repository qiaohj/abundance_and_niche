library(ggplot2)
library(dplyr)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
df<-readRDS("../Tables/box.rda")
ggplot(df) + geom_point(aes(x=width, y=relavent_dist, color=factor(index))) + 
  geom_line(aes(x=width, y=relavent_dist, color=factor(index)))+
  geom_smooth(aes(x=width, y=relavent_dist), method = lm)

ggplot(df) + geom_point(aes(x=width, y=V1_sd, color=factor(index))) + 
  geom_line(aes(x=width, y=V1_sd, color=factor(index)))

ggplot(df) + geom_point(aes(x=V1_sd, y=relavent_dist, color=factor(index))) + 
  geom_line(aes(x=V1_sd, y=relavent_dist, color=factor(index))) + 
  geom_smooth(aes(x=V1_sd, y=relavent_dist), method = lm)

ggplot(df) + geom_point(aes(x=M_V1_ovserved, y=relavent_dist, color=factor(index))) + 
  geom_line(aes(x=M_V1_ovserved, y=relavent_dist, color=factor(index))) + 
  geom_smooth(aes(x=M_V1_ovserved, y=relavent_dist), method = lm)


ggplot(df) + geom_point(aes(x=size_continent, y=relavent_dist, color=factor(index))) + 
  geom_line(aes(x=size_continent, y=relavent_dist, color=factor(index)))
ggplot(df) + geom_point(aes(x=width, y=M_V1_ovserved))
ggplot(df) + geom_point(aes(x=width, y=V1_sd))
ggplot(df) + geom_point(aes(x=width, y=land_percentile))
ggplot(df) + geom_point(aes(x=width, y=M_V2_ovserved))

ggplot(df[which(df$land_percentile>0.7),]) +geom_density(aes(relavent_dist, color=factor(width)))

points<-readRDS("C:/Users/Huijie Qiao/Downloads/Tables/random_points_100.rda")
dim(points)

100, 1000, by=100


plot(points$x, points$y)
