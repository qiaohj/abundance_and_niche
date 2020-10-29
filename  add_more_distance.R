setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
library(dplyr)
library("mixtools")
df<-readRDS("../Tables/box.rda")

df<-as.data.frame(df)
i=925
df$g_dist_g_center_mean<-NA
df$g_dist_g_center_sd<-NA
df$mh_dist_g_center_mean<-NA
df$mh_dist_g_center_sd<-NA

for (i in c(1:nrow(df))){
  print(paste(i, nrow(df)))
  item<-df[i,]
  mve<-readRDS(sprintf("../Models/points/%d_%d/mve.rda", item$semi_ra, item$index))
  points<-readRDS(sprintf("../Models/points/%d_%d/points.rda", item$semi_ra, item$index))
  points$mh_dist<-stats::mahalanobis(points[, c("V_1", "V_2")], center = mve$centroid, 
                                     cov = mve$covariance)
  threshold<-quantile(points$mh_dist, 0.01, na.rm=T)
  points<-points[which(points$mh_dist<=threshold),]
  points$g_dist_g_center<-sqrt((points$x-item$x)^2+(points$y-item$y)^2)
  g_center<-c(item$PC_1, item$PC_2)
  points$mh_dist_g_center <- stats::mahalanobis(points[, c("V_1", "V_2")], 
                                               center = g_center, 
                                               cov = mve$covariance)
  
  df[i, ]$g_dist_g_center_mean <- mean(points$g_dist_g_center)
  df[i, ]$g_dist_g_center_sd <- sd(points$g_dist_g_center)
  df[i, ]$mh_dist_g_center_mean <- mean(points$mh_dist_g_center)
  df[i, ]$mh_dist_g_center_sd <- sd(points$mh_dist_g_center)
}

saveRDS(df, "../Tables/box_with_more_dist.rda")
