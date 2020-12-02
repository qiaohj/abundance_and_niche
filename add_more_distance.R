setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
library(dplyr)
library(Rmisc)
library(mixtools)
df<-readRDS("../Tables/box.rda")

df<-as.data.frame(df)
i=925
df$centers_g_2_centers_n_mean<-NA
df$centers_g_mean<-NA
df$centers_n_2_centers_g_mean<-NA
df$centers_n_mean<-NA

df$centers_g_2_centers_n_sd<-NA
df$centers_g_sd<-NA
df$centers_n_2_centers_g_sd<-NA
df$centers_n_sd<-NA

df$centers_g_2_centers_n_ci<-NA
df$centers_g_ci<-NA
df$centers_n_2_centers_g_ci<-NA
df$centers_n_ci<-NA

df$centers_n_in_centers_g<-NA
df$centers_n_out_centers_g<-NA

df$centers_g_in_centers_n<-NA
df$centers_g_out_centers_n<-NA

for (i in c(1:nrow(df))){
  print(paste(i, nrow(df)))
  item<-df[i,]
  mve<-readRDS(sprintf("../Models/points/%d_%d/mve.rda", item$semi_ra, item$index))
  points<-readRDS(sprintf("../Models/points/%d_%d/points.rda", item$semi_ra, item$index))
  #get 1% points around the center
  points$g_dist_g_center<-sqrt((points$x-item$x)^2+(points$y-item$y)^2)
  points$mh_dist<-stats::mahalanobis(points[, c("V_1", "V_2")], center = mve$centroid, 
                                     cov = mve$covariance)
  
  threshold_g<-quantile(points$g_dist_g_center, 0.01, na.rm=T)
  centers_g<-points%>%dplyr::filter(g_dist_g_center<=threshold_g)
  
  
  threshold_n<-quantile(points$mh_dist, 0.01, na.rm=T)
  centers_n<-points%>%dplyr::filter(mh_dist<=threshold_n)

  df[i, ]$centers_g_2_centers_n_mean <- mean(centers_g$mh_dist)
  df[i, ]$centers_g_mean <- mean(centers_g$g_dist_g_center)
  df[i, ]$centers_n_2_centers_g_mean <- mean(centers_n$g_dist_g_center)
  df[i, ]$centers_n_mean <- mean(centers_n$mh_dist)
  df[i, ]$centers_g_2_centers_n_sd <- sd(centers_g$mh_dist)
  df[i, ]$centers_g_sd <- sd(centers_g$g_dist_g_center)
  df[i, ]$centers_n_2_centers_g_sd <- sd(centers_n$g_dist_g_center)
  df[i, ]$centers_n_sd <- sd(centers_n$mh_dist)
  
  df[i, ]$centers_g_2_centers_n_ci <- CI(centers_g$mh_dist)[1] - CI(centers_g$mh_dist)[2]
  df[i, ]$centers_g_ci <- CI(centers_g$g_dist_g_center)[1] - CI(centers_g$g_dist_g_center)[2]
  df[i, ]$centers_n_2_centers_g_ci <- CI(centers_n$g_dist_g_center)[1] - CI(centers_n$g_dist_g_center)[2]
  df[i, ]$centers_n_ci <- CI(centers_n$mh_dist)[1] - CI(centers_n$mh_dist)[2]
  
  df[i, ]$centers_n_in_centers_g<-nrow(centers_n%>%dplyr::filter(centers_n$g_dist_g_center<=threshold_g))
  df[i, ]$centers_n_out_centers_g<-nrow(centers_n%>%dplyr::filter(centers_n$g_dist_g_center>threshold_g))
  
  df[i, ]$centers_g_in_centers_n<-nrow(centers_g%>%dplyr::filter(centers_g$mh_dist<=threshold_n))
  df[i, ]$centers_g_out_centers_n<-nrow(centers_g%>%dplyr::filter(centers_g$mh_dist>threshold_n))
  
  if (F){
    plot(points$x, points$y)
    points(item$x, item$y, col="red")
  }
}

saveRDS(df, "../Tables/box_with_more_dist.rda")
