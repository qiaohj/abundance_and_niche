setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
library(dplyr)
library("mixtools")
df<-readRDS("../Tables/box.rda")
df$mh_dist<-NA
df$qchisq<-NA
df<-as.data.frame(df)
i=925
for (i in c(1:nrow(df))){
  print(paste(i, nrow(df)))
  item<-df[i,]
  mve<-readRDS(sprintf("../Models/points/%d_%d/mve.rda", item$semi_ra, item$index))
  df[i, ]$mh_dist <- stats::mahalanobis(item[, c("PC_1", "PC_2")], center = mve$centroid, 
                                cov = mve$covariance)
  df[i, ]$qchisq<-stats::qchisq(0.95, length(mve$centroid))
}
df$relavent_mh_dist<-df$mh_dist/stats::qchisq(0.95, length(mve$centroid))
df[132,]
mu <- c(1, 3)
sigma <- matrix(c(1, .3, .3, 1.5), 2, 2)


ellipse(mve$centroid, mve$covariance, npoints = 10000, newplot = TRUE, pch=".")
points(x=mve$centroid[1], y=mve$centroid[2], col="red")
points(x=item[, "PC_1"], y=item[, "PC_2"], col="blue")

stats::mahalanobis(item[, c("PC_1", "PC_2")], center = mve$centroid, 
                              cov = mve$covariance)
stats::mahalanobis(item[, c("PC_2", "PC_1")], center = mve$centroid, 
                   cov = mve$covariance)
lines(x=mve$axis_coordinates[[1]][,1], y=mve$axis_coordinates[[1]][,2])
lines(x=mve$axis_coordinates[[2]][,1], y=mve$axis_coordinates[[2]][,2])

df<-as_tibble(df)
saveRDS(df, "../Tables/box_with_mh_dist.rda")
hist(df$mh_dist)

item<-df[which(df$mh_dist==max(df$mh_dist, na.rm = T)),]
item$mh_dist
item$land_percentile
points<-readRDS(sprintf("../Models/points/%d_%d/points.rda", item$semi_ra, item$index))
                
points_no_NA<-points[complete.cases(points),]
fit <- cov.rob(points_no_NA[, c("V_1", "V_2")], quantile.used=NDquntil(nrow(points_no_NA), 0.95),  method = "mve")
best_ellipse <- ellipsoidhull( as.matrix(points_no_NA[fit$best, c("V_1", "V_2")] ))
plot(points_no_NA$V_1, points_no_NA$V_2, pch=".")
lines(predict(best_ellipse), col="blue")
points(df[1, "PC_1"], df[1, "PC_2"], col="red")
points(fit$center[1], fit$center[2], col="blue")


df_sub<-df%>%filter(land_percentile>0.9)
hist(df_sub$mh_dist)
source("function.r")
df_se<-summarySE(df_sub, "mh_dist", "semi_ra")

ggplot(df_sub) + geom_point(aes(x=mh_dist, y=alt_sd, color=factor(semi_ra))) + xlim(c(0, 6))
ggplot(df_se) + geom_point(aes(x=factor(semi_ra), y=mean))

df_sub %>%
  group_by(semi_ra) %>%
  group_map(~ cor(.x$mh_dist, .x$alt_sd))
ggplot(df_sub) + geom_point(aes(x=mh_dist, y=relavent_dist, color=factor(semi_ra))) 

