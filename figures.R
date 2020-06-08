library(raster)
library(ggplot2)
library(ntbox)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
df<-readRDS("../Tables/box.rda")
df$relavent_dist_x<-abs(df$Niche_PC1-df$PC_1)/df$r_a
df$relavent_dist_y<-abs(df$Niche_PC2-df$PC_2)/df$r_b
df$relavent_dist_y<-abs(df$Niche_PC2-df$PC_2)/df$r_b


hist(df$relavent_dist_x)
hist(df$relavent_dist_y)
plot(df$relavent_dist_x, df$relavent_dist_y)

points<-readRDS("../Models/points/100_1/points.rda")
points_no_NA<-points[complete.cases(points),]
mve<-cov_center(data = points_no_NA, mve = T, 
                level = 0.95, vars = 3:4)
df[1,]
vars<-c(1, 2)
template<-"../Raster/PCs/pc%d.tif"
r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, i))
  r_env[[as.character(i)]]<-r
}
extract(r_env[[1]], df[1, c("x", "y")])
extract(r_env[[2]], df[1, c("x", "y")])

library(car)
library(MASS)
library(cluster)

NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}


fit <- cov.rob(points_no_NA[, c("pc1", "pc2")], quantile.used=NDquntil(nrow(points_no_NA), 0.95),  method = "mve")
best_ellipse <- ellipsoidhull( as.matrix(points_no_NA[fit$best, c("pc1", "pc2")] ))
plot(points$pc1, points$pc2, pch=".")
lines(predict(best_ellipse), col="blue")
points(df[1, "PC_1"], df[1, "PC_2"], col="red")
points(fit$center[1], fit$center[2], col="blue")

mh_dist <- stats::mahalanobis(df[1, c("PC_1", "PC_2")], center = mve$centroid, 
                              cov = mve$covariance)
stats::mahalanobis(df[1, c("PC_2", "PC_1")], center = mve$centroid, 
                   cov = mve$covariance)
stats::mahalanobis(df[1, c("Niche_PC1", "Niche_PC2")], center = mve$centroid, 
                   cov = mve$covariance)

in_Ellipsoid <- mh_dist <= stats::qchisq(0.95, length(mve$centroid))
