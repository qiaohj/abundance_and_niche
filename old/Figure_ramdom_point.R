library(raster)
library(ggplot2)


library(dplyr)
library(raster)
library(Hmisc)
library(stringr)
library(ntbox)
library(hypervolume) 
library(fields)
library(ggfortify)
library(GWmodel)
library(ape)
library(MASS)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
vars<-c(1, 2)
i=1
j=1

k=1
template<-"../Raster/PCs/pc%d.tif"
r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, i))
  r_env[[as.character(i)]]<-r
}

#Filted something
first_round<-readRDS("../Tables/box_with_mh_dist.rda")
my_thd=0.90#  
first_round <- subset(first_round,width<=1001)   #
NUM_width <- length(first_round$width[first_round$index==1])
#There are 29 different width
freq1 <- data.frame(table(first_round$index[first_round$land_percentile>=my_thd]))
good_case <- as.character(freq1$Var1[freq1$Freq==NUM_width])
print(length(good_case))
#first_round$index
dd <- subset(first_round, index %in% good_case)
labels<-paste(dd$x, dd$y)

centers<-readRDS("../Tables/random_points_100.rda")
centers$index<-c(1:nrow(centers))
#centers<-subset(centers, index %in% good_case)     

#plot(centers$x, centers$y)
ratios<-seq(100, 500, by=100)
ratios<-c(seq(10, 90, by=10), ratios)
ratio = ratios[11]
i = 1
res = 1000
#pca<-readRDS("../Models/pca_model.rda")
#x, y, Niche_PC1, Niche_PC2, x_PC1, y_PC2, Niche_Sp_distance (ABS), r_a,r_b,  
#Niche_Sp_distance *2/(r_a+r_b), land size %, width of square, env_heterogenity or SD,
base<-"../Models/points/%d_%d"
r_r<-raster("../Raster/alt_eck4_1km_new.tif")
#result<-readRDS("../Tables/box.rda")
i=1

df_box<-readRDS("../Tables/box.rda")
filter_df_box<-df_box%>%dplyr::filter((semi_ra==ratio)&((alt_sd>=1812)|(alt_sd<=14.3)))
filter_df_box<-df_box%>%dplyr::filter((semi_ra==ratio)&((relavent_dist>=2.545176976)|(relavent_dist<=0.006527766)))

range(df_box[which(df_box$semi_ra==ratio), ]$relavent_dist, na.rm = T)

png(filename="../Figures/overview_relavent_dist.png")
plot(r_r, xlim=c(-11e6,-5e6), ylim=c(2e6, 6e6))
#plot(r_r)
center1<-filter_df_box[1,]
lines(x=c(center1$x-ratio*res, center1$x+ratio*res, center1$x+ratio*res, center1$x-ratio*res, center1$x-ratio*res),
      y=c(center1$y-ratio*res, center1$y-ratio*res, center1$y+ratio*res, center1$y+ratio*res, center1$y-ratio*res))
center2<-filter_df_box[2,]
lines(x=c(center2$x-ratio*res, center2$x+ratio*res, center2$x+ratio*res, center2$x-ratio*res, center2$x-ratio*res),
      y=c(center2$y-ratio*res, center2$y-ratio*res, center2$y+ratio*res, center2$y+ratio*res, center2$y-ratio*res))

dev.off()

center<-filter_df_box[1,]
points<-expand.grid(x=seq(from=center$x - ratio*res, to=center$x + ratio * res, by=res),
                    y=seq(from=center$y - ratio*res, to=center$y + ratio * res, by=res))

for (j in vars){
  center[, sprintf("PC_%d", j)]<-extract(r_env[[as.character(j)]], center[, c("x", "y")])
  print(paste("ratio:", ratio, i, nrow(centers), "Var:", j, sep="/"))
  points[, sprintf("V_%d", j)]<-extract(r_env[[as.character(j)]], points[, c("x", "y")])
}

points$alt<-extract(r_r, points[, c("x", "y")])

points_no_NA<-points[complete.cases(points),]
size_continent<-nrow(points_no_NA)

mve<-cov.mve(points_no_NA[, c("V_1", "V_2")])

center$index<-i
center$semi_ra<-ratio
center$Niche_PC1<-mve$center[1]
center$Niche_PC2<-mve$center[2]
points_no_NA$dist<-sqrt((points_no_NA$V_1-center$Niche_PC1)^2+(points_no_NA$V_2-center$Niche_PC2)^2)
points_no_NA$mhdist<-mahalanobis(points_no_NA[, c("V_1", "V_2")], mve$center, mve$cov)
niche_centers<-points[which(points_no_NA$mhdist<quantile(points_no_NA$mhdist, 0.01)),]


png(filename="../Figures/diagram_min_relavent_dist.png", width=500, height=500)
plot(r_r, xlim=c(center$x-ratio*res,center$x+ratio*res), ylim=c(center$y-ratio*res, center$y+ratio*res))
lines(x=c(center$x-ratio*res, center$x+ratio*res, center$x+ratio*res, center$x-ratio*res, center$x-ratio*res),
      y=c(center$y-ratio*res, center$y-ratio*res, center$y+ratio*res, center$y+ratio*res, center$y-ratio*res))

points(niche_centers$x, niche_centers$y, col="red", pch=".")
points(center$x, center$y)
dev.off()

source("addEllipse.R")
png(filename="../Figures/diagram_min_relavent_dist_E.png", width=500, height=500)
plot(points_no_NA$V_1, points_no_NA$V_2, pch=".", col=alpha("grey", 0.4), 
     xlim=c(-0.3, 1.5), ylim=c(-1, 1.5))
addEllipse(mve$center, mve$cov, col="red", p.interval=0.95)
points(niche_centers$V_1, niche_centers$V_2, col=alpha("pink", 0.2), pch=".")
points(center$PC_1, center$PC_2, pch=2)
points(center$Niche_PC1, center$Niche_PC2, pch=1, col="red")
dev.off()
###End of Figure 1
###Figure 2

center<-filter_df_box[2,]
points<-expand.grid(x=seq(from=center$x - ratio*res, to=center$x + ratio * res, by=res),
                    y=seq(from=center$y - ratio*res, to=center$y + ratio * res, by=res))

for (j in vars){
  center[, sprintf("PC_%d", j)]<-extract(r_env[[as.character(j)]], center[, c("x", "y")])
  print(paste("ratio:", ratio, i, nrow(centers), "Var:", j, sep="/"))
  points[, sprintf("V_%d", j)]<-extract(r_env[[as.character(j)]], points[, c("x", "y")])
}

points$alt<-extract(r_r, points[, c("x", "y")])

points_no_NA<-points[complete.cases(points),]
size_continent<-nrow(points_no_NA)

mve<-cov.mve(points_no_NA[, c("V_1", "V_2")])

center$index<-i
center$semi_ra<-ratio
center$Niche_PC1<-mve$center[1]
center$Niche_PC2<-mve$center[2]
points_no_NA$dist<-sqrt((points_no_NA$V_1-center$Niche_PC1)^2+(points_no_NA$V_2-center$Niche_PC2)^2)
points_no_NA$mhdist<-mahalanobis(points_no_NA[, c("V_1", "V_2")], mve$center, mve$cov)

niche_centers<-points[which(points_no_NA$mhdist<quantile(points_no_NA$mhdist, 0.01)),]

png(filename="../Figures/diagram_max_relavent_dist.png", width=500, height=500)
plot(r_r, xlim=c(center$x-ratio*res,center$x+ratio*res), ylim=c(center$y-ratio*res, center$y+ratio*res))
lines(x=c(center$x-ratio*res, center$x+ratio*res, center$x+ratio*res, center$x-ratio*res, center$x-ratio*res),
      y=c(center$y-ratio*res, center$y-ratio*res, center$y+ratio*res, center$y+ratio*res, center$y-ratio*res))

points(niche_centers$x, niche_centers$y, col="red", pch=".")
points(center$x, center$y)
dev.off()

source("addEllipse.R")
png(filename="../Figures/diagram_max_relavent_dist_E.png", width=500, height=500)
plot(points_no_NA$V_1, points_no_NA$V_2, pch=".", col=alpha("grey", 0.4), 
     xlim=c(-1, 0.2), ylim=c(-1, 3))
addEllipse(mve$center, mve$cov, col="red", p.interval=0.95)
points(niche_centers$V_1, niche_centers$V_2, col=alpha("pink", 0.2), pch=".")
points(center$PC_1, center$PC_2, pch=2)
points(center$Niche_PC1, center$Niche_PC2, pch=1, col="red")

dev.off()
library(ggplot2)
library(dplyr)
df_box_se<-df_box%>%dplyr::group_by(semi_ra)%>%
  dplyr::summarise()
p<-ggplot(df_box, aes(x=semi_ra, y=mh_dist))+geom_point()+geom_smooth()+ylim(0, 10)+theme_bw()
ggsave(p, file="../Figures/radius_mh_dist.png")
p<-ggplot(df_box, aes(x=alt_sd, y=mh_dist))+geom_point()+geom_smooth()+ylim(0, 10)+theme_bw()
ggsave(p, file="../Figures/alt_sd_mh_dist.png")



