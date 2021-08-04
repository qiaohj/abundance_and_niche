library(raster)
library(ggplot2)
library(dplyr)
library(MASS)
library(ggpubr)
library(ggpmisc)
library(Rmisc)
library(vegan)
library(ggpmisc)
library(viridis)

setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
#center points
p_sub<-readRDS("../Tables/random_points_1000_0.8.rda")
if (F){
  mask<-raster("../Raster/alt_eck4_1km_new.tif")
  mask_low<-mask
  res(mask_low)<-c(50000, 50000)
  mask_low<-resample(mask, mask_low)
  plot(mask_low)
  writeRaster(mask_low, "../Raster/alt_eck4_50km.tif", overwrite=T)
}
mask<-raster("../Raster/alt_eck4_50km.tif")
p_mask<-data.frame(rasterToPoints(mask))
colnames(p_mask)[3]<-"v"
map_background<-"white"
map_theme<-theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_rect(fill = map_background, color = NA), 
  panel.background = element_blank(), 
  legend.background = element_rect(fill = map_background, color = NA),
  panel.border = element_blank(),
  legend.position="none"
)



ratio = 400
res = 1000
df_box<-readRDS("../Tables/box_with_more_dist.rda")
df_box<-df_box%>%dplyr::filter(land_percentile>=0.8)
df_box<-df_box%>%dplyr::filter(df_box$semi_ra>10)
df_box$xxx<-df_box$centers_n_in_centers_g/
  (df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)

filter_df_box_1<-df_box%>%dplyr::filter((semi_ra==ratio)&(xxx>=0.5)&(land_percentile==1))
df_box$dist_ttt<-sqrt((df_box$x-filter_df_box_1[3, "x"])^2+
                        (df_box$y-filter_df_box_1[3, "y"])^2)

filter_df_box_2<-df_box%>%dplyr::filter((semi_ra==ratio)&((xxx<=0.05)&(dist_ttt<=2e6)))
filter_df_box_2<-filter_df_box_2%>%dplyr::select(-dist_ttt)
filter_df_box<-rbind(filter_df_box_1[3,], filter_df_box_2[7,])
#hist(filter_df_box_2$dist_ttt)
hist(df_box$xxx)

#View(df_box[which(df_box$index==906),])

range(df_box[which(df_box$semi_ra==ratio), ]$relavent_dist, na.rm = T)

center1<-filter_df_box[1,]
lines1<-data.frame(x=c(center1$x-ratio*res, center1$x+ratio*res, center1$x+ratio*res, center1$x-ratio*res, center1$x-ratio*res),
                   y=c(center1$y-ratio*res, center1$y-ratio*res, center1$y+ratio*res, center1$y+ratio*res, center1$y-ratio*res))
center2<-filter_df_box[2,]
lines2<-data.frame(x=c(center2$x-ratio*res, center2$x+ratio*res, center2$x+ratio*res, center2$x-ratio*res, center2$x-ratio*res),
                   y=c(center2$y-ratio*res, center2$y-ratio*res, center2$y+ratio*res, center2$y+ratio*res, center2$y-ratio*res))
terrain.colors<-c("#ccece6", "#99d8c9", "#66c2a4", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f")
p<-ggplot()+geom_tile(data=p_mask, aes(x=x, y=y, fill=v), alpha=0.8)+
  geom_point(data=p_sub, aes(x=x, y=y), color="#000000", size=0.3, shape=4)+
  scale_fill_gradientn(colours = terrain.colors)+
  #scale_fill_viridis()+
  map_theme
p


p<-p+geom_path(data=lines1, aes(x=x, y=y), color="#000000")+
  geom_point(data=center1, aes(x=x, y=y), size=0.5)+
  geom_text(data=center1, aes(x=x, y=y, label="(1)"), color="#000000", vjust=2.5)


p<-p+geom_path(data=lines2, aes(x=x, y=y), color="#000000")+
  geom_point(data=center2, aes(x=x, y=y), size=0.5)+
  geom_text(data=center2, aes(x=x, y=y, label="(2)"), color="#000000", vjust=2.5)

p
if (F){
  ggsave(p, file="../Figures/center_map_example/center_map_example_top.pdf", width=12, height=6)
  ggsave(p, file="../Figures/center_map_example/center_map_example_top.png", width=12, height=6)
}

vars<-c(1, 2)
template<-"../Raster/PCs/pc%d.tif"
centers<-readRDS("../Tables/random_points_100.rda")
centers$index<-c(1:nrow(centers))

r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, i))
  r_env[[as.character(i)]]<-r
}

p_g_item<-list()
p_n_item<-list()
index=1
source("commonFuns/addEllipse.R")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#800026", "#CC79A7")

for (index in c(1:2)){
  center<-filter_df_box[index,]
  points<-expand.grid(x=seq(from=center$x - ratio*res, to=center$x + ratio * res, by=res),
                      y=seq(from=center$y - ratio*res, to=center$y + ratio * res, by=res))
  
  
  for (j in vars){
    center[, sprintf("PC_%d", j)]<-raster::extract(r_env[[as.character(j)]], center[, c("x", "y")])
    print(paste("ratio:", ratio, i, nrow(centers), "Var:", j, sep="/"))
    points[, sprintf("V_%d", j)]<-raster::extract(r_env[[as.character(j)]], points[, c("x", "y")])
  }
  
  points_no_NA<-points[complete.cases(points),]
  size_continent<-nrow(points_no_NA)
  
  mve<-cov.mve(points_no_NA[, c("V_1", "V_2")])
  
  points_no_NA$g_dist_g_center<-sqrt((points_no_NA$x-center$x)^2+(points_no_NA$y-center$y)^2)
  points_no_NA$mh_dist<-stats::mahalanobis(points_no_NA[, c("V_1", "V_2")], center = mve$center, 
                                           cov = mve$cov)
  
  threshold_g<-quantile(points_no_NA$g_dist_g_center, 0.01, na.rm=T)
  centers_g<-points_no_NA%>%dplyr::filter(g_dist_g_center<=threshold_g)
  
  
  threshold_n<-quantile(points_no_NA$mh_dist, 0.01, na.rm=T)
  centers_n<-points_no_NA%>%dplyr::filter(mh_dist<=threshold_n)
  
  r_r<-raster("../Raster/alt_eck4_1km_new.tif")
  points_no_NA$alt<-raster::extract(r_r, points_no_NA[, c("x", "y")])
  
  
  
  center$index<-index
  center$semi_ra<-ratio
  center$Niche_PC1<-mve$center[1]
  center$Niche_PC2<-mve$center[2]
  lines<-data.frame(x=c(center$x-ratio*res, center$x+ratio*res, center$x+ratio*res, center$x-ratio*res, center$x-ratio*res),
                    y=c(center$y-ratio*res, center$y-ratio*res, center$y+ratio*res, center$y+ratio*res, center$y-ratio*res))
  
  p1<-ggplot()+geom_tile(data=points_no_NA, aes(x=x, y=y, fill=alt), alpha=0.8)+
    scale_fill_gradientn(colours = terrain.colors)+
    geom_path(data=lines, aes(x=x, y=y), color="black")+
    geom_tile(data=centers_n, aes(x=x, y=y), fill="#800026",size=0.2)+
    geom_tile(data=centers_g, aes(x=x, y=y), fill="#999999",size=0.5, alpha=0.8)+
    #scale_fill_viridis()+
    map_theme
  p_g_item[[index]]<-p1
  lines<-data.frame(addEllipse(mve$center, mve$cov, col="#800026", p.interval=0.95))
  p2<-ggplot()+
    geom_point(data=points_no_NA[sample(nrow(points_no_NA), 1000), ], 
               aes(x=V_1, y=V_2), color="black", alpha=0.4, size=0.3)+
    geom_path(data=lines, aes(x=X1, y=X2), color="#800026")+
    geom_point(data=centers_g, aes(x=V_1, y=V_2), color="#999999",size=0.5)+
    geom_point(data=centers_n, aes(x=V_1, y=V_2), color="#800026",size=0.2, alpha=0.8)+
    theme_bw()+xlab("PC 1")+ylab("PC 2")
  p_n_item[[index]]<-p2
  
}
if (F){
  p_g_item[[1]]<-p_g_item[[1]]+coord_fixed()
  p_g_item[[2]]
  p_n_item[[1]]
  p_n_item[[2]]
}
p_g_item[[1]]<-p_g_item[[1]]+coord_fixed()
p_g_item[[2]]<-p_g_item[[2]]+coord_fixed()

parr_1<-ggarrange(p_g_item[[1]], p_n_item[[1]], nrow=1, ncol=2, widths=c(0.8, 1))
parr_2<-ggarrange(p_g_item[[2]], p_n_item[[2]], nrow=1, ncol=2, widths=c(0.8, 1))
parr_g_2<-ggarrange(parr_1, parr_2, nrow=2, ncol=1)
#p<-p+xlim(-7.5e6, -2.5e6)+
#  ylim(-6e6, -2e4)
parr_g_full<-ggarrange(p, parr_g_2, nrow=2, ncol=1, heights=c(8, 12))
if (F){
  ggsave(parr_g_full, file="../Figures/center_map_example/center_map_example.png", 
         width=8, height=8)
  ggsave(parr_g_full, file="../Figures/center_map_example/center_map_example.pdf", 
         width=8, height=8)
}

