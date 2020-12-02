library(raster)
library(ggplot2)
library(dplyr)
library(MASS)
library(ggpubr)
library(ggpmisc)
library(Rmisc)
library(vegan)
library(ggpmisc)

setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
#center points
p_sub<-readRDS("../Tables/random_points_100.rda")
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
map_background<-"#f5f5f2"
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



ratio = 200
res = 1000
df_box<-readRDS("../Tables/box.rda")

filter_df_box<-df_box%>%dplyr::filter((semi_ra==ratio)&((alt_sd>=1812)|(alt_sd<=14.3)))
filter_df_box<-df_box%>%dplyr::filter((semi_ra==ratio)&((relavent_dist>=2.545176976)|(relavent_dist<=0.006527766)))

range(df_box[which(df_box$semi_ra==ratio), ]$relavent_dist, na.rm = T)

center1<-filter_df_box[1,]
lines1<-data.frame(x=c(center1$x-ratio*res, center1$x+ratio*res, center1$x+ratio*res, center1$x-ratio*res, center1$x-ratio*res),
      y=c(center1$y-ratio*res, center1$y-ratio*res, center1$y+ratio*res, center1$y+ratio*res, center1$y-ratio*res))
center2<-filter_df_box[2,]
lines2<-data.frame(x=c(center2$x-ratio*res, center2$x+ratio*res, center2$x+ratio*res, center2$x-ratio*res, center2$x-ratio*res),
                   y=c(center2$y-ratio*res, center2$y-ratio*res, center2$y+ratio*res, center2$y+ratio*res, center2$y-ratio*res))

p<-ggplot()+geom_tile(data=p_mask, aes(x=x, y=y, fill=v), alpha=0.4)+
  scale_fill_gradientn(colours = rev(terrain.colors(10)))+
  geom_point(data=p_sub, aes(x=x, y=y), color="lightblue",size=0.5)+
  map_theme

p<-p+geom_path(data=lines1, aes(x=x, y=y))+
  geom_point(data=center1, aes(x=x, y=y), color="red", size=0.5)+
  geom_text(data=center1, aes(x=x, y=y, label="(2)"), vjust=2)


p<-p+geom_path(data=lines2, aes(x=x, y=y))+
  geom_point(data=center2, aes(x=x, y=y), color="red", size=0.5)+
  geom_text(data=center2, aes(x=x, y=y, label="(1)"), vjust=2)
p

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
source("addEllipse.R")
for (index in c(1:2)){
  center<-filter_df_box[index,]
  points<-expand.grid(x=seq(from=center$x - ratio*res, to=center$x + ratio * res, by=res),
                      y=seq(from=center$y - ratio*res, to=center$y + ratio * res, by=res))
  
  
  for (j in vars){
    center[, sprintf("PC_%d", j)]<-extract(r_env[[as.character(j)]], center[, c("x", "y")])
    print(paste("ratio:", ratio, i, nrow(centers), "Var:", j, sep="/"))
    points[, sprintf("V_%d", j)]<-extract(r_env[[as.character(j)]], points[, c("x", "y")])
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
  points_no_NA$alt<-extract(r_r, points_no_NA[, c("x", "y")])
  
  
  
  center$index<-index
  center$semi_ra<-ratio
  center$Niche_PC1<-mve$center[1]
  center$Niche_PC2<-mve$center[2]
  lines<-data.frame(x=c(center$x-ratio*res, center$x+ratio*res, center$x+ratio*res, center$x-ratio*res, center$x-ratio*res),
                    y=c(center$y-ratio*res, center$y-ratio*res, center$y+ratio*res, center$y+ratio*res, center$y-ratio*res))
  
  p1<-ggplot()+geom_tile(data=points_no_NA, aes(x=x, y=y, fill=alt), alpha=0.8)+
    scale_fill_gradientn(colours = rev(terrain.colors(10)))+
    geom_path(data=lines, aes(x=x, y=y), color="black")+
    geom_tile(data=centers_n, aes(x=x, y=y), fill="red",size=0.2)+
    geom_tile(data=centers_g, aes(x=x, y=y), fill="grey",size=0.5, alpha=0.8)+
    map_theme
  p_g_item[[index]]<-p1
  lines<-data.frame(addEllipse(mve$center, mve$cov, col="red", p.interval=0.95))
  p2<-ggplot()+
    geom_point(data=points_no_NA, aes(x=V_1, y=V_2), color=map_background, alpha=0.4, size=0.3)+
    geom_path(data=lines, aes(x=X1, y=X2), color="red")+
    geom_point(data=centers_g, aes(x=V_1, y=V_2), color="grey",size=0.5)+
    geom_point(data=centers_n, aes(x=V_1, y=V_2), color="red",size=0.2, alpha=0.8)+
    theme_bw()+xlab("PC 1")+ylab("PC 2")
  p_n_item[[index]]<-p2
  
}
p_g_item[[1]]
p_g_item[[2]]
p_n_item[[1]]
p_n_item[[2]]

box_more<-readRDS("../Tables/box_with_more_dist.rda")

p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=centers_g_2_centers_n_mean, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(aes(x=centers_g_2_centers_n_mean, 
                  y=centers_n_2_centers_g_mean/1000),
              method="lm")+
  stat_poly_eq(aes(x=centers_g_2_centers_n_mean, 
                   y=centers_n_2_centers_g_mean/1000,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  xlim(c(0, 25))+
  xlab("Mean Mahalanobis distance of spatial center to niche center in environmental space")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()
ggsave(p3, file="../Figures/Figure6/Fig.6.2.png", width=8, height=5)

parr_1<-ggarrange(p_g_item[[1]], p_n_item[[1]], nrow=1, ncol=2)
parr_2<-ggarrange(p_g_item[[2]], p_n_item[[2]], nrow=1, ncol=2)
parr_g_2<-ggarrange(parr_1, parr_2, nrow=2, ncol=1, labels=c("(1)", "(2)"))
parr_g_full<-ggarrange(p, parr_g_2, nrow=2, ncol=1, heights=c(10, 12))
ggsave(parr_g_full, file="../Figures/Figure6/Fig.6.1.pdf", width=6, height=8)
ggsave(parr_g_full, file="../Figures/Figure6/Fig.6.1.png", width=6, height=8)



p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(aes(x=alt_sd, 
                  y=centers_n_2_centers_g_mean/1000),
              method="lm")+
  stat_poly_eq(aes(x=alt_sd, 
                   y=centers_n_2_centers_g_mean/1000,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  xlab("Standard deviation of elevation")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()

ggsave(p3, file="../Figures/Figure6/Fig.6.3.png", width=8, height=5)

p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(aes(x=alt_sd, 
                  y=centers_n_2_centers_g_mean/1000,
                  color=factor(semi_ra)),
              method="lm")+
  xlab("Standard deviation of elevation")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()

ggsave(p3, file="../Figures/Figure6/Fig.6.3.2.png", width=8, height=5)

p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(method="lm", 
              aes(x=alt_sd, 
                  y=centers_n_2_centers_g_mean/1000))+
  stat_poly_eq(aes(x=alt_sd, 
                   y=centers_n_2_centers_g_mean/1000,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  xlab("Standard deviation of elevation")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()+
  facet_wrap(~semi_ra, ncol=3, scale="free")

ggsave(p3, file="../Figures/Figure6/Fig.6.3.3.png", width=12, height=12)


p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_g_2_centers_n_mean, 
                 color=factor(semi_ra)), size=0.5)+
  xlab("Standard deviation of elevation")+
  ylab("Mean Mahalanobis distance of spatial center to niche center in environmental space")+
  stat_smooth(aes(x=alt_sd, 
                  y=centers_g_2_centers_n_mean),
              method="lm")+
  stat_poly_eq(aes(x=alt_sd, 
                   y=centers_g_2_centers_n_mean,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  ylim(c(0, 25))+
  labs(color = "Radius (km)")+
  theme_bw()

ggsave(p3, file="../Figures/Figure6/Fig.6.4.png", width=8, height=5)


box_more$jaccard<-box_more$centers_g_in_centers_n/
  (box_more$centers_n_out_centers_g+box_more$centers_g_out_centers_n+box_more$centers_g_in_centers_n)
range(box_more$jaccard)
p3<-ggplot(box_more%>%dplyr::filter(box_more$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=semi_ra, 
                 y=jaccard))+
  stat_smooth(aes(x=semi_ra, 
                  y=jaccard),
              method="lm")+
  stat_poly_eq(aes(x=semi_ra, 
                   y=jaccard,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  xlab("Radius (km)")+
  ylab("Jaccard similarity")+
  theme_bw()
ggsave(p3, file="../Figures/Figure6/Fig.6.5.png", width=8, height=5)


box_more_se<-box_more%>%dplyr::group_by(semi_ra)%>%
  dplyr::summarise(mean_jaccard=mean(jaccard),
                   sd_jaccard=sd(jaccard),
                   CI_jaccard=CI(jaccard)[1]-CI(jaccard)[2])
p3<-ggplot(box_more_se)+
  geom_point(aes(x=factor(semi_ra), 
                 y=mean_jaccard))+
  geom_errorbar(aes(x=factor(semi_ra), 
                y=mean_jaccard,
                ymin=mean_jaccard-CI_jaccard,
                ymax=mean_jaccard+CI_jaccard))+
  xlab("Radius (km)")+
  ylab("Mean Jaccard similarity")+
  theme_bw()
ggsave(p3, file="../Figures/Figure6/Fig.6.6.png", width=8, height=5)
