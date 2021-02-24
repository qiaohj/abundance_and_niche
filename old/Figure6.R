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



ratio = 400
res = 1000
df_box<-readRDS("../Tables/box_with_more_dist.rda")
df_box<-df_box%>%dplyr::filter(land_percentile>=0.8)
df_box<-df_box%>%dplyr::filter(df_box$semi_ra>10)
df_box$xxx<-df_box$centers_n_in_centers_g/(df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)



filter_df_box_1<-df_box%>%dplyr::filter((semi_ra==ratio)&(xxx>=0.5)&(land_percentile==1))
df_box$dist_ttt<-sqrt((df_box$x-filter_df_box_1[3, "x"])^2+
                        (df_box$y-filter_df_box_1[3, "y"])^2)

filter_df_box_2<-df_box%>%dplyr::filter((semi_ra==ratio)&((xxx<=0.05)&(dist_ttt<=2e6)))
filter_df_box<-rbind(filter_df_box_1[3,], filter_df_box_2[7,])
hist(filter_df_box_2$dist_ttt)
hist(df_box$xxx)

#View(df_box[which(df_box$index==906),])

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

ggsave(p, file="../Figures/Figure6/Fig.6.0.pdf", width=6, height=4)
ggsave(p, file="../Figures/Figure6/Fig.6.0.png", width=6, height=4)


p<-p+geom_path(data=lines1, aes(x=x, y=y))+
  geom_point(data=center1, aes(x=x, y=y), color="red", size=0.5)+
  geom_text(data=center1, aes(x=x, y=y, label="(1)"), vjust=2)


p<-p+geom_path(data=lines2, aes(x=x, y=y))+
  geom_point(data=center2, aes(x=x, y=y), color="red", size=0.5)+
  geom_text(data=center2, aes(x=x, y=y, label="(2)"), vjust=2)

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

parr_1<-ggarrange(p_g_item[[1]], p_n_item[[1]], nrow=1, ncol=2)
parr_2<-ggarrange(p_g_item[[2]], p_n_item[[2]], nrow=1, ncol=2)
parr_g_2<-ggarrange(parr_1, parr_2, nrow=2, ncol=1, labels=c("(1)", "(2)"))
#p<-p+xlim(-7.5e6, -2.5e6)+
#  ylim(-6e6, -2e4)
parr_g_full<-ggarrange(p, parr_g_2, nrow=2, ncol=1, heights=c(10, 12))
ggsave(parr_g_full, file="../Figures/Figure6/Fig.6.1.pdf", width=6, height=8)
ggsave(parr_g_full, file="../Figures/Figure6/Fig.6.1.png", width=6, height=8)



hist(df_box$centers_g_2_centers_n_mean)

df_box$xxx<-df_box$centers_n_in_centers_g/(df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)
df_box$yyy<-df_box$centers_g_in_centers_n/(df_box$centers_g_in_centers_n
                                               +df_box$centers_g_out_centers_n)


p<-ggplot(df_box)+geom_density(aes(x=xxx, color=factor(semi_ra)))+theme_bw()+
  xlab("")+facet_wrap(~semi_ra)

p1<-ggplot(df_box)+
  #geom_histogram(aes(x=xxx, y=..density.., fill=factor(semi_ra)), bins=10)+
  geom_density(aes(x=xxx, y = ..count.., color=factor(semi_ra)), bw = 0.02) +  
  scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  theme_bw()+
  ylab("Count")+
  xlab("% of overlap in spatial space")+
  labs(color="Width")
p1
#ggsave(p, filename="../Figures/Figure6/Fig.6.1.1.png", width=7, height=4)

p2<-ggplot(df_box)+
  #geom_histogram(aes(x=xxx, y=..density.., fill=factor(semi_ra)), bins=10)+
  geom_density(aes(x=yyy, y = ..count.., color=factor(semi_ra)), bw = 0.02) +  
  scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  ylab("Count")+
  theme_bw()+
  xlab("% of overlap in environmental space")+
  labs(color="Width")+
  theme(axis.title.y = element_blank())
p2
leg<-get_legend(p1)
p<-ggarrange(p1, p2, nrow=1, ncol=2, common.legend=T, legend.grob=leg, legend = "right")
p
ggsave(p, filename="../Figures/Figure6/Fig.6.1.combined.png", width=9, height=4)
#ggsave(p, filename="../Figures/Figure6/Fig.6.1.2.png", width=7, height=4)

p

ggplot(df_box)+geom_histogram(aes(x=yyy))+
  facet_wrap(~semi_ra, scale="free")

cor(df_box$xxx, df_box$yyy)

ggplot(df_box)+geom_histogram(aes(x=centers_n_2_centers_g_mean))+
  facet_wrap(~semi_ra, scale="free")

ggplot(df_box)+geom_histogram(aes(x=centers_g_2_centers_n_mean))+
  facet_wrap(~semi_ra, scale="free")

df_box_f<-df_box%>%dplyr::filter(index %in% df_box[which(df_box$semi_ra==500), "index"])
table(df_box_raw$semi_ra)
table(df_box$semi_ra)
table(df_box_f$semi_ra)

p3<-ggplot(df_box)+
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


lm_eqn <- function(m){
  intercept<-as.numeric(coef(m)[1])
  a<-as.numeric(coef(m)[2])
  b<-as.numeric(coef(m)[3])
  f <- summary(m)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  eq <- substitute(italic(temp) == a + b %.% italic(alt)+ c %.% italic(y)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  return(data.frame(intercept=intercept,
                    a=a,
                    b=b,
                    p=p,
                    r=summary(m)$r.squared)
  )
  
}
df_box$n2g_scale<-df_box$centers_n_2_centers_g_mean/(df_box$width*sqrt(2))
lm_m<-lm(data=df_box, n2g_scale ~ M_V1_ovserved+width)
summary(lm_m)
df_box$g2n_scale<-df_box$centers_g_2_centers_n_mean/(df_box$width*sqrt(2))
lm_m<-lm(data=df_box, g2n_scale ~ M_V1_ovserved+width)
summary(lm_m)

cor(df_box$width, df_box$M_V1_ovserved)
summary(lm_m)
info<-lm_eqn(lm_m)
grid.lines = 30
M_V1_ovserved.pred <- seq(min(df_box$M_V1_ovserved), 
                          max(df_box$M_V1_ovserved), 
                          length.out = grid.lines)
width.pred <- seq(min(df_box$width), 
                  max(df_box$width), 
                  length.out = grid.lines)
xy <- expand.grid(M_V1_ovserved = M_V1_ovserved.pred, 
                  width = width.pred)
dist.pred <- matrix(predict(lm_m, newdata = xy), 
                    nrow = grid.lines, ncol = grid.lines)

library("plot3D")

png(filename="../Figures/Figure6/lm_n2g.png", width=1000, height=1000)

scatter3D(df_box$M_V1_ovserved, 
          df_box$width,
          df_box$n2g_scale, pch = 18, cex = 0.5, 
          theta = 60, phi = 20, ticktype = "detailed",
          xlab = "Moran’s I of PC1", ylab = "Width", zlab = "Center E to Center G",
          surf = list(x = M_V1_ovserved.pred, y = width.pred, z = dist.pred,  
                      facets = NA), 
          main = sprintf("Intercept=%.3f, a=%.3f, b=%.3f, r2=%.3f, p-value=%.3f",
                         info$intercept,
                         info$a,
                         info$b,
                         info$r,
                         info$p))

dev.off()

#df_box_1<-df_box%>%dplyr::filter(centers_g_2_centers_n_mean<=qchisq)
df_box_1<-df_box
lm_m2<-lm(data=df_box, log(centers_g_2_centers_n_mean) ~ M_V1_ovserved)
summary(lm_m2)

lm_m<-lm(data=df_box, log(centers_g_2_centers_n_mean) ~ M_V1_ovserved+width)
summary(lm_m)
info<-lm_eqn(lm_m)
grid.lines = 30
M_V1_ovserved.pred <- seq(min(df_box_1$M_V1_ovserved), 
                          max(df_box_1$M_V1_ovserved), 
                          length.out = grid.lines)
width.pred <- seq(min(df_box_1$width), 
                  max(df_box_1$width), 
                  length.out = grid.lines)
xy <- expand.grid(M_V1_ovserved = M_V1_ovserved.pred, 
                  width = width.pred)
dist.pred <- matrix(predict(lm_m, newdata = xy), 
                    nrow = grid.lines, ncol = grid.lines)



png(filename="../Figures/Figure6/lm_g2n.png", width=1000, height=1000)

scatter3D(df_box_1$M_V1_ovserved, 
          df_box_1$width,
          df_box_1$centers_g_2_centers_n_mean, pch = 18, cex = 0.5, 
          theta = 60, phi = 20, ticktype = "detailed",
          xlab = "Moran’s I of PC1", ylab = "Width", zlab = "Center G to Center N",
          surf = list(x = M_V1_ovserved.pred, y = width.pred, z = dist.pred,  
                      facets = NA), 
          main = sprintf("Intercept=%.3f, a=%.3f, b=%.3f, r2=%.3f, p=%.3f",
                         info$intercept,
                         info$a,
                         info$b,
                         info$r,
                         info$p))

dev.off()

ggplot(df_box%>%dplyr::filter(semi_ra==200))+
  geom_point(aes(x=n2g_scale, y=M_V1_ovserved))

ggplot(df_box%>%dplyr::filter(semi_ra==200))+
  geom_point(aes(x=centers_g_2_centers_n_mean, y=M_V1_ovserved))

ggplot(df_box)+
  geom_point(aes(x=semi_ra, y=M_V1_ovserved))

cor(df_box$semi_ra, df_box$M_V1_ovserved)

hist(df_box$land_percentile)
p4<-ggplot(df_box)+
  geom_point(aes(x=M_V1_ovserved, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(aes(x=M_V1_ovserved, 
                  y=centers_n_2_centers_g_mean/1000),
              method="lm")+
  stat_poly_eq(aes(x=M_V1_ovserved, 
                   y=centers_n_2_centers_g_mean/1000,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  
  xlab("Moran's I autocorrelation index of PC1")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()+
  facet_wrap(~semi_ra, scale="free")
ggsave(p4, file="../Figures/Figure6/Fig.6.2.pc1.png", width=15, height=12)



p5<-ggplot(df_box)+
  geom_point(aes(x=M_V2_ovserved, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(semi_ra)), size=0.5)+
  stat_smooth(aes(x=M_V2_ovserved, 
                  y=centers_n_2_centers_g_mean/1000),
              method="lm")+
  stat_poly_eq(aes(x=M_V2_ovserved, 
                   y=centers_n_2_centers_g_mean/1000,
                   label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label.., 
                                 sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE)+
  
  xlab("Moran's I autocorrelation index of PC2")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "Radius (km)")+
  theme_bw()+
  facet_wrap(~semi_ra, scale="free")
ggsave(p5, file="../Figures/Figure6/Fig.6.2.pc2.png", width=15, height=12)
df_box[which(df_box$M_V1_ovserved==max(df_box$M_V1_ovserved)),]




p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
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

p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
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

p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
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


p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
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


df_box$jaccard<-df_box$centers_g_in_centers_n/
  (df_box$centers_n_out_centers_g+df_box$centers_g_out_centers_n+df_box$centers_g_in_centers_n)
range(df_box$jaccard)
p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
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


df_box_se<-df_box%>%dplyr::group_by(semi_ra)%>%
  dplyr::summarise(mean_jaccard=mean(jaccard),
                   sd_jaccard=sd(jaccard),
                   CI_jaccard=CI(jaccard)[1]-CI(jaccard)[2])
p3<-ggplot(df_box_se)+
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
