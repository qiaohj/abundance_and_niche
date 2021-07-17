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
##continue with centers_map_example.R

df_box<-readRDS("../Tables/box_with_more_dist.rda")
df_box<-df_box%>%dplyr::filter(land_percentile>=0.8)
df_box<-df_box%>%dplyr::filter(df_box$semi_ra>10)
df_box$xxx<-df_box$centers_n_in_centers_g/
  (df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)
ratio = 400
res = 1000
filter_df_box_1<-df_box%>%dplyr::filter((semi_ra==ratio)&(xxx>=0.5)&(land_percentile==1))
df_box$dist_ttt<-sqrt((df_box$x-filter_df_box_1[3, "x"])^2+
                        (df_box$y-filter_df_box_1[3, "y"])^2)


hist(df_box$centers_g_2_centers_n_mean)

df_box$xxx<-df_box$centers_n_in_centers_g/
  (df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)
df_box$yyy<-df_box$centers_g_in_centers_n/
  (df_box$centers_g_in_centers_n+df_box$centers_g_out_centers_n)
#plot(df_box$xxx<-df_box$centers_n_in_centers_g, df_box$xxx<-df_box$centers_n_in_centers_n)
#df_box<-df_box%>%dplyr::filter(semi_ra %in% c(30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500))
p<-ggplot(df_box)+geom_density(aes(x=xxx, color=factor(semi_ra)))+theme_bw()+
  xlab("")+facet_wrap(~semi_ra)
binwidth<-0.1
p1<-ggplot(df_box)+
  geom_histogram(aes(x=xxx, fill=factor(width)), bins=20)+
  geom_density(aes(x=xxx, y = ..density..* nrow(df_box) * binwidth, 
                   color=factor(width)), bw = binwidth) +  
  scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  theme_bw()+
  ylab("Count")+
  xlab("Overlap of spatial space center and environmental centroid (square root transform)")+
  labs(color="width", fill="width")
p1
#ggsave(p, filename="../Figures/Figure6/Fig.6.1.1.png", width=7, height=4)

p2<-ggplot(df_box)+
  #geom_histogram(aes(x=yyy, fill=factor(width)), bins=20)+
  geom_density(aes(x=yyy, y = ..density..* nrow(df_box) * binwidth, 
                   color=factor(width)), bw = binwidth) +  
  #scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  scale_x_continuous(breaks=seq(0, 0.9, 0.1), labels = paste(seq(0, 0.9, 0.1)*100, "%"))+
  ylab("Count")+
  theme_bw()+
  xlab("Overlap of spatial space center and environmental centroid")+
  labs(color="width", fill="width")+
  theme(axis.title.y = element_blank())
p2
leg<-get_legend(p1)
p<-ggarrange(p1, p2, nrow=1, ncol=2, common.legend=T, legend.grob=leg, legend = "right")
p
ggsave(p2, filename="../Figures/overlap/overlap.combined.png", width=7, height=5)
ggsave(p2, filename="../Figures/overlap/overlap.combined.pdf", width=7, height=5)

range(df_box$centers_g_2_centers_n_mean, na.rm = T)
p_x<-ggplot(df_box)+
  geom_histogram(aes(x=centers_g_2_centers_n_mean, fill=factor(width)), bins=50)+
  geom_vline(xintercept=stats::qchisq(0.95, 2), linetype=2)+
  geom_vline(xintercept=stats::qchisq(0.01, 2), linetype=2)+
  labs(fill="width", x="Mean Mahalanobis distance of spatial center to niche center in environmental space", y="Count")+
  scale_x_sqrt(breaks=seq(0, 6, 1)^2)+
  theme_bw()

p_x_facet<-ggplot(df_box)+
  geom_histogram(aes(x=centers_g_2_centers_n_mean, fill=factor(width)), bins=50)+
  labs(fill="width", x="Mean Mahalanobis distance of spatial center to niche center in environmental space", y="Count")+
  scale_x_sqrt(breaks=seq(0, 6, 1)^2)+
  theme_bw()+
  facet_wrap(~width)
p_x_facet
ggsave(p_x_facet, filename="../Figures/overlap/distance_in_n_facet.png", width=10, height=7)


p_y<-ggplot(df_box)+
  geom_histogram(aes(x=centers_n_2_centers_g_mean/1000, fill=factor(width)), bins=50)+
  labs(fill="width", x="Mean Euclidean distance of spatial center to niche center in spatial space (KM)", y="Count")+
  scale_x_sqrt(breaks=seq(0, 25, 5)^2)+
  theme_bw()
p_y
p_y_facet<-ggplot(df_box)+
  geom_histogram(aes(x=centers_n_2_centers_g_mean/1000, fill=factor(width)), bins=50)+
  labs(fill="width", x="Mean Euclidean distance of spatial center to niche center in spatial space (KM)", y="Count")+
  scale_x_sqrt(breaks=seq(0, 25, 5)^2)+
  theme_bw()+
  facet_wrap(~width)
p_y_facet
ggsave(p_y_facet, filename="../Figures/overlap/distance_in_g_facet.png", width=10, height=7)

df_box$n2g_scale<-as.vector(df_box$centers_n_2_centers_g_mean/(df_box$width*sqrt(2))/1000)
p_y2<-ggplot(df_box)+
  geom_histogram(aes(x=n2g_scale, fill=factor(width)), bins=50)+
  labs(fill="width", 
       x="scaled mean Euclidean distance of spatial center to niche center in spatial space", y="Count")+
  #scale_x_sqrt(breaks=seq(0, 25, 5)^2)+
  theme_bw()
p_y2_facet<-ggplot(df_box)+
  geom_histogram(aes(x=n2g_scale, fill=factor(width)), bins=50)+
  labs(fill="width", 
       x="scaled mean Euclidean distance of spatial center to niche center in spatial space", y="Count")+
  #scale_x_sqrt(breaks=seq(0, 25, 5)^2)+
  theme_bw()+
  facet_wrap(~width)
p_y2_facet
ggsave(p_y2_facet, filename="../Figures/overlap/distance_in_g_facet_scaled.png", width=10, height=7)


leg<-get_legend(p_x)
pp<-ggarrange(p_x, p_y, nrow=2, ncol=1, common.legend = T, legend = "right", legend.grob = leg)
ggsave(pp, filename="../Figures/overlap/histgram.combined.png", width=7, height=7)
ggsave(pp, filename="../Figures/overlap/histgram.combined.pdf", width=7, height=7)

pp<-ggarrange( p_y, p_y2, p_x, p2, nrow=2, ncol=2, 
               common.legend = T, legend = "right", legend.grob = leg)
ggsave(pp, filename="../Figures/overlap/histgram.combined_all.png", width=14, height=7)
ggsave(pp, filename="../Figures/overlap/histgram.combined_all.pdf", width=14, height=7)


df_box$cuts<-">10% and <50%"
df_box[which(df_box$yyy<=0.1), "cuts"]<-"<=10%"
df_box[which(df_box$yyy>=0.5), "cuts"]<-">=50%"
table(df_box$cuts)
df_box$yyy_count<-floor(df_box$yyy * 10)
data.frame(table(df_box$yyy_count))

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
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()
p3
ggsave(p3, file="../Figures/overlap/lm_distance.png", width=8, height=5)
ggsave(p3, file="../Figures/overlap/lm_distance.pdf", width=8, height=5)


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
df_box$n2g_scale<-as.vector(df_box$centers_n_2_centers_g_mean/(df_box$width*sqrt(2))/1000)
#df_box$n2g_scale<-scale(df_box$centers_n_2_centers_g_mean)
df_box$m_scale<-as.vector(base::scale(df_box$M_V1_ovserved, center=T, scale=T))
df_box$width_scale<-as.vector(scale(df_box$width))

lm_m1<-lm(data=df_box, n2g_scale ~ m_scale+width_scale)
summary(lm_m1)

lm_m1_2<-lm(data=df_box, n2g_scale ~ m_scale)
summary(lm_m1_2)

df_box$g2n_scale<-df_box$centers_g_2_centers_n_mean/(df_box$width*sqrt(2))
df_box$g2n_scale<-scale(df_box$centers_g_2_centers_n_mean)
lm_m2<-lm(data=df_box, g2n_scale ~ m_scale+width_scale)
summary(lm_m2)

cor(df_box$width, df_box$M_V1_ovserved)
info<-lm_eqn(lm_m2)
grid.lines = 30
m_scale.pred <- seq(min(df_box$m_scale), 
                          max(df_box$m_scale), 
                          length.out = grid.lines)
width_scale.pred <- seq(min(df_box$width_scale), 
                  max(df_box$width_scale), 
                  length.out = grid.lines)
xy <- expand.grid(m_scale = m_scale.pred, 
                  width_scale = width_scale.pred)
dist.pred <- matrix(predict(lm_m2, newdata = xy), 
                    nrow = grid.lines, ncol = grid.lines)

library("plot3D")

png(filename="../Figures/overlap/lm_n2g.png", width=1000, height=1000)

item<-df_box[which(df_box$g2n_scale<=5),]
scatter3D(item$m_scale, 
          item$width_scale,
          item$g2n_scale, pch = 18, cex = 0.5, 
          theta = 60, phi = 20, ticktype = "detailed",
          xlab = "scaled Moran’s I of PC1", ylab = "scaled width", zlab = "Center E to Center G",
          surf = list(x = m_scale.pred, y = width_scale.pred, z = dist.pred,  
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
info<-lm_eqn(lm_m1)
grid.lines = 30

dist.pred <- matrix(predict(lm_m1, newdata = xy), 
                    nrow = grid.lines, ncol = grid.lines)

library(plot3D)

png(filename="../Figures/overlap/lm_g2n.png", width=1000, height=1000)

scatter3D(df_box_1$m_scale, 
          df_box_1$width_scale,
          df_box_1$n2g_scale, pch = 18, cex = 1, 
          theta = 60, phi = 20, ticktype = "detailed",
          xlab = "scaled Moran’s I of PC1", ylab = "scaled width", zlab = "scaled Center G to Center N",
          surf = list(x = m_scale.pred, y = width_scale.pred, z = dist.pred,  
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
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()+
  facet_wrap(~semi_ra, scale="free")
ggsave(p4, file="../Figures/overlap/lm.pc1.png", width=15, height=12)



p5<-ggplot(df_box)+
  geom_point(aes(x=M_V2_ovserved, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()+
  facet_wrap(~semi_ra, scale="free")
ggsave(p5, file="../Figures/overlap/lm.pc2.png", width=15, height=12)
df_box[which(df_box$M_V1_ovserved==max(df_box$M_V1_ovserved)),]




p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()

ggsave(p3, file="../Figures/overlap/elevation_lm.png", width=8, height=5)

p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(width)), size=0.5)+
  stat_smooth(aes(x=alt_sd, 
                  y=centers_n_2_centers_g_mean/1000,
                  color=factor(width)),
              method="lm")+
  xlab("Standard deviation of elevation")+
  ylab("Mean Euclidean distance of spatial center to niche center in spatial space (KM)")+
  labs(color = "width")+
  theme_bw()

ggsave(p3, file="../Figures/overlap/elevation_lm_group.png", width=8, height=5)

p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_n_2_centers_g_mean/1000, 
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()+
  facet_wrap(~semi_ra, ncol=3, scale="free")

ggsave(p3, file="../Figures/overlap/elevation_lm_part.png", width=12, height=12)


p3<-ggplot(df_box%>%dplyr::filter(df_box$centers_g_2_centers_n_mean<25))+
  geom_point(aes(x=alt_sd, 
                 y=centers_g_2_centers_n_mean, 
                 color=factor(width)), size=0.5)+
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
  labs(color = "width")+
  theme_bw()

ggsave(p3, file="../Figures/overlap/elevation_g2n_group.png", width=8, height=5)


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
