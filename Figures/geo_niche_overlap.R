
##continue with centers_map_example.R

hist(df_box$centers_g_2_centers_n_mean)

df_box$xxx<-df_box$centers_n_in_centers_g/
  (df_box$centers_n_in_centers_g+df_box$centers_n_out_centers_g)
df_box$yyy<-df_box$centers_g_in_centers_n/
  (df_box$centers_g_in_centers_n+df_box$centers_g_out_centers_n)
#plot(df_box$xxx<-df_box$centers_n_in_centers_g, df_box$xxx<-df_box$centers_n_in_centers_n)

p<-ggplot(df_box)+geom_density(aes(x=xxx, color=factor(semi_ra)))+theme_bw()+
  xlab("")+facet_wrap(~semi_ra)
binwidth<-0.05
p1<-ggplot(df_box)+
  geom_histogram(aes(x=xxx, fill=factor(semi_ra)), bins=20)+
  geom_density(aes(x=xxx, y = ..density..* nrow(df_box) * binwidth, 
                   color=factor(semi_ra)), bw = binwidth) +  
  scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  theme_bw()+
  ylab("Count")+
  xlab("Overlap of spatial space center and environmental centroid (square root transform)")+
  labs(color="radius", fill="radius")
p1
#ggsave(p, filename="../Figures/Figure6/Fig.6.1.1.png", width=7, height=4)

p2<-ggplot(df_box)+
  geom_histogram(aes(x=yyy, fill=factor(semi_ra)), bins=20)+
  geom_density(aes(x=yyy, y = ..density..* nrow(df_box) * binwidth, 
                   color=factor(semi_ra)), bw = binwidth) +  
  #scale_x_sqrt(breaks=seq(0, 0.9, 0.1)^2)+
  ylab("Count")+
  theme_bw()+
  xlab("Overlap of spatial space center and environmental centroid")+
  labs(color="radius", fill="radius")+
  theme(axis.title.y = element_blank())
p2
leg<-get_legend(p1)
p<-ggarrange(p1, p2, nrow=1, ncol=2, common.legend=T, legend.grob=leg, legend = "right")
p
ggsave(p2, filename="../Figures/overlap/overlap.combined.png", width=7, height=5)
ggsave(p2, filename="../Figures/overlap/overlap.combined.pdf", width=7, height=5)

p
df_box$cuts<-">10% and <50%"
df_box[which(df_box$yyy<=0.1), "cuts"]<-"<=10%"
df_box[which(df_box$yyy>=0.5), "cuts"]<-">=50%"
table(df_box$cuts)
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
