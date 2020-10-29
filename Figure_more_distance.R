library(dplyr)
library(ggplot2)
library(raster)
library(Rmisc)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
df_box<-readRDS("../Tables/box_with_more_dist.rda")
df_box<-df_box[which(df_box$mh_dist_g_center_mean<100),]
p<-ggplot(df_box, aes(x=g_dist_g_center_mean, y=mh_dist_g_center_mean, color=semi_ra))+
  geom_point()+theme_bw()+ylim(0, 100)
ggsave(p, filename = "../Figures/g_mh_sc_p.png")

range(df_box$mh_dist_g_center_mean, na.rm = T)
df_box_se<-df_box%>%dplyr::group_by(semi_ra)%>%
  dplyr::summarise(mean_mh_dist_g_center_mean=mean(mh_dist_g_center_mean, na.rm=T),
                   sd_mh_dist_g_center_mean=sd(mh_dist_g_center_mean, na.rm=T),
                   CI_mh_dist_g_center_mean=CI(mh_dist_g_center_mean)[2]-CI(mh_dist_g_center_mean)[3],
                   mean_g_dist_g_center_mean=mean(g_dist_g_center_mean, na.rm=T),
                   sd_g_dist_g_center_mean=sd(g_dist_g_center_mean, na.rm=T),
                   CI_g_dist_g_center_mean=CI(g_dist_g_center_mean)[2]-CI(g_dist_g_center_mean)[3],)
p<-ggplot(df_box_se, aes(x=semi_ra, y=mean_g_dist_g_center_mean))+
  #geom_errorbar(aes(ymin=mean_g_dist_g_center_mean-CI_g_dist_g_center_mean, 
  #                  ymax=mean_g_dist_g_center_mean+CI_g_dist_g_center_mean), width=2,
  #              position=position_dodge(.9)) +
  geom_point()+geom_line()+theme_bw()
ggsave(p, filename = "../Figures/g_semi_ra.png")

p<-ggplot(df_box_se, aes(x=semi_ra, y=mean_mh_dist_g_center_mean))+
  #geom_errorbar(aes(ymin=mean_mh_dist_g_center_mean-CI_mh_dist_g_center_mean, 
  #                  ymax=mean_mh_dist_g_center_mean+CI_mh_dist_g_center_mean), width=2,
  #              position=position_dodge(.9)) +
  geom_point()+geom_line()+theme_bw()
ggsave(p, filename = "../Figures/mh_semi_ra.png")

p<-ggplot(df_box, aes(x=g_dist_g_center_mean, y=mh_dist_g_center_mean, color=semi_ra))+
  geom_point()+theme_bw()+ylim(0, 100)+facet_wrap(~semi_ra, ncol=3)
ggsave(p, filename = "../Figures/g_mh_sc_p_wrap.png")

df_box_se$mean_mh_dist_g_center_mean
df_box$mh_dist_g_center_mean

t_test_df<-NULL
ra1=10
ra2=100
for (ra1 in unique(df_box$semi_ra)){
  for (ra2 in unique(df_box$semi_ra)){
    t_re<-t.test(df_box[which(df_box$semi_ra==ra1), "mh_dist_g_center_mean"],
                 df_box[which(df_box$semi_ra==ra2), "mh_dist_g_center_mean"])
    item<-data.frame(ra1=ra1, ra2=ra2, p.value=t_re$p.value)
    if (is.null(t_test_df)){
      t_test_df<-item
    }else{
      t_test_df<-rbind(t_test_df, item)
    }
  }
}
p<-ggplot(t_test_df, aes(factor(ra1), factor(ra2), fill= p.value)) + 
  geom_tile()
ggsave(p, filename = "../Figures/g_mh_t.test.png")

library(pheatmap)
ma<-matrix(t_test_df$p.value, nrow=length(unique(df_box$semi_ra)))
colnames(ma)<-unique(df_box$semi_ra)
rownames(ma)<-unique(df_box$semi_ra)

pheatmap(ma, display_numbers = T, legend=F, cluster_cols=F, cluster_row=F,
         color=colorRampPalette(c("white", "red"))(100),
         filename = "../Figures/g_mh_t.test_with_value.png")
