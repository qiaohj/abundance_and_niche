setwd("/media/huijieqiao/WD12T/Experiments/abundance_and_niche/abundance_and_niche")
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
library(dplyr)
library(raster)
library(Hmisc)
library(stringr)
library(ntbox)
library(hypervolume) 

if (F){
  groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")
  
  df<-data.frame()
  g<-groups[1]
  min_max<-NULL
  for (g in groups){
    print(g)
    df_item<-readRDS(sprintf("../Tables/Distribution_Area/%s.rda", g))
    cuts<-quantile(df_item$N, seq(from=0, to=1, by=0.1))
    df_item$cuts<-as.numeric(cut2(df_item$N, cuts=cuts))
    
    range<-range(df_item[which(!(as.numeric(df_item$cuts) %in% c(1, length(cuts)-1))), "N"])
    if (is.null(min_max)){
      min_max<-range
      df<-df_item
    }else{
      min_max[1]<-if(min_max[1]>range[1]) range[1] else min_max[1]
      min_max[2]<-if(min_max[2]<range[2]) range[2] else min_max[2]
      df<-bind_rows(df, df_item)
    }
    
  }
  
  df_sub<-df%>% filter((N>min_max[1])&(N<min_max[2]))
  
  ggplot(df_sub) + geom_density(aes(x = N, group=factor(group), color=factor(group))) +
    scale_x_log10()
  mask<-raster("../Raster/Bioclim2.0/ECK4/mask_1k.tif")
  p<-as_tibble(rasterToPoints(mask))
  
  
  
  saveRDS(p, "../Tables/all_points.rda")
  
  df_sub<-df%>% filter((N>25)&(N<2930^2))
  
  
  df_sub$cuts<-cut2(df_sub$N, cuts=ratios^2)
  table(df_sub$cuts)
  #beginCluster()
}





p<-readRDS("../Tables/all_points.rda")
p$N<-0

ratios<-seq(100, 1000, by=(1000-100)/9)
i=1
threshold<-0.99
rep<-1000
p_sub<-NULL
ratio<-ratios[length(ratios)] * 1000
ratio<-100 * 1000
while(i<=rep){
  index <-round(runif(1, 1, nrow(p)))
  center<-p[index, ]
  N<-center$N
  if (N==0){
    sub<-p %>% filter(dplyr::between(x, center$x-ratio, center$x+ratio) & 
                        dplyr::between(y, center$y-ratio, center$y+ratio))
    N<-nrow(sub)
    center$N<-N
    p[index, "N"]<-center$N
  }
  minN<-((ratio/1000 * 2)+1)^2*threshold
  print(paste(i, rep, N, minN))
  if (N>=minN){
    if (is.null(p_sub)){
      p_sub<-center
    }else{
      p_sub<-bind_rows(center, p_sub)
    }
    i <- i+1
  }
}
saveRDS(p, "../Tables/all_points.rda")
saveRDS(p_sub, "../Tables/random_points_100.rda")
plot(p_sub$x, p_sub$y)

vars<-c(1, 4, 12, 15)
i=1
j=1

k=1
template<-"../Raster/Bioclim2.0/ECK4/wc2.0_bio_1km_%s_eck4.tif"
r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, str_pad(i, 2, pad = "0")))
  r_env[[as.character(i)]]<-r
}
for (ratio in ratios/2){
}

env_train<-p[, c(4:7)]
env_train<-env_train[complete.cases(env_train),]
mve<-ellipsoid_omr(env_train, env_train, 4)


mve<-cov_center(data = env_train, mve = T, 
                level = 0.95, vars = 1:ncol(env_train))
