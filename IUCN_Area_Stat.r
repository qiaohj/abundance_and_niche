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


