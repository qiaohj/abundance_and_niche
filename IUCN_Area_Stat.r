setwd("/media/huijieqiao/WD12T/Experiments/abundance_and_niche/abundance_and_niche")

library(dplyr)
library(raster)
library(Hmisc)
library(stringr)
library(ntbox)
library(hypervolume)


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

ratios<-seq(5, 2930, by=(2930-5)/10)
df_sub<-df%>% filter((N>25)&(N<2930^2))

df_sub$cuts<-cut2(df_sub$N, cuts=ratios^2)
table(df_sub$cuts)
#beginCluster()
vars<-c(1, 4, 12, 15)
i=1
j=1
rep=100
k=1
base<-"/media/huijieqiao/WD12T/Experiments/IUCN_FIX/Raster/IUCN_Range_By_Species_eck4"
template<-"../Raster/Bioclim2.0/ECK4/wc2.0_bio_1km_%s_eck4.tif"
r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, str_pad(i, 2, pad = "0")))
  r_env[[as.character(i)]]<-r
}
for (j in c(1:length(unique(df_sub$cuts)))){
  sub<-df_sub[which(df_sub$cuts==unique(df_sub$cuts)[j]),]
  sub<-sub[sample(nrow(sub), rep), ]
  for (k in c(1:nrow(sub))){
    rec<-sub[k,]
    r<-raster(sprintf("%s/%s/%s.tif", base, rec$group, rec$species))
    p<-data.frame(rasterToPoints(r))
    for (i in vars){
      
      
      print(paste("CUTS:", j, length(unique(df_sub$cuts)),
                  " SPECIES:", k, nrow(sub), " ENV:", i, length(vars), sep="/"))
      p[,sprintf("V_%d", i)]<-extract(r_env[[as.character(i)]], p[, c("x", "y")])
    }
  }
}
#endCluster()
env_train<-p[, c(4:7)]
env_train<-env_train[complete.cases(env_train),]
mve<-ellipsoid_omr(env_train, env_train, 4)


mve<-cov_center(data = env_train, mve = T, 
           level = 0.95, vars = 1:ncol(env_train))
