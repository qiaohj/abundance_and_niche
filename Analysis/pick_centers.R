library(dplyr)

setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")

p<-readRDS("../Tables/all_points.rda")
p$N<-0

ratios<-seq(100, 1000, by=(1000-100)/9)
i=1
threshold<-0.80
rep<-1000
p_sub<-NULL

ratio<-ratios[length(ratios)] * 1000
ratio<-500 * 1000
minN<-((ratio/1000 * 2)+1)^2*threshold
N_p<-nrow(p)

if (F){
  p_sub<-readRDS("../Tables/random_points_100.rda")
  box_more_raw<-readRDS("../Tables/box_with_more_dist.rda")
  box_filtered<-box_more_raw%>%dplyr::filter(land_percentile>=0.8)
  xxx<-data.frame(table(box_filtered$index))
  xxx<-xxx[which(xxx$Freq==14),]
  box_filtered<-box_filtered[which(box_filtered$index %in% xxx$Var1),]
  p_sub<-p_sub%>%dplyr::filter((x %in% box_filtered$x)&(y %in% box_filtered$y))
  p_sub$N<-0
  i=1
  for (i in c(1:nrow(p_sub))){
    
    print(i)
    center<-p_sub[i,]
    target<-sprintf("../Tables/P_I/%f_%f.rda", center$x, center$y)
    if (file.exists(target)){
      next()
    }
    saveRDS(NULL, target)
    sub<-p %>% dplyr::filter(dplyr::between(x, center$x-ratio, center$x+ratio) & 
                               dplyr::between(y, center$y-ratio, center$y+ratio))
    N<-nrow(sub)
    if (N>=minN){
      center$N<-N
      saveRDS(center, target)
    }
  }
 
}
p<-p[sample(nrow(p), nrow(p)),]
for (index in c(1:nrow(p))){
  #index <-round(runif(1, 1, N_p))
  print(paste(index, nrow(p), minN))
  center<-p[index, ]
  target<-sprintf("../Tables/P_I/%f_%f.rda", center$x, center$y)
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  
  sub<-p %>% dplyr::filter(dplyr::between(x, center$x-ratio, center$x+ratio) & 
                             dplyr::between(y, center$y-ratio, center$y+ratio))
  N<-nrow(sub)

  print(paste(N, minN))
  if (N>=minN){
    center$N<-N
    saveRDS(center, target)
  }
}

asdf

fs<-list.files("../Tables/P_Ix", full.names = T)
f<-fs[1]
p_sub<-NULL
for (f in fs){
  df<-readRDS(f)
  if (!is.null(df)){
    if (is.null(p_sub)){
      p_sub<-df
    }else{
      p_sub<-bind_rows(p_sub, df)
    }
  }
}
dim(p_sub)
p_sub<-p_sub[sample(nrow(p_sub), 5000),]
dim(p_sub)
saveRDS(p_sub, "../Tables/random_points_5000_0.8.rda")

#saveRDS(p, "../Tables/all_points.rda")
#saveRDS(p_sub, "../Tables/random_points_100.rda")

library(raster)
library(rgdal)
p_sub<-readRDS("../Tables/random_points_5000_0.8.rda")
p_sub[which(abs((p_sub$x+10123184))<10),]

mask<-raster("../Raster/Bioclim2.0/ECK4/mask_1k.tif")
p_shape<-SpatialPointsDataFrame(p_sub[, c("x", "y")], p_sub, 
                                proj4string = CRS(proj4string(mask)))
plot(p_shape)
writeOGR(obj=p_shape, dsn="../Shape/Random_Centers", layer="Random_Centers_5000", 
         driver="ESRI Shapefile", overwrite_layer = T)

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
