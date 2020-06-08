


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

library(raster)
library(rgdal)
p_sub<-readRDS("../Tables/random_points_100.rda")
mask<-raster("../Raster/Bioclim2.0/ECK4/mask_1k.tif")
p_shape<-SpatialPointsDataFrame(p_sub[, c("x", "y")], p_sub, proj4string = CRS(proj4string(mask)))
plot(p_shape)
writeOGR(obj=p_shape, dsn="../Shape/Random_Centers", layer="Random_Centers", driver="ESRI Shapefile")

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
