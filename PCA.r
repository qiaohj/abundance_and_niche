library(raster)

template<-"../Raster/Bioclim2.0/ECK4/wc2.0_bio_1km_%s_eck4.tif"

if (F){
  i=1
  #vars<-c(1, 4, 12, 15)
  for (i in c(1:19)){
    print(i)
    r<-raster(sprintf(template, str_pad(i, 2, pad = "0")))
    if (i==1){
      p<-data.frame(rasterToPoints(r))
      p<-p[sample(nrow(p), 100000), ]
    }else{
      p[,sprintf("V_%d", i)]<-extract(r, p[, c("x", "y")])
    }
  }
  
  head(p)
  dim(p)
  p_com<-p[complete.cases(p),]
  dim(p_com)
  pca<-prcomp(p_com[, c(3:21)], scale = T)
  saveRDS(pca, "../Models/pca_model.rda")
  pca_pred<-predict(pca, p_com[, c(3:21)])
  head(pca$x)
  i=1
}
library('raster')
library('RStoolbox')
library("stringr")
rasters<-c()
for (i in c(1:19)){
  rasters<-c(rasters, sprintf(template, str_pad(i, 2, pad = "0")))
}
rasters <- stack(rasters)
beginCluster()
#pca1 <- rasterPCA(rasters)
pca_100000 <- rasterPCA(rasters, nSamples = 100000)  # sample 100000 random grid cells
pca_100000_spca <- rasterPCA(rasters, nSamples = 100000, spca=T)  # sample 100000 random grid cells

saveRDS(pca_100000, "../Models/pca.rda")
pca_100000<-readRDS("../Models/pca.rda")
for (i in c(1:19)){
  print(i)
  writeRaster(pca_100000$map[[i]], sprintf("../Raster/PCs/pc%d.tif", i), overwrite=T)
}
endCluster()
