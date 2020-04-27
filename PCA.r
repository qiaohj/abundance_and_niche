library(raster)
library("stringr")
template<-"../Raster/Bioclim2.0/ECK4/wc2.0_bio_1km_%s_eck4.tif"
i=1
vars<-c(1, 4, 12, 15)
for (i in c(1:19)){
  print(i)
  r<-raster(sprintf(template, str_pad(i, 2, pad = "0")))
  if (i==1){
    p<-data.frame(rasterToPoints(r))
  }else{
    p[,sprintf("V_%d", i)]<-extract(r, p[, c("x", "y")])
  }
}
