library("raster")
library("stringr")
setwd("/media/huijieqiao/WD12T/Experiments/abundance_and_niche/abundance_and_niche")
i=1
mask<-raster("../Raster/Bioclim2.0/ECK4/mask_1k.tif")
for (i in c(10:19)){
  print(i)
  r<-raster(sprintf("../Raster/Bioclim2.0/RAW/wc2.0_bio_30s_%s.tif", str_pad(i, 2, pad = "0")))
  r_r<-projectRaster(r,
                     crs = crs(mask),
                     res = res(mask))
  writeRaster(r_r, sprintf("../Raster/Bioclim2.0/ECK4/wc2.0_bio_1km_%s_eck4.tif", str_pad(i, 2, pad = "0")), format="GTiff", 
              datatypeCharacter="FLT4S", NAflag=-9999)
}



