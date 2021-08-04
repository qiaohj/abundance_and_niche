library(raster)
library(spatstat)
library(rgdal)
library(ntbox)
library(ggplot2)
library(viridis)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")

if (F){
  r<-raster("../Raster/bSTFLx_CONUS_HabMap_2001v1.tif")
  alt<-raster("../Raster/alt_eck4_1km_new.tif")
  win<-readOGR("../Shape/bSTFLx_CONUS_Range_2001v1/bSTFLx_CONUS_Range_2001v1.shp")
  win_eck4 <- spTransform(win, crs(alt))
  
  r_box<-crop(alt, extent(win_eck4))
  r_box<-mask(r_box, win_eck4)
  alt_aea<- projectRaster(r_box, crs = crs(r), res=1000)
  writeRaster(alt_aea, "../Raster/alt_aea.tif")
  
  
  
  pc1<-raster("../Raster/PCs/pc1_example.tif")
  pc1_aea<- projectRaster(pc1, crs = crs(r), res=1000)
  writeRaster(pc1_aea, "../Raster/PCs/pc1_aea.tif")
  
  pc2<-raster("../Raster/PCs/pc2_example.tif")
  pc2_aea<- projectRaster(pc2, crs = crs(r), res=1000)
  writeRaster(pc2_aea, "../Raster/PCs/pc2_aea.tif")
  
  
  pc1<-raster("../Raster/PCs/pc1_aea.tif")
  pc2<-raster("../Raster/PCs/pc2_aea.tif")
  r_eck4_p_table_with_xy<-data.frame(rasterToPoints(pc1))
  r_eck4_p_table_with_xy$pc1<-extract(pc1, r_eck4_p_table_with_xy[, c("x", "y")])
  r_eck4_p_table_with_xy$pc2<-extract(pc2, r_eck4_p_table_with_xy[, c("x", "y")])
  
  r_eck4_p_table_with_xy<-r_eck4_p_table_with_xy[complete.cases(r_eck4_p_table_with_xy),]
  geo_center_x<-mean(r_eck4_p_table_with_xy$x)
  geo_center_y<-mean(r_eck4_p_table_with_xy$y)
  
  mve<-cov_center(data = r_eck4_p_table_with_xy, mve = T, 
                  level = 0.95, vars = 4:5)
  
  niche_center_pc1<-mve$centroid[1]
  niche_center_pc2<-mve$centroid[2]
  
  r_eck4_p_table_with_xy$dist_geo_center<-sqrt((r_eck4_p_table_with_xy$x-geo_center_x)^2+(r_eck4_p_table_with_xy$y-geo_center_y)^2)
  r_eck4_p_table_with_xy$dist_niche_center<-stats::mahalanobis(r_eck4_p_table_with_xy[, c("pc1", "pc2")], center = mve$centroid, 
                                                               cov = mve$covariance)
  
  plot(r_eck4_p_table_with_xy$dist_geo_center, r_eck4_p_table_with_xy$dist_niche_center)
  cor(r_eck4_p_table_with_xy$dist_niche_center, r_eck4_p_table_with_xy$dist_geo_center)
  map_background<-"#f5f5f2"
  map_theme<-theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = map_background, color = NA), 
    panel.background = element_blank(), 
    legend.background = element_rect(fill = map_background, color = NA),
    panel.border = element_blank(),
    legend.position="none"
  )
  centers<-data.frame(x=c(geo_center_x, niche_center_pc1), 
                      y=c(geo_center_y, niche_center_pc2), 
                      type=c("geo_center", "niche_center"))
  saveRDS(centers, "../Tables/example_centers.rda")
  
  geo_center_threshold<-quantile(r_eck4_p_table_with_xy$dist_geo_center, 0.01)
  
  niche_center_threshold<-quantile(r_eck4_p_table_with_xy$dist_niche_center, 0.01)
  r_eck4_p_table_with_xy$is_geo_center<-ifelse(r_eck4_p_table_with_xy$dist_geo_center<geo_center_threshold, T, F)
  r_eck4_p_table_with_xy$is_niche_center<-ifelse(r_eck4_p_table_with_xy$dist_niche_center<niche_center_threshold, T, F)
  terrain.colors<-c("#ccece6", "#99d8c9", "#66c2a4", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f")
  r_eck4_p_table_with_xy$alt<-extract(alt_aea, r_eck4_p_table_with_xy[, c("x", "y")])
  p<-ggplot(r_eck4_p_table_with_xy)+geom_tile(aes(x=x, y=y, fill=alt))+
    geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_geo_center),], aes(x=x, y=y), fill="#0072B2", alpha=1)+
    #geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_density_center),], aes(x=x, y=y), fill="#56B4E9", alpha=0.8)+
    geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_niche_center),], aes(x=x, y=y), fill="#CC79A7", alpha=0.8)+
    scale_fill_gradientn(colours = terrain.colors)+
    #scale_fill_viridis()+
    map_theme
  p
  ggsave(p, filename="../Figures/example_g.png", width = 8, height=8)
  saveRDS(r_eck4_p_table_with_xy, "../Tables/example_eck4_table.rda")
  
}

mask<-raster("../Raster/alt_eck4_1km_new.tif")
r_eck4<-raster("../Shape/bSTFLx_CONUS_HabMap_2001v1/bSTFLx_CONUS_HabMap_2001v1_eck4.tif")
win<-read_sf("../Shape/bSTFLx_CONUS_Range_2001v1/bSTFLx_CONUS_Range_2001v1.shp")

vars<-c(1, 2)
template<-"../Raster/PCs/pc%d.tif"
r_env<-list()
for (i in vars){
  r<-raster(sprintf(template, i))
  r_env[[as.character(i)]]<-r
}


geo_center<-(extent(r_eck4)[1]+extent(r_eck4)[2])/2
geo_center<-(extent(r_eck4)[1]+extent(r_eck4)[2])/2



if (F){
  w  <- as.owin(win)
  p_df<-data.frame(rasterToPoints(r))
  ppoints<-SpatialPointsDataFrame(p_df[, c("x", "y")], p_df)
  ppoints_sf<-st_as_sf(ppoints)
  ppoints_p<-as.ppp(ppoints_sf)
  Q <- quadratcount(ppoints, nx= 6, ny=3)
  
  K2 <- density(r, sigma=1000) # Using a 1000m bandwidth
  plot(K2, main=NULL, las=1)
  contour(K2, add=TRUE)
}

