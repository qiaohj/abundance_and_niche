library(raster)
library(spatstat)
library(rgdal)
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")

if (F){
  mask<-raster("../Raster/alt_eck4_1km_new.tif")
  r<-raster("../Shape/bSTFLx_CONUS_HabMap_2001v1/bSTFLx_CONUS_HabMap_2001v1.tif")
  
  r_eck4<- projectRaster(r, crs = crs(mask), res=30)
  writeRaster(r_eck4, "../Shape/bSTFLx_CONUS_HabMap_2001v1/bSTFLx_CONUS_HabMap_2001v1_eck4.tif")
  
  win<-readOGR("../Shape/bSTFLx_CONUS_Range_2001v1/bSTFLx_CONUS_Range_2001v1.shp")
  win_eck4 <- spTransform(win,
                                crs(r_eck4))
  
  r<-raster("../Raster/PCs/pc1.tif")
  plot(r)
  plot(win_eck4, add=T)
  r_box<-crop(r, extent(win_eck4))
  r_box<-mask(r_box, win_eck4)
  writeRaster(r_box, "../Raster/PCs/pc1_example.tif")
  
  
  r<-raster("../Raster/PCs/pc2.tif")
  plot(r)
  plot(win_eck4, add=T)
  r_box<-crop(r, extent(win_eck4))
  r_box<-mask(r_box, win_eck4)
  writeRaster(r_box, "../Raster/PCs/pc2_example.tif")
  
  mask<-raster("../Raster/alt_eck4_10km.tif")
  ppp<-data.frame(rasterToPoints(mask))
  values(mask)[!is.na(values(mask))]<-c(1:nrow(ppp))
  r_eck4<-raster("../Shape/bSTFLx_CONUS_HabMap_2001v1/bSTFLx_CONUS_HabMap_2001v1_eck4.tif")
  r_eck4_p<-data.frame(rasterToPoints(r_eck4))
  r_eck4_p$index<-extract(mask, r_eck4_p[, c("x", "y")])
  saveRDS(r_eck4_p, "../Tables/example_eck4_p.rda")
  pp_index<-data.frame(rasterToPoints(mask))
  colnames(pp_index)[3]<-"index"
  r_eck4_p_table<-data.frame(table(r_eck4_p$index))
  colnames(r_eck4_p_table)<-c("index", "Freq")
  r_eck4_p_table_with_xy<-merge(r_eck4_p_table, pp_index, by="index")
  points<-SpatialPointsDataFrame(r_eck4_p_table_with_xy[, c("x", "y")], r_eck4_p_table_with_xy,
                                 proj4string=crs(win_eck4))
  is_over<-over(points, win_eck4)
  r_eck4_p_table_with_xy$over<-is_over$SeasonCode
  hist(r_eck4_p_table_with_xy$Freq)
  r_eck4_p_table_with_xy<-r_eck4_p_table_with_xy[which(!is.na(r_eck4_p_table_with_xy$over)),]
  ggplot(r_eck4_p_table_with_xy)+
    geom_tile(aes(x=x, y=y, fill=Freq))+theme_bw()
  colnames(r_eck4_p_table_with_xy)[2]<-"density"
  
  pc1<-raster("../Raster/PCs/pc1.tif")
  pc2<-raster("../Raster/PCs/pc2.tif")
  
  r_eck4_p_table_with_xy$pc1<-extract(pc1, r_eck4_p_table_with_xy[, c("x", "y")])
  r_eck4_p_table_with_xy$pc2<-extract(pc2, r_eck4_p_table_with_xy[, c("x", "y")])
  
  r_eck4_p_table_with_xy<-r_eck4_p_table_with_xy[complete.cases(r_eck4_p_table_with_xy),]
  geo_center_x<-mean(r_eck4_p_table_with_xy$x)
  geo_center_y<-mean(r_eck4_p_table_with_xy$y)
  
  density_center_x<-sum(r_eck4_p_table_with_xy$x * r_eck4_p_table_with_xy$density)/sum(r_eck4_p_table_with_xy$density)
  density_center_y<-sum(r_eck4_p_table_with_xy$y * r_eck4_p_table_with_xy$density)/sum(r_eck4_p_table_with_xy$density)
  mve<-cov_center(data = r_eck4_p_table_with_xy, mve = T, 
                  level = 0.95, vars = 6:7)
  
  niche_center_pc1<-mve$centroid[1]
  niche_center_pc2<-mve$centroid[2]
  
  r_eck4_p_table_with_xy$dist_geo_center<-sqrt((r_eck4_p_table_with_xy$x-geo_center_x)^2+(r_eck4_p_table_with_xy$y-geo_center_y)^2)
  r_eck4_p_table_with_xy$dist_density_center<-sqrt((r_eck4_p_table_with_xy$x-density_center_x)^2+(r_eck4_p_table_with_xy$y-density_center_y)^2)
  r_eck4_p_table_with_xy$dist_niche_center<-stats::mahalanobis(r_eck4_p_table_with_xy[, c("pc1", "pc2")], center = mve$centroid, 
                                                               cov = mve$covariance)
  
  plot(r_eck4_p_table_with_xy$dist_geo_center, r_eck4_p_table_with_xy$dist_niche_center)
  plot(r_eck4_p_table_with_xy$dist_geo_center, r_eck4_p_table_with_xy$dist_density_center)
  plot(r_eck4_p_table_with_xy$dist_density_center, r_eck4_p_table_with_xy$dist_niche_center)
  cor(r_eck4_p_table_with_xy$dist_niche_center, r_eck4_p_table_with_xy$dist_density_center)
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
  centers<-data.frame(x=c(geo_center_x, density_center_x, niche_center_pc1), 
                      y=c(geo_center_y, density_center_y, niche_center_pc2), 
                      type=c("geo_center", "density_center", "niche_center"))
  saveRDS(centers, "../Tables/example_centers.rda")
  
  geo_center_threshold<-quantile(r_eck4_p_table_with_xy$dist_geo_center, 0.01)
  density_center_threshold<-quantile(r_eck4_p_table_with_xy$dist_density_center, 0.01)
  niche_center_threshold<-quantile(r_eck4_p_table_with_xy$dist_niche_center, 0.01)
  r_eck4_p_table_with_xy$is_geo_center<-ifelse(r_eck4_p_table_with_xy$dist_geo_center<geo_center_threshold, T, F)
  r_eck4_p_table_with_xy$is_density_center<-ifelse(r_eck4_p_table_with_xy$dist_density_center<density_center_threshold, T, F)
  r_eck4_p_table_with_xy$is_niche_center<-ifelse(r_eck4_p_table_with_xy$dist_niche_center<niche_center_threshold, T, F)
  
  p<-ggplot(r_eck4_p_table_with_xy)+geom_tile(aes(x=x, y=y, fill=density))+
    geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_geo_center),], aes(x=x, y=y), fill="#0072B2", alpha=1)+
    geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_density_center),], aes(x=x, y=y), fill="#56B4E9", alpha=0.8)+
    geom_tile(data=r_eck4_p_table_with_xy[which(r_eck4_p_table_with_xy$is_niche_center),], aes(x=x, y=y), fill="#CC79A7", alpha=0.6)+
    scale_fill_gradient(low="#999999", high="#E69F00")+
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

