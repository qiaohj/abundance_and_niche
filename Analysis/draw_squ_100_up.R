
setwd("/Volumes/Disk2/Experiments/abundance_and_niche/abundance_and_niche")
library(dplyr)
library(raster)
library(Hmisc)
library(stringr)
library(ntbox)
library(hypervolume) 
library(fields)
library(ggfortify)
library(GWmodel)
library(ape)
centers<-readRDS("../Tables/random_points_100.rda")
i=1
for (i in c(1:nrow(centers))){
  f<-sprintf("../Tables/box_%d.rda", i)
  if (file.exists(f)){
    next()
  }
  result<-NULL
  saveRDS(result, f)
  for (ratio in ratios){
    if (!is.null(result)){
      check<-result %>% filter((index==i)&(semi_ra==ratio))
      if (nrow(check)>0){
        print("skip")
        next()
      }
    }
    
    print(paste(ratio, i, sep="/"))
    folder<-sprintf(base, ratio, i)
    dir.create(folder, showWarnings = F)
    center<-centers[i,]
    points<-expand.grid(x=seq(from=center$x - ratio*res, to=center$x + ratio * res, by=res),
                        y=seq(from=center$y - ratio*res, to=center$y + ratio * res, by=res))
    
    size<-nrow(points)
    for (j in vars){
      center[, sprintf("PC_%d", j)]<-extract(r_env[[as.character(j)]], center[, c("x", "y")])
      print(paste("ratio:", ratio, i, nrow(centers), "Var:", j, sep="/"))
      points[, sprintf("V_%d", j)]<-extract(r_env[[as.character(j)]], points[, c("x", "y")])
    }
    saveRDS(points, sprintf("%s/points.rda", folder))
    points_no_NA<-points[complete.cases(points),]
    size_continent<-nrow(points_no_NA)
    mve<-cov_center(data = points_no_NA, mve = T, 
                    level = 0.95, vars = 3:4)
    saveRDS(mve, sprintf("%s/mve.rda", folder))
    center$index<-i
    center$semi_ra<-ratio
    center$Niche_PC1<-mve$centroid[1]
    center$Niche_PC2<-mve$centroid[2]
    center$Niche_Sp_distance<-sqrt((center$PC_1-center$Niche_PC1)^2+(center$PC_2-center$Niche_PC2)^2)
    center$r_a<-mve$SemiAxis_length[1]
    center$r_b<-mve$SemiAxis_length[2]
    center$relavent_dist<-center$Niche_Sp_distance*2/(center$r_a+center$r_b)
    center$size_continent<-size_continent
    center$size<-size
    center$land_percentile<-center$size_continent/center$size
    center$width<-sqrt(size)
    
    #center$SD<-(sd(points_no_NA$V_1) + sd(points_no_NA$V_2))/2
    #rgl::plot3d(rgl::ellipse3d(mve$covariance,
    #                           centre = mve$centroid,
    #                           level = 0.95),
    #            alpha=0.4,col="blue")
    
    distance<-points %>% mutate(dist = sqrt((V_1-mve$centroid[1])^2 + 
                                              (V_2-mve$centroid[2])^2))
    r<-raster(matrix(distance$dist, nrow=sqrt(size)), crs=crs(r_env[[as.character(1)]]))
    NAvalue(r)<--9999
    extent(r)<-c(min(distance$x)-500, max(distance$x)+500, min(distance$y)-500, max(distance$y)+500)
    writeRaster(r, sprintf("%s/distribution.tif", folder), overwrite=T, NAflag=-9999)
    #raster(sprintf("%s/distribution.tif", folder))
    #saveRDS(r, sprintf("%s/r.rda", folder))
    #rr<-readRDS(sprintf("%s/r.rda", folder))
    #plot(rr)
    #
    #autoplot(pca, data = data.frame(mve$centroid), loadings = TRUE,
    #         loadings.colour = 'blue',
    #         loadings.label = TRUE, loadings.label.size = 3)
    
    #plot(r)
    #points(center$x, center$y)
    
    sample_size<-1000
    if (nrow(points_no_NA)<sample_size){
      sample_size<-nrow(points_no_NA)
    }
    
    points_sample<-points_no_NA[sample(nrow(points_no_NA), sample_size),]
    saveRDS(points_sample, sprintf("%s/points_sample.rda", folder))
    dists <- as.matrix(dist(cbind(points_sample$x, points_sample$y)))
    
    dists.inv <- 1/dists
    diag(dists.inv) <- 0
    
    #dists.inv[1:sample_size, 1:sample_size]
    center$V1_sd<-sd(points$V_1, na.rm=T)
    center$V2_sd<-sd(points$V_2, na.rm=T)
    
    m<-Moran.I(points_sample$V_1, dists.inv, na.rm = TRUE)
    center$M_V1_ovserved<-m$observed
    center$M_V1_expected<-m$expected
    center$M_V1_M_sd<-m$sd
    center$M_V1_p.value<-m$p.value
    m<-Moran.I(points_sample$V_2, dists.inv, na.rm = TRUE)
    center$M_V2_ovserved<-m$observed
    center$M_V2_expected<-m$expected
    center$M_V2_M_sd<-m$sd
    center$M_V2_p.value<-m$p.value
    if (is.null(result)){
      result<-center
    }else{
      result<-bind_rows(result, center)
    }
  }
  saveRDS(result, f)
}
