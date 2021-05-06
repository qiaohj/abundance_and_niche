EXP_res<-31
a<-1

runExp_with_mask <- function(dispersal_ability_fun=dispersal_ability_exp,
                   max_dispersal_ability, 
                   background_object,
                   flag_gif=FALSE,flag_gif2=FALSE,
                   label="temp",
                   label_exp ="temp",
                   runs=10){
  if (!("is_valid" %in% colnames(background_object$outbg))){
    background_object$outbg$is_valid<-1
  }
  dir.create(paste0("fig_",label_exp))
  background <- background_object$outbg
  this_extent = background_object$res
  label=paste0(label,
               "range_",this_extent,"_",
               background_object$type,"_",
               background_object$toRatio,"_",
               "disp_",max_dispersal_ability)
  print(label)
  
  # define distance matrix again
  distances<-as.matrix(dist(expand.grid(x=c(1:this_extent), y=c(1:this_extent))))
  
  # calculate the max sum of the center point
  # mid_position <- ceiling(nrow(background)/2)
  # individual<-background[mid_position, "individual"]
  # background$distance<-distances[mid_position,]
  # neighbor_index<-as.numeric(rownames(distances)[distances[mid_position,]<=max_dispersal_ability])
  # neighbors<-background[neighbor_index,]
  # neighbors$dispersal_probability<-dispersal_ability_fun(neighbors$distance, max_dispersal_ability)
  # fixed_weight <- sum(neighbors$dispersal_probability)
  
  # a new way to calculate weight
  distancesMAX <-as.matrix(dist(expand.grid(x=c(1:(max_dispersal_ability+1)), y=c(1:(max_dispersal_ability+1) ))))
  dim(distancesMAX)
  mid_pos_max <- ceiling(dim(distancesMAX)[1]/2)
  dis_max <- distancesMAX[mid_pos_max,]
  fixed_weight <- sum( dispersal_ability_fun(dis_max, max_dispersal_ability) )
  
  all_data <- vector()
  plot(raster(matrix(background$individual, nrow=this_extent)),main="start")
  
  j=1 # do  j simulations
  for (j in c(1:runs)){
    #print(j)
    #r<-raster(matrix(background$individual, nrow=this_extent))
    #plot(r,main=j)
    background$new_individual<-0
    
    background$individual<-background$individual * background$reproduction_rate
    
    plot(raster(matrix(background$individual, nrow=this_extent)),main=paste0("time:",j,"; do:reproduce"))
    
    #save individuals
    background$time = j
    background$stage= "reproduce"
    all_data <- rbind(all_data,background[c("x",
                                            "y",
                                            "individual",
                                            "reproduction_rate",
                                            "new_individual",
                                            "time",
                                            "stage")])
    
    # a loop of each cell
    i=1
    for (i in c(1:nrow(background))){
      #cat("cell",i,"\n")
      individual<-background[i, "individual"]   # loop for every cell
      background$distance<-distances[i,]        # pre-calulated distance matrix
      neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability]) # draw a circle, find neighbors index
      neighbors<-background[neighbor_index,] # find nei.
      neighbors$dispersal_probability <- dispersal_ability_fun(neighbors$distance, max_dispersal_ability) # cal prob based on distance
      neighbors$dispersal_probability<-neighbors$dispersal_probability/fixed_weight
      flag_qiao=TRUE
      if(flag_qiao){
        background[neighbor_index, "new_individual"] <- background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual  # prob * pop size = in
      } else {
        # consider individuals move out
        #within a max distance, (sum of prob - prob_stay_at_home)* old_individual
        temp_index <- which(neighbors$distance==0)
        move_out <- neighbors$individual[temp_index] * (1-neighbors$dispersal_probability[temp_index]) # this is not related neighbor; only depends on dispersal ability 
        # condiser individuals come in (including stay_at_home)
        come_in <- sum(neighbors$dispersal_probability[-temp_index]*neighbors$individual [-temp_index])
        background[i, "new_individual"] =  background[i, "individual"] + come_in - move_out
      }
    }
    background$individual<-background$new_individual
    background[which(is.na(background$is_valid)), "individual"]<-0
    plot(raster(matrix(background$individual, nrow=this_extent)),main=paste0("time:",j,"; do:disperse"))
    
    background_object$outbg
    background$time = j
    background$stage= "dispersal"
    all_data <- rbind(all_data,background[names(all_data)])
  }# end of runs
  
  if(flag_gif){
    #ggplot(subset(all_data,time==100), aes(x = x, y = y, fill = individual)) +  geom_raster() + 
    #  scale_fill_gradient(low = "white",high = "red")+coord_quickmap()
    library(gganimate)
    all_data_sub <- all_data
    #all_data_sub <- subset(all_data_sub,time==1 | time==10 |  time==20 | time==100)
    setDT(all_data_sub)
    
    temp_max <- all_data_sub[,.(max_value_index=which.max(individual)),
                             by= time]
    #all_data <- merge(all_data,temp_max,by="time")
    
    #all_data_sub$individual_log <- log10(all_data_sub$individual)
    p <- ggplot(  all_data_sub,
                  aes(x = x, y = y, fill = individual )) +
      geom_raster() +  scale_fill_gradientn(colours = c("white","red","black") )+ #colours = c("white","red","black")
      #colours = terrain.colors(10) 
      geom_point()+ # show the highest point/s...............
      coord_quickmap()+
      labs(title = 'Year: {current_frame}', x = 'aaaaa', y = 'aaaaaaa') + view_follow()+
      transition_manual(time)
    pa <- gganimate::animate(
      plot = p, 
      nframes = length(unique(all_data_sub$time)), 
      fps = 1#, # 10 
      #end_pause = 2
    )
    pa
    anim_save(paste0("fig_",label_exp,"/",label,".gif"), animation = pa)
  }
  
  if(flag_gif2){
    library(gifski)
    
    png(paste0( "fig_",label_exp,"/dispersal",label,"%03d.png"))
    par(ask = FALSE)
    t=1
    for(t in 1:runs){
      temp_data <- subset(all_data,time==t & stage== "dispersal")
      max_index <- which(temp_data$individual==max(temp_data$individual))
      max_data <- temp_data[max_index,]
      temp_data <- temp_data[order(-temp_data$y),]
      r<-raster(matrix(temp_data$individual, nrow=this_extent,byrow=TRUE))
      extent(r) <- c(0,this_extent,0,this_extent)
      plot(r,main=paste0("time:",t))
      points(max_data$x-0.5,max_data$y-0.5,col="red",pch=20)
    }
    dev.off()
    png_files <- sprintf(paste0( "fig_",label_exp,"/dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( "fig_",label_exp,"/",label,"_dispersal.gif"))  )
    unlink(png_files)
    
    png(paste0( "fig_",label_exp,"/reproduce",label,"%03d.png"))
    par(ask = FALSE)
    t=1
    for(t in 1:runs){
      temp_data <- subset(all_data,time==t & stage== "reproduce")
      max_index <- which(temp_data$individual==max(temp_data$individual))
      max_data <- temp_data[max_index,]
      temp_data <- temp_data[order(-temp_data$y),]
      r<-raster(matrix(temp_data$individual, nrow=this_extent,byrow=TRUE))
      extent(r) <- c(0,this_extent,0,this_extent)
      plot(r,main=paste0("time:",t))
      points(max_data$x-0.5,max_data$y-0.5,col="red",pch=20)
    }
    dev.off()
    png_files <- sprintf(paste0( "fig_",label_exp,"/reproduce",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( "fig_",label_exp,"/",label,"_reproduce.gif"))  )
    unlink(png_files)
    
    png(paste0( "fig_",label_exp,"/reproduce_dispersal",label,"%03d.png"))
    par(ask = FALSE)
    t=1
    for(t in 1:runs){
      temp_data <- subset(all_data,time==t & stage== "reproduce")
      temp_data2 <- subset(all_data,time==t & stage== "dispersal")
      temp_data$individual <- (temp_data$individual+temp_data2$individual)/2
      max_index <- which(temp_data$individual==max(temp_data$individual))
      max_data <- temp_data[max_index,]
      temp_data <- temp_data[order(-temp_data$y),]
      r<-raster(matrix(temp_data$individual, nrow=this_extent,byrow=TRUE))
      extent(r) <- c(0,this_extent,0,this_extent)
      plot(r,main=paste0("time:",t))
      points(max_data$x-0.5,max_data$y-0.5,col="red",pch=20)
    }
    dev.off()
    png_files <- sprintf(paste0( "fig_",label_exp,"/reproduce_dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( "fig_",label_exp,"/",label,"_reproduce_dispersal.gif"))  )
    unlink(png_files)
  }# end of plot
  
  flag_curve=TRUE
  if(flag_curve){
    #png(paste0( "fig/reproduce_dispe_CURVE_",label,"%03d.png"))
    #par(ask = FALSE)
    t=1
    all_diag <- vector()
    for(t in 1:runs){
      temp_data1 <- subset(all_data,time==t & stage== "reproduce")
      temp_data2 <- subset(all_data,time==t & stage== "dispersal")
      temp_data <- temp_data1
      temp_data$individual <- (temp_data$individual+temp_data2$individual)/2
      
      temp_data <- temp_data[order(-temp_data$y),]
      temp_data1 <- temp_data1[order(-temp_data1$y),]
      temp_data2 <- temp_data2[order(-temp_data2$y),]
      
      my_diag <- function(aaa,bbb){
        temp_matrix <- matrix(aaa, nrow=bbb,byrow=TRUE)
        diag_cell <- vector()
        for(iii in 1:  bbb ){
          temp_cell <- temp_matrix[bbb+1-iii,iii]
          diag_cell <- c(diag_cell,temp_cell)
        }
        return(diag_cell)
      }
      
      one_diag0 <- data.frame(pos=1:this_extent,
                              value= my_diag(aaa=temp_data$individual,bbb=this_extent),
                              time=t,
                              type="mean")
      one_diag1 <- data.frame(pos=1:this_extent,
                              value= my_diag(aaa=temp_data1$individual,bbb=this_extent),
                              time=t,
                              type="reproduce")
      one_diag2 <- data.frame(pos=1:this_extent,
                              value= my_diag(aaa=temp_data2$individual,bbb=this_extent),
                              time=t,
                              type="dispersal")
      all_diag <- rbind(all_diag,one_diag0,one_diag1,one_diag2)
    }
    
    p0 <- ggplot(all_diag,aes(pos,value,col=time,group=time))+geom_point()+geom_line()+facet_wrap(~type)
    #p0
    ggsave(plot=p0,paste0("fig_",label_exp,"/curve_",label,".jpg"),width = 20,height = 20,units = "cm",dpi=300)
    #dev.off()
    #png_files <- sprintf(paste0( "fig/reproduce_dispersal",label,"%03d.png"), 1:runs)
    #gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( "fig/",label,"_reproduce_dispersal.gif"))  )
    #unlink(png_files)
    
    
    all_diag$extent <- this_extent
    all_diag$bg_type <- background_object$type
    all_diag$bg_ratio <- background_object$toRatio
    all_diag$max_dispersal_ability <- max_dispersal_ability
    saveRDS(all_diag,paste0("fig_",label_exp,"/curve_",label,".rds"))
  }  
  #plot the mean of the last 10 runs?
  
  
  all_data$extent <- this_extent
  all_data$bg_type <- background_object$type
  all_data$bg_ratio <- background_object$toRatio
  all_data$max_dispersal_ability <- max_dispersal_ability
  
  saveRDS(all_data,paste0("fig_",label_exp,"/raster_",label,".rds"))
  
  return(all_data)
}



#rect
background_object= generateBG(res=EXP_res,
                              type="uniform",
                              toRatio=1)
background_object$outbg$is_valid<-NA
background_object$outbg[which(background_object$outbg$x>10), "is_valid"]<-1
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "individual"]<-0
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "reproduction_rate"]<-0

plot(raster(matrix(background_object$outbg$is_valid,nrow=EXP_res)))

results0 <- runExp_with_mask(max_dispersal_ability=a,
                   background_object= background_object,
                   label="FINAL_",runs=20,flag_gif2=TRUE,
                   label_exp ="simple")

#circle
background_object= generateBG(res=EXP_res,
                              type="uniform",
                              toRatio=1)
background_object$outbg$is_valid<-NA
distance<-sqrt((background_object$outbg$x-16)^2+(background_object$outbg$y-16)^2)
background_object$outbg[(distance<15), "is_valid"]<-1
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "individual"]<-0
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "reproduction_rate"]<-0
plot(raster(matrix(background_object$outbg$is_valid,nrow=EXP_res)))

results0 <- runExp_with_mask(max_dispersal_ability=a,
                             background_object= background_object,
                             label="FINAL_",runs=20,flag_gif2=TRUE,
                             label_exp ="simple")

#triangle
background_object= generateBG(res=EXP_res,
                              type="uniform",
                              toRatio=1)
x_coord <- c(0.5, 30.5, 30.5, 0.5)
y_coord <- c(15.5, 0.5, 30.5, 15.5)
xym <- cbind(x_coord, y_coord)
xym
library(sp)
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

mask<-raster(matrix(background_object$outbg$individual,nrow=EXP_res), xmn = 0, xmx = 31, ymn = 0, ymx = 31)
res(mask)
mask_cut<-mask(mask, sps)
plot(mask)
points(xym[,1], xym[,2], col="black")
plot(mask_cut, add=T, col="blue")

mask_cut_p<-data.frame(rasterToPoints(mask_cut))

background_object$outbg$is_valid<-extract(mask_cut, background_object$outbg[, c("x", "y")])
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "individual"]<-0
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "reproduction_rate"]<-0
plot(raster(matrix(background_object$outbg$is_valid,nrow=EXP_res)))

results0 <- runExp_with_mask(max_dispersal_ability=a,
                             background_object= background_object,
                             label="FINAL_",runs=20,flag_gif2=TRUE,
                             label_exp ="simple")


#any polygon
background_object= generateBG(res=EXP_res,
                              type="uniform",
                              toRatio=1)

#writeRaster(mask, "~/temp/mask.tif")

x_coord <- runif(5) * 30.5
y_coord <- runif(5) * 30.5
xym <- cbind(x_coord, y_coord)
xym
library(sp)
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

mask<-raster(matrix(background_object$outbg$individual,nrow=EXP_res), xmn = 0, xmx = 31, ymn = 0, ymx = 31)
res(mask)
mask_cut<-mask(mask, sps)
plot(mask)
points(xym[,1], xym[,2], col="black")
plot(mask_cut, add=T, col="blue")

mask_cut_p<-data.frame(rasterToPoints(mask_cut))

background_object$outbg$is_valid<-extract(mask_cut, background_object$outbg[, c("x", "y")])
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "individual"]<-0
background_object$outbg[which(is.na(background_object$outbg$is_valid)), "reproduction_rate"]<-0
plot(raster(matrix(background_object$outbg$is_valid,nrow=EXP_res)))

results0 <- runExp_with_mask(max_dispersal_ability=a,
                             background_object= background_object,
                             label="FINAL_",runs=20,flag_gif2=TRUE,
                             label_exp ="simple")

