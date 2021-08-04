library(raster)

dispersal_ability_exp<-function(dist, radius,ratio=1){ 
  probability <- exp(-1*(ratio)* (dist)/radius)   
  probability[probability<0]<-0  
  probability [ dist>radius] <- 0 
  return(probability)
}

dispersal_ability_linear<-function(dist, radius){
  probability<-(1-dist/radius)/2
  probability[probability<0]<-0
  return(probability)
}
decay_linear<-function(dist, max_dist,slope=1){
  y = 1*slope - slope* dist /max_dist
  return(y)
}

decay_linear_2center<-function(dist, max_dist,slope=1){
  y = 1*slope - slope* dist /max_dist
  y[y<0] = 0
  return(y)
}

generateBG <- function(res=11,type="gradient",toRatio=1,mask_Type=NULL){
  
  matrix<-expand.grid(x=c(1:res), y=c(1:res))
  distances<-as.matrix(dist(matrix))
  raw_pop <- res^2
  if(type=="uniform"){
    uniform_background<-matrix
    uniform_background$individual<-rep(raw_pop, res^2)
    uniform_background$reproduction_rate<-rep(1, res^2) 
    outbg <- uniform_background
  }
  if(type=="gradient"){
    individual<-raw_pop
    i=1
    gradient_background<-matrix
    gradient_background$distance<-distances[i,]
    gradient_background$dispersal_probability<-
      dispersal_ability_linear(gradient_background$distance, ceiling(max(gradient_background$distance)))
    gradient_background$dispersal_probability<-
      gradient_background$dispersal_probability/sum(gradient_background$dispersal_probability)
    gradient_background$individual<-
      gradient_background$dispersal_probability*individual
    gradient_background$reproduction_rate<-
      gradient_background$individual / 
      (max(gradient_background$individual)-min(gradient_background$individual))+1    
    gradient_background$reproduction_rate <- gradient_background$reproduction_rate^(toRatio)
    outbg <- gradient_background    
  }
  if(type=="gradientSimple"){
    individual<-raw_pop
    i=1
    gradient_background<-matrix
    gradient_background$distance<-distances[i,]
    gradient_background$reproduction_rate <- 
      decay_linear(dist = gradient_background$distance,
                   max_dist = max(gradient_background$distance),
                   slope = toRatio) +1
    outbg <- gradient_background
  }
  if(type=="two_centers_4"){
    individual<-raw_pop
    i=1
    gradient_background<-matrix
    gradient_background$distance<-distances[i,]
    gradient_background$reproduction_rate <-
      (  decay_linear_2center(dist = gradient_background$distance,
                              max_dist = max(gradient_background$distance)/2 ,
                              slope = toRatio) +
           
           decay_linear_2center(dist = distances[res,]  ,
                                max_dist = max(  distances[nrow(distances),]    )/2 ,
                                slope = toRatio)  )     
    gradient_background$reproduction_rate  = gradient_background$reproduction_rate +1  
    outbg <- gradient_background
    plot(raster(matrix(outbg$reproduction_rate,nrow=res)))
  }
  if(type=="two_centers_3"){
    individual<-raw_pop
    i=1
    gradient_background<-matrix
    gradient_background$distance<-distances[i,]
    gradient_background$reproduction_rate <-
    (  decay_linear_2center(dist = gradient_background$distance,
                 max_dist = max(gradient_background$distance)/2 ,
                 slope = toRatio) +
         
         decay_linear_2center(dist = distances[nrow(distances),]  ,
                 max_dist = max(  distances[nrow(distances),]    )/2 ,
                 slope = toRatio)  ) 
    
    gradient_background$reproduction_rate  = gradient_background$reproduction_rate +1  
    outbg <- gradient_background
    plot(raster(matrix(outbg$reproduction_rate,nrow=res)))
  }
  if(type=="two_centers_2"){
    left_gradient_background<-matrix
    left_gradient_background$distance<-distances[1,]
    dispersal_probability<-
      dispersal_ability_linear(distances[1,], ceiling(max(distances[1,])/2))+
      dispersal_ability_linear(distances[res,], ceiling(max(distances[res,])/2))
    dispersal_probability<-
      dispersal_probability/sum(dispersal_probability)
    left_gradient_background$dispersal_probability<-dispersal_probability
    left_gradient_background$individual<-
      left_gradient_background$dispersal_probability*raw_pop
    left_gradient_background$reproduction_rate<-
      left_gradient_background$individual /
      (max(left_gradient_background$individual)-min(left_gradient_background$individual))+0.7  #????? why
    left_gradient_background$reproduction_rate <- left_gradient_background$reproduction_rate^(toRatio)
    outbg <- left_gradient_background
  }  
  if(type=="random" ){
    random_background <- matrix
    random_background$distance <- distances[1,]
    set.seed(1)
    random_background$reproduction_rate = runif(n=nrow(random_background),min=0,max=toRatio)+1
    outbg <- random_background
  }  
  outbg$individual <- raw_pop   
  if (!is.null(mask_Type) ){
  if(mask_Type=="rectangle"){
    mask_Para = 2* (res-1)/3 +1
    outbg$is_valid<-NA
    outbg[which(outbg$x < mask_Para), "is_valid"]<-1
    outbg[which(is.na(outbg$is_valid)), "individual"]<-0
    outbg[which(is.na(outbg$is_valid)), "reproduction_rate"]<-0
  }
  if(mask_Type=="circle"){
    outbg$is_valid<-NA
    distance<-sqrt((outbg$x - ceiling(res/2)  )^2 + (outbg$y - ceiling(res/2)   )^2)
    outbg[(distance  <  floor(res/2)  ), "is_valid"]<-1
    outbg[which(is.na(outbg$is_valid)), "individual"]<-0
    outbg[which(is.na(outbg$is_valid)), "reproduction_rate"]<-0
  }
  if(mask_Type=="triangle"){	
	x_coord <- c(ceiling(res/2)-0.5,
						0.5, 
						res-0.5, 
						ceiling(res/2)-0.5     )
    y_coord <- c( res-0.5,
					0.5, 
					0.5,
					res-0.5)	
    xym <- cbind(x_coord, y_coord)
    library(sp)
    p = Polygon(xym)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    mask<-raster(matrix(outbg$individual,nrow=res), xmn = 0, xmx = res, ymn = 0, ymx = res)
    mask_cut<-mask(mask, sps)
    mask_cut_p<-data.frame(rasterToPoints(mask_cut))
    outbg$is_valid<-extract(mask_cut, outbg[, c("x", "y")])
    outbg[which(is.na(outbg$is_valid)), "individual"]<-0
    outbg[which(is.na(outbg$is_valid)), "reproduction_rate"]<-0    
  }
  if(mask_Type=="random"){
    x_coord <- runif(5) * (res-0.5)
    y_coord <- runif(5) * (res-0.5)
    xym <- cbind(x_coord, y_coord)
    library(sp)
    p = Polygon(xym)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    mask<-raster(matrix(outbg$individual,nrow=res), xmn = 0, xmx = res, ymn = 0, ymx = res)
    mask_cut<-mask(mask, sps)
    mask_cut_p<-data.frame(rasterToPoints(mask_cut))
    outbg$is_valid<-extract(mask_cut, outbg[, c("x", "y")])
    outbg[which(is.na(outbg$is_valid)), "individual"]<-0
    outbg[which(is.na(outbg$is_valid)), "reproduction_rate"]<-0
  }
  }  
  plot(stack(   raster(matrix(outbg$individual,nrow=res)),
                raster(matrix(outbg$reproduction_rate,nrow=res))   ))  
  results <- list(outbg=outbg,
                  res=res,
                  type=type,
                  toRatio=toRatio)
  return (results)
}

runExp <- function(dispersal_ability_fun=dispersal_ability_exp,
                   max_dispersal_ability, 
                   background_object,
                   flag_gif=FALSE,flag_gif2=FALSE,
                   label="temp",
                   label_exp ="temp",
                   runs=10,
                   dispersalMode = "passitve"){
  outFolder = paste0("experiments/","fig_",label_exp,"_",dispersalMode)
  dir.create(  outFolder  )
  background <- background_object$outbg
  this_extent = background_object$res
  label=paste0(label,
               "range_",this_extent,"_",
               background_object$type,"_",
               background_object$toRatio,"_",
               "disp_",max_dispersal_ability)
  print(label)
  
  distances<-as.matrix(dist(expand.grid(x=c(1:this_extent), y=c(1:this_extent))))  
  distancesMAX <-as.matrix(dist(expand.grid(x=c(1:(max_dispersal_ability+1)), y=c(1:(max_dispersal_ability+1) ))))
  dim(distancesMAX)
  mid_pos_max <- ceiling(dim(distancesMAX)[1]/2)
  dis_max <- distancesMAX[mid_pos_max,]
  fixed_weight <- sum( dispersal_ability_fun(dis_max, max_dispersal_ability) )  
  all_data <- vector()
  plot(raster(matrix(background$individual, nrow=this_extent)),main="start")  
  j=1 
  for (j in c(1:runs)){
    background$new_individual<-0    
    background$individual<-background$individual * background$reproduction_rate    
    plot(raster(matrix(background$individual, nrow=this_extent)),main=paste0("time:",j,"; do:reproduce"))    
    background$time = j
    background$stage= "reproduce"
    all_data <- rbind(all_data,background[c("x",
                                            "y",
                                            "individual",
                                            "reproduction_rate",
                                            "new_individual",
                                            "time",
                                            "stage")])    
    i=1
    for (i in c(1:nrow(background))){
      individual<-background[i, "individual"]  
      background$distance<-distances[i,]        
      neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability]) 
      neighbors<-background[neighbor_index,] 
      neighbors$dispersal_probability <- dispersal_ability_fun(neighbors$distance, max_dispersal_ability)      
    if(dispersalMode == "active"){
	  suitability_weight = neighbors$reproduction_rate / sum(neighbors$reproduction_rate)	  
	  neighbors$dispersal_probability  = suitability_weight  * (neighbors$dispersal_probability/fixed_weight)
      neighbors$dispersal_probability  = neighbors$dispersal_probability /sum(neighbors$dispersal_probability)	
    }
      if(dispersalMode == "passitve"){
      neighbors$dispersal_probability<-neighbors$dispersal_probability/fixed_weight
      }	  
        background[neighbor_index, "new_individual"] <- background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual 
    }    
    background$individual<-background$new_individual    
    plot(raster(matrix(background$individual, nrow=this_extent)),main=paste0("time:",j,"; do:disperse"))    
    background$time = j
    background$stage= "dispersal"
    all_data <- rbind(all_data,background[names(all_data)])
  }
  
  if(flag_gif){
    library(gganimate)
    all_data_sub <- all_data
    setDT(all_data_sub)    
    temp_max <- all_data_sub[,.(max_value_index=which.max(individual)),
                             by= time]    
    p <- ggplot(  all_data_sub,
                  aes(x = x, y = y, fill = individual )) +
      geom_raster() +  scale_fill_gradientn(colours = c("white","red","black") )+ 
      geom_point()+ 
      coord_quickmap()+
      labs(title = 'Year: {current_frame}', x = 'aaaaa', y = 'aaaaaaa') + view_follow()+
      transition_manual(time)
    pa <- gganimate::animate(
      plot = p, 
      nframes = length(unique(all_data_sub$time)), 
      fps = 1
    )
    pa
    anim_save(paste0(outFolder,"/",label,".gif"), animation = pa)
  }
  
  if(flag_gif2){
    library(gifski)    
    png(paste0( outFolder,"/dispersal",label,"%03d.png"))
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
    png_files <- sprintf(paste0(outFolder,"/dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0(outFolder,"/",label,"_dispersal.gif"))  )
    unlink(png_files)    
    png(paste0( outFolder,"/reproduce",label,"%03d.png"))
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
    png_files <- sprintf(paste0( outFolder,"/reproduce",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( outFolder,"/",label,"_reproduce.gif"))  )
    unlink(png_files)    
    png(paste0(outFolder,"/reproduce_dispersal",label,"%03d.png"))
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
    png_files <- sprintf(paste0( outFolder,"/reproduce_dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0(outFolder,"/",label,"_reproduce_dispersal.gif"))  )
    unlink(png_files)
  }
  
  flag_curve=TRUE
  if(flag_curve){
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
    ggsave(plot=p0,paste0(outFolder,"/curve_",label,".jpg"),width = 20,height = 20,units = "cm",dpi=300)    
    all_diag$extent <- this_extent
    all_diag$bg_type <- background_object$type
    all_diag$bg_ratio <- background_object$toRatio
    all_diag$max_dispersal_ability <- max_dispersal_ability
    saveRDS(all_diag,paste0(outFolder,"/curve_",label,".rds"))
  } 
  all_data$extent <- this_extent
  all_data$bg_type <- background_object$type
  all_data$bg_ratio <- background_object$toRatio
  all_data$max_dispersal_ability <- max_dispersal_ability  
  saveRDS(all_data,paste0(outFolder,"/raster_",label,".rds"))  
  return(all_data)
}

runExp_mask <- function(dispersal_ability_fun=dispersal_ability_exp,
                   max_dispersal_ability, 
                   background_object,
                   flag_gif=FALSE,flag_gif2=FALSE,
                   label="temp",
                   label_exp ="temp",
                   runs=10){
  outFolder = paste0("experiments/","fig_",label_exp)
  dir.create(  outFolder  )  
  if (!("is_valid" %in% colnames(background_object$outbg))){
    background_object$outbg$is_valid<-1
  }
  background <- background_object$outbg
  this_extent = background_object$res
  label=paste0(label,
               "range_",this_extent,"_",
               background_object$type,"_",
               background_object$toRatio,"_",
               "disp_",max_dispersal_ability)
  print(label)  
  distances<-as.matrix(dist(expand.grid(x=c(1:this_extent), y=c(1:this_extent))))  
  distancesMAX <-as.matrix(dist(expand.grid(x=c(1:(max_dispersal_ability+1)), 
                                            y=c(1:(max_dispersal_ability+1) ))))
  dim(distancesMAX)
  mid_pos_max <- ceiling(dim(distancesMAX)[1]/2)
  dis_max <- distancesMAX[mid_pos_max,]
  fixed_weight <- sum( dispersal_ability_fun(dis_max, max_dispersal_ability) )
  
  all_data <- vector()
  plot(raster(matrix(background$individual, nrow=this_extent)),main="start")  
  j=1
  for (j in c(1:runs)){
    background$new_individual<-0    
    background$individual<-background$individual * background$reproduction_rate    
	background[which(is.na(background$is_valid)), "individual"]<-0     
    background$time = j
    background$stage= "reproduce"
    all_data <- rbind(all_data,background[c("x",
                                            "y",
                                            "individual",
                                            "reproduction_rate",
                                            "new_individual",
                                            "time",
                                            "stage")])
    i=1
    for (i in c(1:nrow(background))){
      if(   is.na( background$is_valid[i])     )   {
        next
      }
      individual<-background[i, "individual"]   
      background$distance<-distances[i,]        
      background$distance[is.na(background$is_valid)] = 1E100
      
      
      neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability])
      neighbors<-background[neighbor_index,] 
      neighbors$dispersal_probability <- dispersal_ability_fun(neighbors$distance, max_dispersal_ability) 
      neighbors$dispersal_probability<-neighbors$dispersal_probability/fixed_weight
      flag_qiao=TRUE
      if(flag_qiao){
        background[neighbor_index, "new_individual"] <- background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual  
      } else {
        temp_index <- which(neighbors$distance==0)
        move_out <- neighbors$individual[temp_index] * (1-neighbors$dispersal_probability[temp_index]) 
        come_in <- sum(neighbors$dispersal_probability[-temp_index]*neighbors$individual [-temp_index])
        background[i, "new_individual"] =  background[i, "individual"] + come_in - move_out
      }
    }
    background$individual<-background$new_individual	
    background[which(is.na(background$is_valid)), "individual"]<-0 
    plot(raster(matrix(log(background$individual), nrow=this_extent)),
         main=paste0("time:",j,"; do:disperse"))
    background$time = j
    background$stage= "dispersal"
    all_data <- rbind(all_data,background[names(all_data)])
  }
  
  if(flag_gif){
    library(gganimate)
    all_data_sub <- all_data
    setDT(all_data_sub)    
    temp_max <- all_data_sub[,.(max_value_index=which.max(individual)),
                             by= time]    
    p <- ggplot(  all_data_sub,
                  aes(x = x, y = y, fill = individual )) +
      geom_raster() +  scale_fill_gradientn(colours = c("white","red","black") )+ 
      geom_point()+ 
      coord_quickmap()+
      labs(title = 'Year: {current_frame}', x = 'aaaaa', y = 'aaaaaaa') + view_follow()+
      transition_manual(time)
    pa <- gganimate::animate(
      plot = p, 
      nframes = length(unique(all_data_sub$time)), 
      fps = 1
    )
    pa
    anim_save(paste0(outFolder,"/",label,".gif"), animation = pa)
  }  
  if(flag_gif2){
    library(gifski)    
    png(paste0( outFolder,"/dispersal",label,"%03d.png"))
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
    png_files <- sprintf(paste0(outFolder,"/dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0(outFolder,"/",label,"_dispersal.gif"))  )
    unlink(png_files)    
    png(paste0( outFolder,"/reproduce",label,"%03d.png"))
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
    png_files <- sprintf(paste0( outFolder,"/reproduce",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0( outFolder,"/",label,"_reproduce.gif"))  )
    unlink(png_files)    
    png(paste0(outFolder,"/reproduce_dispersal",label,"%03d.png"))
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
    png_files <- sprintf(paste0(outFolder,"/reproduce_dispersal",label,"%03d.png"), 1:runs)
    gif_file <- gifski(png_files=png_files,gif_file=paste0(paste0(outFolder,"/",label,"_reproduce_dispersal.gif"))  )
    unlink(png_files)
  }
  
  flag_curve=TRUE
  if(flag_curve){
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
    ggsave(plot=p0,paste0(outFolder,"/curve_",label,".jpg"),width = 20,height = 20,units = "cm",dpi=300)    
    all_diag$extent <- this_extent
    all_diag$bg_type <- background_object$type
    all_diag$bg_ratio <- background_object$toRatio
    all_diag$max_dispersal_ability <- max_dispersal_ability
    saveRDS(all_diag,paste0(outFolder,"/curve_",label,".rds"))
  }   
  all_data$extent <- this_extent
  all_data$bg_type <- background_object$type
  all_data$bg_ratio <- background_object$toRatio
  all_data$max_dispersal_ability <- max_dispersal_ability  
  saveRDS(all_data,paste0(outFolder,"/raster_",label,".rds"))  
  return(all_data)
}

  