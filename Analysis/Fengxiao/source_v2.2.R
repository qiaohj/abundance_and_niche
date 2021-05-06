
library(raster)
library(ggplot2)

# load functions: dispersal_ability_exp,generateBG 
source("code/2020_10_27_myfunctions.R")

#generateBG(res=11,type="gradientSimple",toRatio=1)
#generateBG(res=11,type="gradientSimple",toRatio=2)
#dispersal_ability_fun=dispersal_ability_exp;background_object <- generateBG(res=7,type="gradient",toRatio=1);max_dispersal_ability=100;flag_gif=TRUE;label="test_";runs=10
#max_dispersal_ability=100;background=generateBG(res=11,type="gradient",toRatio=0.1);label="";runs=20;flag_gif2=TRUE;
#max_dispersal_ability=100;background_object=generateBG(res=11,type="gradientSimple",toRatio=1);label="";runs=20;flag_gif2=TRUE;dispersal_ability_fun=dispersal_ability_exp
#max_dispersal_ability=1;background_object= generateBG(res=31,type="uniform",toRatio=1);label="FINAL_";runs=20;flag_gif2=TRUE;dispersal_ability_fun=dispersal_ability_exp

runExp <- function(dispersal_ability_fun=dispersal_ability_exp,
         max_dispersal_ability, 
         background_object,
         flag_gif=FALSE,flag_gif2=FALSE,
         label="temp",
         label_exp ="temp",
         runs=10){
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
    plot(raster(matrix(background$individual, nrow=this_extent)),main=paste0("time:",j,"; do:disperse"))
    
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

#
flag_test=FALSE
if(flag_test){
  # test a uniform dispersal
  
# experiment with uniform start
results0 <- runExp(max_dispersal_ability=1,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)

results0 <- runExp(max_dispersal_ability=2,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)

results0 <- runExp(max_dispersal_ability=5,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)

results0 <- runExp(max_dispersal_ability=11,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)

results0 <- runExp(max_dispersal_ability=20,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)

results0 <- runExp(max_dispersal_ability=100,background_object = generateBG(res=11,type="uniform",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)


# experiment with gradient start
results0 <- runExp(max_dispersal_ability=0.1,background_object= generateBG(res=11,type="gradient",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# lead to a stable point between N-center and Sp-center

results0 <- runExp(max_dispersal_ability=1,background_object= generateBG(res=11,type="gradientSimple",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# lead to a stable point between N-center and Sp-center

results0 <- runExp(max_dispersal_ability=2,background_object= generateBG(res=11,type="gradient",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# lead to a stable point in between. But it is two diff curve with the same peak location.

results0 <- runExp(max_dispersal_ability=5,background_object= generateBG(res=11,type="gradientSimple",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# lead to a stable point pretty close to Sp-center.

results0 <- runExp(max_dispersal_ability=11,background_object= generateBG(res=11,type="gradient",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# seems further away from Sp-center, than disp=5.

results0 <- runExp(max_dispersal_ability=20 ,background_object= generateBG(res=11,type="gradient",toRatio=1),label="test_",runs=20,flag_gif2=TRUE)# gets closer to N-center, than disp=11

results0 <- runExp(max_dispersal_ability=100,background_object= generateBG(res=11,type="gradientSimple",toRatio=1),label="test_",runs=20,flag_gif2=TRUE) # I can see two peaks, one peak shows up at the end of dispersal, another peak shows up at the end of reproduction. But as dispersal get larger, the effect of dispersal becomes smaller (actually, the relative difference[ability to receive new individuals] between cells becomes smaller), then, by taking the mean of the two curves, the abundance center becomes the niche center


results0 <- runExp(max_dispersal_ability=100,background_object= generateBG(res=11,type="gradientSimple",toRatio=0.1),label="test_",runs=20,flag_gif2=TRUE) #by make env gradient smaller (close to uniform), the abundance center goes toward the Sp-center 

results0 <- runExp(max_dispersal_ability=100,background_object= generateBG(res=11,type="gradientSimple",toRatio=10 ),label="test_",runs=20,flag_gif2=TRUE) # by imposing a strong env-gradient, the abundance is almost always at N-center, when taking the mean.

}


### formal experiment!!
EXP_res <- 31
#temp_radius <- floor(EXP_res/2)*sqrt(2) # distance from corner to center.
temp_radius <- floor(EXP_res/2)+1   #~~~~~~~~~~~~~~~~~~~  change to 16 
EXP_disp <- c(1,2,temp_radius*0.5,temp_radius,
              temp_radius*1.2,temp_radius*sqrt(2) ,
              temp_radius*1.6,
              temp_radius*2,temp_radius*5) #,temp_radius*10
print(EXP_disp)

flag_rerun = FALSE
if (flag_rerun){
  for( a in EXP_disp){
    # for uniform env-gradient
    results0 <- runExp(max_dispersal_ability=a,
                       background_object= generateBG(res=EXP_res,
                                                     type="uniform",
                                                     toRatio=1),
                       label="FINAL_",runs=20,flag_gif2=TRUE,
                       label_exp ="simple")
    
    for(one_ratio in c(1,0.1,5,8,10) ){
      print(paste("this disp: ",a ,";   this ratio:",one_ratio) )
      results0 <- runExp(max_dispersal_ability=a,
                         background_object= generateBG(res=EXP_res,
                                                       type="gradientSimple",
                                                       toRatio=one_ratio),
                         label="FINAL_",runs=20,flag_gif2=TRUE,
                         label_exp ="gradientSimple")
    }
  }
}


# additional simulation of two niche center, random
for(eee in c("random","two_centers_3","two_centers_4") ){
  a = 1
  for( a in EXP_disp){
    for(one_ratio in c(1,0.1,5,8,10) ){
      results0 <- runExp(max_dispersal_ability=a,
                         background_object= generateBG(res=EXP_res,
                                                       type=eee,
                                                       toRatio=one_ratio),
                         label="FINAL_",runs=20,flag_gif2=TRUE,
                         label_exp =eee)
    }
  }
}

eee = "random"
for( a in c(1,2)){
  for(one_ratio in c(1,0.1,5,8,10) ){
    results0 <- runExp(max_dispersal_ability=a,
                       background_object= generateBG(res=EXP_res,
                                                     type=eee,
                                                     toRatio=one_ratio),
                       label="FINAL_",runs=20,flag_gif2=TRUE,
                       label_exp =eee)
  }
}
# a=1
# for( a in EXP_disp){
#   results0 <- runExp(max_dispersal_ability=a,
#                      background_object= generateBG(res=EXP_res,
#                                                    type="random",
#                                                    toRatio=1),
#                      label="FINAL_",runs=20,flag_gif2=TRUE,
#                      label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=0.1),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=0.5),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=2),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=10),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
# results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=0),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=5),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
#   results0 <- runExp(max_dispersal_ability=a,background_object= generateBG(res=EXP_res,type="random",toRatio=8),label="FINAL_",runs=20,flag_gif2=TRUE,label_exp ="random")
# 
# }



###################  end
