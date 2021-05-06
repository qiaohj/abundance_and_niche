library(raster)

dispersal_ability_exp<-function(dist, radius,ratio=1){ # check Jorge's value!!
  probability <- exp(-1*(ratio)* (dist)/radius)   # folliwng Jorge's equation
  probability[probability<0]<-0
  
  # added by xiao
  probability [ dist>radius] <- 0 # make sure this is the same logic used in the loop
  return(probability)
}

# not used
# dispersal_ability_uniform<-function(dist, radius,ratio=1){
#   probability <- rep(1,length(dist))
#   probability[probability<0]<-0
#   
#   # added by xiao
#   probability [ dist>radius] <- 0 # make sure this is the same logic used in the loop
#   return(probability)
# }



# x=1:10;max_x = 10
# plot(x, decay_linear(x,max_x,slope=10))
# points(x, decay_linear(x,max_x,slope=2),col="red")
# points(x, decay_linear(x,max_x,slope=1),col="blue")
# points(x, decay_linear(x,max_x,slope=0.5),col="green")

# # used in two center experiment (qiao)
dispersal_ability_linear<-function(dist, radius){
  probability<-(1-dist/radius)/2
  probability[probability<0]<-0
  return(probability)
}
# (xiao)
decay_linear<-function(dist, max_dist,slope=1){
  y = 1*slope - slope* dist /max_dist
  return(y)
}

decay_linear_2center<-function(dist, max_dist,slope=1){
  y = 1*slope - slope* dist /max_dist
  y[y<0] = 0
  return(y)
}

# x = 1:100
# y1 = dispersal_ability_linear(x,100)
# y2 = decay_linear(x,100,1)
# plot(x,y1)
# points(x,y2,col="red")

generateBG <- function(res=11,type="gradient",toRatio=1){
  
  #res<-31  # spatial extent
  matrix<-expand.grid(x=c(1:res), y=c(1:res))
  #define the distance between each pair of cell
  distances<-as.matrix(dist(matrix))
  raw_pop <- res^2
  if(type=="uniform"){
    uniform_background<-matrix
    uniform_background$individual<-rep(raw_pop, res^2)
    uniform_background$reproduction_rate<-rep(1, res^2) # 1 means no growth
    #uniform_background$reproduction_rate<-rep(1.1, res^2)
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
      (max(gradient_background$individual)-min(gradient_background$individual))+1#0.7  #????? why
    
    # make the ratio larger or smaller
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
    #plot(raster(matrix(outbg$individual,nrow=res)))
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
    #plot(raster(matrix(outbg$individual,nrow=res)))
    plot(raster(matrix(outbg$reproduction_rate,nrow=res)))
    # #对角，通过调整标红的两个值（圆心坐标），来实现不同的gradient，通过标蓝的数值，来调整圆的半径
    # diagonal_gradient_background<-matrix
    # diagonal_gradient_background$distance<-distances[1,]
    # dispersal_probability<-
    #   dispersal_ability_linear(distances[1,], ceiling(max(distances[1,])/2))+
    #   dispersal_ability_linear(distances[nrow(distances),], ceiling(max(distances[nrow(distances),])/2))
    # dispersal_probability<-
    #   dispersal_probability/sum(dispersal_probability)
    # diagonal_gradient_background$dispersal_probability<-dispersal_probability
    # diagonal_gradient_background$individual<-
    #   diagonal_gradient_background$dispersal_probability*raw_pop
    # diagonal_gradient_background$reproduction_rate<-
    #   diagonal_gradient_background$individual /
    #   (max(diagonal_gradient_background$individual)-min(diagonal_gradient_background$individual))+0.7  #????? why
    # 
    # diagonal_gradient_background$reproduction_rate <- diagonal_gradient_background$reproduction_rate^(toRatio)
    # 
    # outbg <- diagonal_gradient_background

  }
  if(type=="two_centers_2"){
    #这里我弄了一个左上，左下两个中心，你可以修改res为任何一个值，来改变第二个gradient的中心。
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
    #plot(raster(matrix(outbg$individual,nrow=res)))
    #plot(raster(matrix(outbg$reproduction_rate,nrow=res)))
  }  
  if(type=="random" ){
    random_background <- matrix
    random_background$distance <- distances[1,]
    set.seed(1)
    #random_background$individual <- runif(res^toRatio)
    #random_background$individual <- random_background$individual * res^2 / sum(random_background$individual)
    #random_background$reproduction_rate<-random_background$individual
    random_background$reproduction_rate = runif(n=nrow(random_background),min=0,max=toRatio)+1
    
    #sum(random_background$individual)
    #r<-raster(matrix(random_background$individual, nrow=res))
    #plot(r)
    #plot(r==max(values(r)))
    #r<-raster(matrix(random_background$reproduction_rate, nrow=res))
    #plot(r)
    outbg <- random_background
  }
  
  
  
  #************** make sure starting pop size are the same.**************
  outbg$individual <- raw_pop
  plot(stack(   raster(matrix(outbg$individual,nrow=res)),
                raster(matrix(outbg$reproduction_rate,nrow=res))   ))
  print( max(outbg$reproduction_rate)/min(outbg$reproduction_rate) )
  
  results <- list(outbg=outbg,
                  res=res,
                  type=type,
                  toRatio=toRatio)
  return (results)
}

#generateBG(res=11,type="random",toRatio=10)

#generateBG(res=11,type="two_centers_1",toRatio=10)
#generateBG(res=11,type="two_centers_2",toRatio=10) # to do, make it consistent with simple Gradient.

generateBG(res=11,type="two_centers_3",toRatio=0.1)
#generateBG(res=11,type="two_centers_4",toRatio=20)


#generateBG(res=11,type="two_centers_2",toRatio=10)

#generateBG(res=11,type="gradientSimple",toRatio=1)
