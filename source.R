#setwd("~/Experiments/Abundance_and_Niche/Script/abundance_and_niche")
library(dplyr)
library(raster)



max_dispersal_ability<-floor(res/2) # max dispersal distance

# define three ways of dispersl
dispersal_ability_linear<-function(dist, radius){
  probability<-(1-dist/radius)/2
  probability[probability<0]<-0
  return(probability)
}
dispersal_ability_ln<-function(dist, radius){
  probability<-(-log((dist+0.01)/radius))
  probability[probability<0]<-0
  return(probability)
}
dispersal_ability_fix<-function(dist, radius){
  probability<-radius
  return(probability)
}

#dispersal_ability_ln(c(0:4), 4)
#dispersal_ability_linear(c(0:4), 4)
#dispersal_ability_linear(c(0:4), 8)

dispersal_ability_ln(c(0:5), 5)
res<-31  # spatial extent
matrix<-expand.grid(x=c(1:res), y=c(1:res))
#define the distance between each pair of cell
rep<-100
distances<-as.matrix(dist(matrix))

uniform_background<-matrix
uniform_background$individual<-rep(1, res^2)
uniform_background$reproduction_rate<-rep(1, res^2)

individual<-res^2
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
  (max(gradient_background$individual)-min(gradient_background$individual))+0.7  #????? why
gradient_background$reproduction_rate <- gradient_background$reproduction_rate^2

diagonal_gradient_background<-matrix
diagonal_gradient_background$distance<-distances[1,]
dispersal_probability<-
  dispersal_ability_linear(distances[1,], ceiling(max(distances[1,])/2))+
  dispersal_ability_linear(distances[nrow(distances),], ceiling(max(distances[nrow(distances),])/2))
dispersal_probability<-
  dispersal_probability/sum(dispersal_probability)
diagonal_gradient_background$dispersal_probability<-dispersal_probability
diagonal_gradient_background$individual<-
  diagonal_gradient_background$dispersal_probability*individual
diagonal_gradient_background$reproduction_rate<-
  diagonal_gradient_background$individual /
  (max(diagonal_gradient_background$individual)-min(diagonal_gradient_background$individual))+0.7  #????? why
diagonal_gradient_background$reproduction_rate <- diagonal_gradient_background$reproduction_rate^2
r<-raster(matrix(diagonal_gradient_background$reproduction_rate, nrow=res))
plot(r)

left_gradient_background<-matrix
left_gradient_background$distance<-distances[1,]
dispersal_probability<-
  dispersal_ability_linear(distances[1,], ceiling(max(distances[1,])/2))+
  dispersal_ability_linear(distances[res,], ceiling(max(distances[res,])/2))
dispersal_probability<-
  dispersal_probability/sum(dispersal_probability)
left_gradient_background$dispersal_probability<-dispersal_probability
left_gradient_background$individual<-
  left_gradient_background$dispersal_probability*individual
left_gradient_background$reproduction_rate<-
  left_gradient_background$individual /
  (max(left_gradient_background$individual)-min(left_gradient_background$individual))+0.7  #????? why
left_gradient_background$reproduction_rate <- left_gradient_background$reproduction_rate^2
r<-raster(matrix(left_gradient_background$reproduction_rate, nrow=res))
plot(r)


hist(gradient_background$reproduction_rate)
hist(gradient_background$reproduction_rate^2)


max(gradient_background$reproduction_rate)/min(gradient_background$reproduction_rate)

r<-raster(matrix(gradient_background$individual, nrow=res))
plot(r)
r<-raster(matrix(gradient_background$reproduction_rate, nrow=res))
plot(r)

#gradient_background_ln<-matrix
#gradient_background_ln$distance<-distances[i,]
#gradient_background_ln$dispersal_probability<-
#  dispersal_ability_ln(gradient_background_ln$distance, ceiling(max(gradient_background_ln$distance)))
#gradient_background_ln$dispersal_probability<-
#  gradient_background_ln$dispersal_probability/sum(gradient_background_ln$dispersal_probability)
#gradient_background_ln$individual<-
#  gradient_background_ln$dispersal_probability*individual
#r<-raster(matrix(gradient_background_ln$individual, nrow=res))
#plot(r)

#random_background<-matrix
#random_background$distance<-distances[i,]
#random_background$individual<-runif(res^2)
#random_background$individual<-random_background$individual * res^2 / sum(random_background$individual)
#random_background$reproduction_rate<-gradient_background$reproduction_rate
#sum(random_background$individual)
#r<-raster(matrix(random_background$individual, nrow=res))
#plot(r)
#r<-raster(matrix(random_background$reproduction_rate, nrow=res))
#plot(r)


background<-left_gradient_background
#background<-gradient_background
#background<-gradient_background_ln
#background<-random_background

dispersal_ability_fun<-dispersal_ability_ln
#dispersal_ability_fun<-dispersal_ability_linear
#max_dispersal_ability<-ceiling(max(background$distance));dispersal_ability_fun<-dispersal_ability_fix


#max_dispersal_ability=20
#max_dispersal_ability=30
#max_dispersal_ability=1
#max_dispersal_ability=3
max_dispersal_ability=5
#max_dispersal_ability=10
#max_dispersal_ability=13
#max_dispersal_ability=16
i=1
#background$distance<-distances[100,]
#r<-raster(matrix(distances[100,], nrow=res))
  plot(raster(matrix(distances[res * (res/2) ,], nrow=res)))

for (j in c(1:rep)){
  print(j)
  r<-raster(matrix(background$individual, nrow=res))
  plot(r)
  background$new_individual<-0
  
  background$distance<-distances[res * (res/2),]
  neighbor_index<-as.numeric(rownames(distances)[distances[res * (res/2),]<=max_dispersal_ability])
  neighbors<-background[neighbor_index,]
  neighbors$dispersal_probability<-dispersal_ability_fun(neighbors$distance, max_dispersal_ability)
  max_pro<-sum(neighbors$dispersal_probability)
  
  for (i in c(1:nrow(background))){
    individual<-background[i, "individual"]
    background$distance<-distances[i,]
    neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability])
    neighbors<-background[neighbor_index,]
    neighbors$dispersal_probability<-dispersal_ability_fun(neighbors$distance, max_dispersal_ability)
    neighbors$dispersal_probability<-neighbors$dispersal_probability/max(neighbors$dispersal_probability)
    background[neighbor_index, "new_individual"]<-
      background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual
  }
  background$individual<-background$new_individual
  background$individual<-background$individual * background$reproduction_rate
  x<-readline(prompt="Press Enter to continue, Press X to exit: ")
  #x=""
  if (toupper(x)=="X"){
    break()
  }
}


#background[neighbors, "individual"]<-1
r<-raster(matrix(background$individual, nrow=res))
plot(r)
plot( r==max(values(r))  )
plot(log10(r))
