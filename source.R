setwd("~/Experiments/Abundance_and_Niche/Script/abundance_and_niche")
library(dplyr)
library(raster)
library(ggplot2)
library(rasterVis)
res<-31
matrix<-expand.grid(x=c(1:res), y=c(1:res))
max_dispersal_ability<-floor(res/2)
#max_dispersal_ability<-floor(res/10)
#dispersal_ability<-c(0.5, 0.25, 0.15, 0.07, 0.03)
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

dispersal_ability_ln(c(0:4), 4)
#plot(dispersal_ability_linear, type="l")
rep<-100
distances<-as.matrix(dist(matrix))

uniform_background<-matrix
uniform_background$individual<-rep(1, res^2)
uniform_background$reproduction_rate<-1

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
  (max(gradient_background$individual)-min(gradient_background$individual)) * 2
r<-raster(matrix(gradient_background$individual, nrow=res))
plot(r)
r<-raster(matrix(gradient_background$reproduction_rate, nrow=res))
plot(r)
gradient_background$individual<-rep(1, res^2)


gradient_background_ln<-matrix
gradient_background_ln$distance<-distances[i,]
gradient_background_ln$dispersal_probability<-
  dispersal_ability_ln(gradient_background_ln$distance, ceiling(max(gradient_background_ln$distance)))
gradient_background_ln$dispersal_probability<-
  gradient_background_ln$dispersal_probability/sum(gradient_background_ln$dispersal_probability)
gradient_background_ln$individual<-
  gradient_background_ln$dispersal_probability*individual
r<-raster(matrix(gradient_background_ln$individual, nrow=res))
plot(r)

random_background<-matrix
random_background$distance<-distances[i,]
random_background$individual<-runif(res^2)
random_background$individual<-random_background$individual * res^2 / sum(random_background$individual)
random_background$reproduction_rate<-gradient_background$reproduction_rate
sum(random_background$individual)
r<-raster(matrix(random_background$individual, nrow=res))
plot(r)
r<-raster(matrix(random_background$reproduction_rate, nrow=res))
plot(r)


#background<-uniform_background
background<-gradient_background
#background<-gradient_background_ln
#background<-random_background

#dispersal_ability_fun<-dispersal_ability_ln
dispersal_ability_fun<-dispersal_ability_linear
#dispersal_ability_fun<-dispersal_ability_fix
#max_dispersal_ability<-ceiling(max(background$distance));dispersal_ability_fun<-dispersal_ability_fix


for (j in c(1:rep)){
  print(j)
  r<-raster(matrix(background$individual, nrow=res))
  persp(r, exp=1,phi=35, xlab="Longitude", ylab="Latitude", zlab="Elevation")
  #plot(r)
  background$new_individual<-0
  for (i in c(1:nrow(background))){
    individual<-background[i, "individual"]
    background$distance<-distances[i,]
    neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability])
    neighbors<-background[neighbor_index,]
    neighbors$dispersal_probability<-dispersal_ability_fun(neighbors$distance, max_dispersal_ability)
    neighbors$dispersal_probability<-neighbors$dispersal_probability/sum(neighbors$dispersal_probability)
    background[neighbor_index, "new_individual"]<-
      background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual
  }
  background$individual<-background$new_individual
  background$individual<-background$individual * background$reproduction_rate
  
  #x<-readline(prompt="Press Enter to continue, Press X to exit: ")
  x=""
  if (toupper(x)=="X"){
    break()
  }
}


#background[neighbors, "individual"]<-1
r<-raster(matrix(background$individual, nrow=res))
plot(r)


a<-rnorm(1000, -5, 1)
b<-rnorm(1000, 5, 1)
hist(a, xlim=c(-10, 10), border="blue")
hist(b, add=T, col=NA, border="red")

hist(a*b, col=NA, border="black")
hist(a+b, col=NA, border="black")
