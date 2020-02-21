setwd("~/Experiments/Abundance_and_Niche/Script/abundance_and_niche")
library(dplyr)
library(raster)
res<-9
matrix<-expand.grid(x=c(1:res), y=c(1:res))
max_dispersal_ability<-4
#dispersal_ability<-c(0.5, 0.25, 0.15, 0.07, 0.03)
dispersal_ability<-function(dist, radius){
  probability<-(1-dist/radius)/2
  probability[probability<0]<-0
  return(probability)
}

plot(dispersal_ability, type="l")
rep<-100
distances<-as.matrix(dist(matrix))

uniform_background<-matrix
uniform_background$individual<-rep(1, res^2)

individual<-81
i=1
gradient_background<-matrix
gradient_background$distance<-distances[i,]
gradient_background$dispersal_probability<-
  dispersal_ability(gradient_background$distance, 13)
gradient_background$dispersal_probability<-
  gradient_background$dispersal_probability/sum(gradient_background$dispersal_probability)
gradient_background$individual<-
  gradient_background$dispersal_probability*individual
r<-raster(matrix(gradient_background$individual, nrow=res))
plot(r)

#background<-uniform_background
background<-gradient_background
for (j in c(1:rep)){
  print(j)
  background$new_individual<-0
  for (i in c(1:nrow(background))){
    individual<-background[i, "individual"]
    background$distance<-distances[i,]
    neighbor_index<-as.numeric(rownames(distances)[distances[i,]<=max_dispersal_ability])
    neighbors<-background[neighbor_index,]
    neighbors$dispersal_probability<-dispersal_ability(neighbors$distance, 4)
    neighbors$dispersal_probability<-neighbors$dispersal_probability/sum(neighbors$dispersal_probability)
    background[neighbor_index, "new_individual"]<-
      background[neighbor_index, "new_individual"]+neighbors$dispersal_probability*individual
  }
  background$individual<-background$new_individual
  r<-raster(matrix(background$individual, nrow=res))
  plot(r)
  x<-readline(prompt="Press Enter to continue, Press X to exit.")
  if (toupper(x)=="X"){
    break()
  }
}


#background[neighbors, "individual"]<-1
r<-raster(matrix(background$individual, nrow=res))
plot(r)
