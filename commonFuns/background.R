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
