library(raster)
library(dplyr)
setwd("/media/huijieqiao/WD12T/Experiments/abundance_and_niche/abundance_and_niche")
base<-"/media/huijieqiao/WD12T/Experiments/IUCN_FIX/Raster/IUCN_Range_By_Species_eck4"
group<-"Amphibians"
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
files<-list.files(sprintf("%s/%s", base, group), pattern = "*.tif")
f<-files[1]
result<-NULL
for (f in files){
  print(f)
  r<-raster(sprintf("%s/%s/%s", base, group, f))
  n<-length(values(r)[!is.na(values(r))])
  item<-as_tibble(data.frame(group=group, species=gsub("\\.tif", "", f), N=n, stringsAsFactors = F))
  if (is.null(result)){
    result<-item
  }else{
    result<-bind_rows(result, item)
  }
}

saveRDS(result, sprintf("../Tables/Distribution_Area/%s.rda", group))
