library(raster)
library(ggplot2)
source("code/0_2021_05_30_myfunctions.R")

EXP_res <- 31
temp_radius <- floor(EXP_res/2)+1   
EXP_disp <- c(1,2,temp_radius*0.5,temp_radius,
              temp_radius*1.2,temp_radius*sqrt(2) ,
              temp_radius*1.6,
              temp_radius*2,temp_radius*5) 
print(EXP_disp)

flag_rerun = FALSE
if (flag_rerun){
  for( a in EXP_disp){
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

for( a in EXP_disp){
  results0 <- runExp_mask(max_dispersal_ability=a,
                     background_object= generateBG(res=EXP_res,
                                                   type="uniform",
                                                   toRatio=1,
                                                   mask_Type="rectangle"),
                     label="FINAL_",runs=20,flag_gif2=TRUE,
                     label_exp ="gradientSimple_rectangle")
  
  for(one_ratio in c(1,0.1,5,8,10) ){
    print(paste("this disp: ",a ,";   this ratio:",one_ratio) )
    results0 <- runExp_mask(max_dispersal_ability=a,
                       background_object= generateBG(res=EXP_res,
                                                     type="gradientSimple",
                                                     toRatio=one_ratio,
                                                     mask_Type="rectangle"),
                       label="FINAL_",runs=20,flag_gif2=TRUE,
                       label_exp ="gradientSimple_rectangle")
  }
}

for( a in EXP_disp){
  # for uniform env-gradient
  results0 <- runExp_mask(max_dispersal_ability=a,
                          background_object= generateBG(res=EXP_res,
                                                        type="uniform",
                                                        toRatio=1,
                                                        mask_Type="triangle"),
                          label="FINAL_",runs=20,flag_gif2=TRUE,
                          label_exp ="gradientSimple_triangle")
  
  for(one_ratio in c(1,0.1,5,8,10) ){
    print(paste("this disp: ",a ,";   this ratio:",one_ratio) )
    results0 <- runExp_mask(max_dispersal_ability=a,
                            background_object= generateBG(res=EXP_res,
                                                          type="gradientSimple",
                                                          toRatio=one_ratio,
                                                          mask_Type="triangle"),
                            label="FINAL_",runs=20,flag_gif2=TRUE,
                            label_exp ="gradientSimple_triangle")
  }
}

for( a in EXP_disp){
  results0 <- runExp_mask(max_dispersal_ability=a,
                          background_object= generateBG(res=EXP_res,
                                                        type="uniform",
                                                        toRatio=1,
                                                        mask_Type="circle"),
                          label="FINAL_",runs=20,flag_gif2=TRUE,
                          label_exp ="gradientSimple_circle")
  
  for(one_ratio in c(1,0.1,5,8,10) ){
    print(paste("this disp: ",a ,";   this ratio:",one_ratio) )
    results0 <- runExp_mask(max_dispersal_ability=a,
                            background_object= generateBG(res=EXP_res,
                                                          type="gradientSimple",
                                                          toRatio=one_ratio,
                                                          mask_Type="circle"),
                            label="FINAL_",runs=20,flag_gif2=TRUE,
                            label_exp ="gradientSimple_circle")
  }
}

for( a in EXP_disp){
  results0 <- runExp(max_dispersal_ability=a,
                     background_object= generateBG(res=EXP_res,
                                                   type="uniform",
                                                   toRatio=1),
                     label="FINAL_",runs=20,flag_gif2=TRUE,
                     label_exp ="simple",
                     dispersalMode = "active")
  
  for(one_ratio in c(1,0.1,5,8,10) ){
    print(paste("this disp: ",a ,";   this ratio:",one_ratio) )
    results0 <- runExp(max_dispersal_ability=a,
                       background_object= generateBG(res=EXP_res,
                                                     type="gradientSimple",
                                                     toRatio=one_ratio),
                       label="FINAL_",runs=20,flag_gif2=TRUE,
                       label_exp ="gradientSimple",
                       dispersalMode = "active")
  }
}



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

