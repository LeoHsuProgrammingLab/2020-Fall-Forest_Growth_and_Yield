
##################
### DBH growth ###
##################


### DBH increment
a0 <- -6.801 ; a1 <- 0.213 ; a2 <- -0.002 ; a3 <- -0.030 ; a4 <- 0.195 ; rmse <- 0.203
DBH_inc <- function(DBH,Age,H,D){
  mH <- mean(H)
  CI <- 10000/(mH*D^0.5)
  lnG <- a0 + a1*DBH + a2*DBH^2 + a3*Age + a4*CI + rnorm(1,0,rmse)
  return(exp(lnG))
}
### Mortality
mor_rate <- function(DBH_this,DBH_next){
  mor <- 1-(DBH_next/DBH_this)^(-1.19+rnorm(1,0,sqrt(0.1)))
  return(mor)
}
### Height increment
All_k <- 0.0001*c(336,351,378,376,339,558,460,
                  331,454,574,453,449,551)
k_N <- mean(All_k)
###M
All_M <- 0.01*c(2671,2516,2542,2601,2525,2384,2357,
                2693,2394,2253,2483,2269,2294)
M_N <- mean(All_M)
#MSE for Height Prsdiction (Wang et al., 2012)
##Use result from Model 4
MSE_W <- 1.64
#Height Growth Function
H_Growth <- function(H, Age, k = k_N, M = M_N, MSE = MSE_W){
  H1 <- M*((1-exp(-k*(Age+1)))) #去掉rnorm
  Growth <- H1 - H
  return(Growth)
}

#DBH=test$dbh.cm;H=mean(test$height);D=nrow(test);startage=20;endage=100
DBH = test$DBH; H = mean(test$H); D = nrow(test); startage = 10; endage = 100; Area = 0.153
library(tidyr)
library(dplyr)

CjGY <- function(DBH,H,D,startage,endage, Area){
  
  forest <- data.frame("ID"=seq(1:D),"DBH"=DBH,"Height"=H,"Age"=rep(startage,D))
  stand_attribute <- data.frame("Age"=startage,"QMD"=sqrt(mean(DBH^2, na.rm = T)),"Stand_height"=H,"Den"=D/Area,"Stand_vol"=(sum(0.00007854*DBH^2, na.rm = T)*H*0.45))
  
  H_mean_now <- mean(H)
  DBH_now <- DBH
  D_now <- D
  
  for (i in startage:(endage-1)) {
    #Calculating Growth
    DBH_next <- DBH_now + DBH_inc(DBH = DBH_now, H = H_mean_now, D = D_now/Area, Age = i)
    H_mean_next <- H_mean_now + H_Growth(H = H_mean_now, Age = i)
    MorRate_now <- mor_rate(DBH_this = DBH_now, DBH_next = DBH_next)
    Death_index <- runif(D)
    Death_boolean <- Death_index <= MorRate_now
    DBH_next[Death_boolean] <- NA
    D_next <- D_now - sum(Death_boolean, na.rm = T)
    
    #append data
    df <- data.frame("ID"=1:D, "DBH"=DBH_next, "Height"=rep(H_mean_next, D), "Age"=rep(i+1, D))
    forest <- forest %>% bind_rows(df)
    stand_attribute[nrow(stand_attribute)+1,] <- c(i+1, sqrt(mean(DBH_next^2, na.rm = T)), H_mean_next, D_next/Area, sum(0.00007854*DBH_next^2, na.rm = T)*H_mean_next*0.45)
    
    #update
    H_mean_now <- H_mean_next
    DBH_now <- DBH_next
    D_now <- D_next
    
    if (D_now <= 0) {
      print('stop because nonpositive stand density')
      break
    }
  }
  #output
  return(list(forest, stand_attribute))
}

test <- read.csv("C:/D/ForestGandY/JGZ1936.csv")

run_test <- CjGY(test$DBH, mean(test$H), nrow(test), 10, 100, 0.153)

run_test[1]
plot(run_test[[2]]$Age, run_test[[2]]$QMD)
