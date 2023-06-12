#Height Growth Equation (T. Nishizono et al., 2014)
##Calculating Parameters (Take means for equation parameters estimated for different sites)
###k
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
H_Growth <- function(H0, EndAge, k = k_N, M = M_N, MSE = MSE_W){
  H <- M*(1-exp(-k*EndAge)) + rnorm(1, mean = 0, sd = sqrt(MSE))
  Growth <- H - H0
  return(Growth)
}

#Example
H_Growth(10, 11)
