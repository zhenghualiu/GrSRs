#Temperature niche
#Tem, temperature at time t
#T0, optimal growth temperature
#sd0, niche breadth of temperature 
#Number of prophage in host
#Maximal phages contributing to host niche breadth

#For Bacteria 
theta1 <- function(Tem,T0,sd0)
{
  Theta1 = exp( -(Tem - T0)^2 / (2 * sd0^2 ) )
  return(Theta1)
}

#For lysogeny
theta2 <- function(Tem,n,T0,sd0,nc)
{
  k <- sapply(n,function(x){min(nc,x)})
  CS <- k_f(nc)
  k <- CS[k]
  Theta2 = exp( -(Tem - T0)^2 / (2 * (sd0 + k)^2 ) )
  return(Theta2)
}

k_f <- function(nc)
{
  l <- seq(1,nc)
  k <- 1-(l-1)/nc
  CS <- cumsum(k)
  return(CS)
}