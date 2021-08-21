# growth stree
#Temperature niche
#Tem, temperature at time t
#T0, optimal growth temperature
#sd0, niche breadth of temperature 
#Number of prophage in host
#Maximal phages contributing to host niche breadth

stress1 <- function(Tem,T0,sd0)
{
  str = exp( -(Tem - T0)^2 / (2 * sd0^2 ) ) 
  return(str)
}

stress2 <- function(Tem,n,T0,sd0,nc)
{
  k = sapply(n,function(x){min(nc,x)})
  CS = k_f(nc)
  k = CS[k]
  str2 = exp( -(Tem - T0)^2 / (2 * (sd0 + k)^2 ) )
  return(str2)
}
