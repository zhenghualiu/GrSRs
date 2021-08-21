#Immunity rate
#Tem, temperature at time t
#K, average of phage infecting one host
#T1, optimal temperature of immunity activity
#sd1, temperature breadth of immunity activity
#K0, defensive capacity

I <- function(Tem,Ic,T1,sd1)
{
  immunity = Ic * exp(-(Tem-T1)^2/(2*sd1^2))
  return(immunity)
}
