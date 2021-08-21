#Growth rate 
#N, population size
#Vm, maximal growth rate
#C, environment capacity
gr <- function(N,Vm,C) 
{
  growth_rate = Vm * (1 - N/C )
  return(growth_rate)
}
