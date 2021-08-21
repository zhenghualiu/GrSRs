#Induction rate
#Tem, temperature at time t
#deltaC, induction rate constant

delta <- function(Tem,deltaC)
{
  Delta = deltaC * Tem^(0.5)
  return(Delta)
}
