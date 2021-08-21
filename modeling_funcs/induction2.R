#Induction rate2

cauchy2 <- function(Tem,cauchy0,sc,T2,T3)
{
  f1 = sigmod_f(Tem,sc,T2,"heat")
  f2 = sigmod_f(Tem,sc,T3,"cold")
  Cauchy = cauchy0 * max(f1,f2)
  return(Cauchy)
}

sigmod_f <- function(Tem,sc,Tc,type)
{
  if(type == "heat") res = 1/(1+exp(-sc * (Tem - Tc) ))
  if(type == "cold") res = 1/(1+exp(-sc * (Tc - Tem) ))
  if(res >= 0.999) res = 1
  if(res <= 0.001) res = 0
  return(res)
}