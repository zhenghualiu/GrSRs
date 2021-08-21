#Induction rate
#Tem, temperature at time t
#gr, growth rate at time t
#cauchy0, induction rate constant
#T2, temperature switch of lytic at heat condition
#T3, temperature switch of lytic at cold condition
#f0, switch activity
#e, basal switch activity
#k, constant
#Tem=30;gr=0.1;cauchy0=4.1 * 10^(-8);T2=40;T3=0;f0=0.77;e=0.033;k=0.55;mu0=0.062;H = 1.8
#Induction 
cauchy <- function(Tem,gr,cauchy0,T2,T3,f0,e,k,mu0,H)
{
  f1 = f(Tem,gr,T2,f0,e,k,mu0,H,"heat")
  f2 = f(Tem,gr,T3,f0,e,k,mu0,H,"cold")
  Cauchy = cauchy0 * ( 1 - min(f1,f2) )
  return(Cauchy)
}

#Hill function
f <- function(Tem,gr,Tc,f0,e,k,mu0,H,type)
{
  f_ = f0 * (e + (1 - e) / (1 + (mu(Tem,gr,k,Tc,type)/mu0)^(-H)) )
  return(f_)
}

#switch activity at given temperature for heat or cold lysis
mu <- function(Tem,gr,k,Tc,type="heat")
{
  if(type == "heat") Mu = max(gr/(gr + exp( k * (Tem - Tc) )),0)
  if(type == "cold") Mu = max(gr/(gr + exp( k * (Tc - Tem) )),0)
  return(Mu)
}