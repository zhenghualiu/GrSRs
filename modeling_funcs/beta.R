# Burst size, particles
# n, number of prophage in host
# Basal burst size, particles

beta <- function(n,beta0)
{
  Beta <- n * beta0 
  return(Beta)
}