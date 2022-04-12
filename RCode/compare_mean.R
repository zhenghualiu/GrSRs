compare_mean <- function(y,x,method="kruskal.test")
{
  if(method=="kruskal.test") tmp <- kruskal.test(y~x)
  res <- tmp$p.value
  names(res) <- "P"
  return(res)
}