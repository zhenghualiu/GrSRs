cal_lm_p <- function(y,x)
{
  tmp <- lm(y~x)
  tmp <- summary(tmp)
  try(N <- min(sum(!is.na(y)),sum(!is.na(x))))
  res <- c(tmp$coefficients[2,1],
           tmp$coefficients[2,4],
		   tmp$coefficients[2,2],
           tmp$adj.r.squared,
		   N
           )
  names(res) <- c("Coef","P","StdErr","Adj_R","N")
  return(res)
}