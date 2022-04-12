cal_cor <- function(y,x,method="pearson")
{
  res <- NA
  tmp <- NA
  try(tmp <- cor.test(y,x,method=method),silent =T)
  if(!is.na(tmp)){
  try(StdErr <- (tmp$conf.int[2] -mean(tmp$conf.int))/1.96)
  if(length(StdErr)==0) StdErr=0
  try(N <- min(sum(!is.na(y)),sum(!is.na(x))))
  try(res <- c(tmp$estimate,tmp$p.value,StdErr,N))
  try(names(res) <- c("Cor","P","StdErr","N"))
  }
  return(res)
}

cal_cor2 <- function(y,x,v="estimate")
{
  res <- NA
  tmp <- NA
  try(tmp <- cor.test(y,x),silent =T)
  if(!is.na(tmp)){res=eval(parse(text=paste("tmp$",v)))}
  return(res)
}
