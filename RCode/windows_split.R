windows_split <- function(x,N_col,wd)
{
  res <- NULL
  x[is.na(x[,N_col]),N_col] <- 0
  for(i in wd)
  {
    if(i==0){
      row_id <- which(x[,N_col] ==i)
      if(length(row_id) == 0) next
      tmp <- x[row_id,]
      tmp$cut_off <- i
      res <- rbind.data.frame(tmp,res)
    }else{
    tmp <- x[x[,N_col] >=i,]
    tmp$cut_off <- i
    res <- rbind.data.frame(tmp,res)
    }
  }
  return(res)
}