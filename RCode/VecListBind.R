VecListBind <- function(x)
{
  Ln <- names(x)
  if(is.null(Ln)) Ln <- seq(1,length(x))
  ObjName <- unique(unlist(lapply(x,names)))
  res <- as.data.frame(matrix(NA,nrow = length(ObjName),ncol=length(x)))
  rownames(res) <- ObjName
  for(i in 1:length(x))
  {
    res[names(x[[i]]),i] <- x[[i]]
  }
  colnames(res) <- Ln
  return(res)
}

VecListBind2 <- function(x)
{
  Ln <- names(x)
  if(is.null(Ln)) Ln <- seq(1,length(x))
  ObjName <- unique(unlist(x))
  res <- as.data.frame(matrix(NA,nrow = length(x),ncol=length(ObjName)))
  rownames(res) <- Ln
  colnames(res) <- ObjName
  for(i in 1:nrow(res))
  {
   tmp <- table(x[[i]])
   res[i,names(tmp)] <- as.numeric(tmp)
  }
  return(res)
}