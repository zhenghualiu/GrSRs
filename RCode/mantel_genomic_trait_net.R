mantel_genomic_trait_net <- function(Metadata,NetMatrix,used.traits,group.trait=NULL,
                                     method.1="euclidean",method.2="gower",
                                     Zscores=F)
{
  library(plyr)
  if(is.null(group.trait))
  {
    traits.m <- Metadata[,used.traits]
    traits.m[is.na(traits.m)] <- 0
    del_rows <- which(rowSums(traits.m) == 0)
    if(length(del_rows) !=0) {
      traits.m <- traits.m[-del_rows,]
      NetMatrix <- NetMatrix[-del_rows,]
    }
    res <- c(r=NA,p=NA,N=NA)
   if(nrow(traits.m)>=10){
     if(Zscores) traits.m <- zscores(traits.m)
     tmp <- mantel(vegdist(traits.m,method=method.1),
     vegdist(NetMatrix,method=method.2,binary = T) )
   res <- c(r=tmp$statistic,p=tmp$signif,N=nrow(traits.m))
   }
  }
  
  if(!is.null(group.trait))
  {
    group_by_trait <- Metadata[,group.trait]
    group_by_trait.uni <- unique(group_by_trait)
    group_by_trait.uni <- group_by_trait.uni[!is.na(group_by_trait.uni)]
    res <- lapply(group_by_trait.uni,function(x){
      traits.m <- Metadata[group_by_trait == x,used.traits]
      traits.m[is.na(traits.m)] <- 0
      del_rows <- which(rowSums(traits.m) == 0)
      NetMatrix2 <- NetMatrix[rownames(traits.m),]
      if(length(del_rows) !=0) {
        traits.m <- traits.m[-del_rows,]
        NetMatrix2 <- NetMatrix2[-del_rows,]
      }
      tmp2 <- c(r=NA,p=NA,N=NA)
      if(nrow(traits.m)>=10){
      if(Zscores) traits.m <- zscores(traits.m)
      tmp <- mantel(vegdist(traits.m,method=method.1),
                    vegdist(NetMatrix2,method=method.2,binary = T) )
      tmp2 <- c(r=tmp$statistic,p=tmp$signif,N=nrow(traits.m))}
      return(tmp2)
    })
    names(res) <- group_by_trait.uni
    res <- ldply(res,.id ="Groups")
  }
  return(res)
  
}

zscores <- function(M){
  res <- apply(M,2,function(x){
    res <- (x-min(x))/(max(x)-min(x))})
  return(res)
}
