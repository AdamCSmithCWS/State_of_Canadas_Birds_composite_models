specid_rename <- function(x){
  y <- vector("character",length(x))
  for(i in 1:length(x)){
    if(grepl("speciesId",x[i]) |
       grepl("species_id",x[i])){
      y[i] <- "speciesID"
    }else{
      y[i] <- x[i]
    }
  }
  return(y)
}
