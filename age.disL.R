# set age-specific likelihoods of disease state (e.g., all juveniles "P"; randomise rest)
# a = age classes vector; d = disease states vector

age.disL <- function (a, d) {
  agedisLmat <- as.data.frame(matrix(1, nrow=length(a), ncol=length(d)))
  colnames(agedisLmat) <- c(d)
  rownames(agedisLmat) <- c(a)
  return(agedisLmat)} # end age.disL function