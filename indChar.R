# randomise initial individual characteristics
# a = age proportion vector; ac = age classes vector; s = sex ratio, d = initial disease proportions vector; ds = disease states vector
# n = population size, disL = age-specific likelihoods of disease state matrix
indChar <- function(a, as, s, d, ds, disL, n) {
  library(tcltk)
  age.dat <- as.data.frame(t(rmultinom(n, size = 1, prob = a)))
  age.vec <- rep(0, n)
  pb <- txtProgressBar(min=1, max=n, style=3)
  for (i in 1:n) {
    sub.age <- which(age.dat[i,] == 1)
    age.vec[i] <- age.classes[sub.age]
    setTxtProgressBar(pb, i)} # end i loop
  close(pb)
  sex.vec <- ifelse(rbinom(n,1,s) == 1, "F", "M")
  dis.state.dat <- as.data.frame(t(rmultinom(n, size = 1, prob = d)))
  dis.state.vec <- rep(0, n)
  pb2 <- txtProgressBar(min=1, max=n, style=3)
  for (i in 1:n) {
    sub.dis <- which(dis.state.dat[i,] == 1)
    dis.state.vec[i] <- disease.states[sub.dis]
    setTxtProgressBar(pb2, i)} # end i loop
  close(pb2)
  ind.vec <- 1:n
  initdat <- data.frame(ind.vec, sex.vec, age.vec, dis.state.vec)
  colnames(initdat) <- c("ind", "sex", "age", "dis")
  pb3 <- txtProgressBar(min=min(a), max=length(a), style=3)
  for (j in 1:length(a)) {
    dat.sub <- subset(initdat, age == as[j])
    dis.prop.update <- disL[j, ] * d / sum(disL[j, ] * d)
    dis.state.new <- as.data.frame(t(rmultinom(dim(dat.sub)[1], size = 1, prob = dis.prop.update)))
    dis.state.new.vec <- rep(0, dim(dat.sub)[1])
    for (i in 1:dim(dat.sub)[1]) {
      sub.dis1 <- which(dis.state.new[i,] == 1)
      dis.state.new.vec[i] <- ds[sub.dis1]} # end i loop
    initdat[dat.sub$ind,4] <- dis.state.new.vec
    setTxtProgressBar(pb3, j)} # end j loop
  close(pb3)
  return(initdat)
} # end indChar function
