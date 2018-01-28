# initdat = initial data matrix (derived from indChar); disparams = disease parameter list (from disParams); days = number of days to project;
# n = population size

OutbreakR.init <- function(initdat, disparams, days, n) {
  ## transmission state projection matrix
  projdat <- as.data.frame(matrix(0,nrow=n,ncol=days+1))
  colnames(projdat)[1] <- "ST"
  colnames(projdat)[2:(days+1)] <- paste("D",seq(1:days),sep="")
  projdat[,1] <- initdat$dis
  
  # pre-susceptibles remaining so indefinitely
  nP <- which(projdat$ST == "P")
  lnP <- length(nP)
  n.perm.P <- ceiling(lnP*disparams[1])
  perm.P <- nP[sort(ceiling(runif(n.perm.P, min=1, max=lnP)))]
  if(length(perm.P) > 0) {
    projdat[perm.P, (2:(days+1))] <- rep("P",days)} # end if
  
  # time P->S
  nP2 <- which(projdat$ST == "P")
  lnP2 <- length(nP2)
  P.max.tim <- ifelse(disparams[3] > days, days, disparams[3])
  P.dur.tim <- ceiling(runif(lnP2,disparams[2], P.max.tim))
  P.dur.mat <- as.data.frame(matrix(0,nrow=lnP2,ncol=P.max.tim))
  for (i in 1:lnP2) {
    P.dur.vec <- c(rep("P",(P.dur.tim[i]-1)), "S")
    if (length(P.dur.vec) == 0) {
      P.dur.vec <- "S"} # end if
    P.dur.mat[i,1:(length(P.dur.vec))] <- P.dur.vec} # end i loop
  projdat[nP2, (2:(P.max.tim+1))] <- P.dur.mat
  
  # remaining I indefinitely
  nI <- which(projdat$ST == "I")
  lnI <- length(nI)
  if (ceiling(disparams[10]*lnI) > 0) {
    perm.I <- sort(sample(nI,ceiling(disparams[10]*lnI), replace=F))
    projdat[perm.I, (2:(days+1))] <- rep("I",days)
  }
  
  # S->E
  nS <- which(projdat$ST == "S")
  lnS <- length(nS)
  nSexp <- ceiling(lnI * disparams[4] * lnS) # number of susceptibles exposed
  nE.int <- ceiling(nSexp * disparams[5]) # new exposed from internal exposure
  nE.ext <- ceiling(lnS * disparams[6] * disparams[7]) # new exposed from outside source
  nE.new <- nE.int + nE.ext # total number of new exposed
  if (nE.new > 0) {
    E.new <- sort(sample(nS, nE.new, replace=F))
    projdat$D1[E.new] <- "E" # newly exposed on DAY 1
    E1.dur1 <- ceiling(runif(nE.new, disparams[8], disparams[9]))
    E1.dur <- ifelse(E1.dur1 > days, days, E1.dur1)
    E1.max.tim <- ifelse(disparams[9] > days, days, disparams[9])
    E1.dur.mat <- as.data.frame(matrix(0,nrow=nE.new,ncol=E1.max.tim))
    for (i in 1:nE.new) {
      E1.dur.vec <- rep("E",E1.dur[i])
      if (E1.dur[i] < days) {
        E1.dur.vec[E1.dur[i]] <- "I"} # end if
      if (length(E1.dur.vec) == 0) {
        E1.dur.vec <- "E"} # end if
      E1.dur.mat[i,1:(length(E1.dur.vec))] <- E1.dur.vec } # end for loop
    projdat[E.new, (3:(E1.max.tim+1))] <- E1.dur.mat[-1]} # end if
  
  # assign remaining individuals on DAY1 as status on ST
  D1.0 <- which(projdat$D1 == 0)
  projdat[D1.0, 2] <- as.character(projdat[D1.0, 1])
  
  # time remaining infectious
  nI2 <- which(projdat$D1 == "I")
  lnI2 <- length(nI2)
  I.max.tim <- ifelse(disparams[12] > days, days, disparams[12])
  I.dur.tim1 <- ceiling(runif(lnI2, disparams[11], disparams[12]))
  I.dur.tim <- ifelse(I.dur.tim1 > days, (days), I.dur.tim1)
  I.days.rem1 <- I.dur.tim1 - (days - 1)
  I.days.rem <- ifelse(I.days.rem1 < 0, 0, I.days.rem1)
  I.dur.mat <- as.data.frame(matrix(0,nrow=lnI2,ncol=I.max.tim))
  
  for (i in 1:lnI2) {
    I.dur.vec <- rep("I",I.dur.tim[i])
    if (length(I.dur.vec) == 0) {
      I.dur.vec <- "I"} # end if
    I.dur.mat[i,1:(length(I.dur.vec))] <- I.dur.vec} # end i loop
  projdat[nI2, (2:(I.max.tim+1))] <- I.dur.mat
  
  ## project through all days
  for (d in 3:(dim(projdat)[2])) {
    nS <- which(projdat[,(d-1)] == "S" & projdat[,d] == 0)
    lnS <- length(nS)
    nSexp <- ceiling(lnI* disparams[4] * lnS) # number of susceptibles exposed
    nE.int <- ceiling(nSexp * disparams[5]) # new exposed from internal exposure
    nE.ext <- ceiling(lnS * disparams[6] * disparams[7]) # new exposed from outside source
    nE.new <- nE.int + nE.ext # total number of new exposed
    if (nE.new > 0) {
      E.new <- sort(sample(nS, nE.new, replace=F))
      projdat[E.new, d] <- "E"} # newly exposed on DAY d
    
    if (nE.new > 0) {
      E1.dur1 <- ceiling(runif(nE.new, disparams[8], disparams[9]))
      E1.dur <- ifelse(E1.dur1 > (days+1-d), (days+1-d), E1.dur1)
      for (i in 1:nE.new) {
        E1.dur.vec <- rep("E", E1.dur[i])
        E1.dur.vec[length(E1.dur.vec)] <- "I"
        if (length(E1.dur.vec) < (days+1-d)) {
          projdat[E.new[i], ((d+1):(E1.dur[i]+d))] <- E1.dur.vec} # end if
        if (length(E1.dur.vec) == 0) {
          projdat[E.new[i], dim(projdat)[2]] <- "E"} # end if
      } # i loop
    }  # end if
    
    # assign remaining individuals on DAY d as status on d-1
    D.0 <- which(projdat[d] == 0)
    if (d <= dim(projdat)[2]) {
      projdat[D.0, d] <- as.character(projdat[D.0, (d-1)])} # end if
    print(paste("day = ", (d-1), sep=""))  
  } # end d loop
  
  # resistance/recovery
  # start recovery
  proj.tmpI <- as.data.frame(ifelse(projdat == "I", 1, 0))
  I.sub <- (ifelse(apply(proj.tmpI, 1, sum) > 0, 1, 0))
  ind.vec <- 1:n
  I.tot <- ind.vec[(which(I.sub == 1))]
  nI.tot <- length(I.tot) # total number of individuals infected
  nRpt <- ceiling(nI.tot*disparams[13]) # number of potential recovereds
  if (nRpt > 0) {
    Rpt <- sort(sample(I.tot, nRpt, replace=F)) # identify potential recovereds
    nI.days <- I.1st <- R.min.start <- R.start <- R.start2 <- rep(0,nRpt)
    for (r in 1:nRpt) {
      nI.days[r] <- length(which(projdat[Rpt[r], 2:(dim(projdat)[2])] == "I")) # number of days infected over total projection interval
      I.1st[r] <- which(projdat[Rpt[r], 1:(dim(projdat)[2])] == "I")[1] # day when first infected
      R.min.start[r] <- ifelse((I.1st[r] + disparams[11]) > dim(projdat)[2], dim(projdat)[2], (I.1st[r] + disparams[11])) # minimum day recovery could start
      R.start[r] <- R.min.start[r] + ceiling(runif(1, R.min.start[r], dim(projdat)[2])) # actual day of recovery
      R.start2[r] <- ifelse(R.start[r] > dim(projdat)[2], dim(projdat)[2], R.start[r])
      projdat[Rpt[r], R.start2[r]] <- "R"} # end r loop
    R.start.dat <- data.frame(Rpt, R.start2)
    
    # fill R days after start R day
    R.dur <- ceiling(runif(nRpt, disparams[16], disparams[17])) # duration of recovered time between min and max
    R.end <- R.start.dat[,2] + R.dur # end day of recovery
    R.remain <- ifelse((R.end - dim(projdat)[2]) < 0, 0, (R.end - dim(projdat)[2])) # days remaining infected beyond final projection day (holdover for next run)
    R.fin <- ifelse(R.end > dim(projdat)[2], dim(projdat)[2], R.end)
    for (s in 1:nRpt) {
      projdat[Rpt[s], (R.start.dat[s,2]:R.fin[s])] <- "R"
      if (R.fin[s] == dim(projdat)[2]) {
        projdat[Rpt[s], dim(projdat)[2]] <- projdat[Rpt[s], (dim(projdat)[2]-1)]
      } # end if
    } # end s loop
    
    # permanently resistant
    R.remain.dat <- data.frame(R.start.dat[,1], R.remain)
    colnames(R.remain.dat) <- c("recovered","R.remain")
    if (length(which(R.remain == 0)) > 0) {
      R.remain.nonzero <- R.remain.dat[-which(R.remain == 0), ]} # end if
    if (length(which(R.remain == 0)) == 0) {
      R.remain.nonzero <- R.remain.dat} # end if
    nR.perm <- round((dim(R.remain.nonzero)[1]*disparams[15]), 0)
    if (nR.perm > 0) {
      R.perm <- sort(sample(1:dim(R.remain.nonzero)[1], nR.perm, dim(R.remain.nonzero)[1]))
      R.remain.nonzero[R.perm,2] <- (disparams[18]+1)*365 # set time remaining to full lifespan (in days)
      R.perm.ind <- R.remain.nonzero[R.perm,1]}
    
    # returning to susceptible/pre-susceptible
    R.fin.dat <- data.frame(R.start.dat[,1], R.fin)
    colnames(R.fin.dat) <- c("recovered", "endR")
    R.end.before <- which(R.fin.dat[,2] < dim(projdat)[2])
    if (length(R.end.before) > 0) {
      endR.before <- R.fin.dat[R.end.before,] # which ones recovered before end of entire projection interval?
      S.return <- sort(sample(endR.before[,1], ceiling(dim(endR.before)[1]*disparams[14]))) # which individuals become susceptible?
      S.return.dat <- endR.before[endR.before[,1] %in% S.return,]
      P.return <- endR.before[!(endR.before[,1] %in% S.return),1] # which become pre-susceptible?
      P.return.dat <- endR.before[!(endR.before[,1] %in% S.return),]
      if (dim(S.return.dat)[1] > 0) {
        for (q in 1:length(S.return)) {
          projdat[S.return.dat[q,1], (S.return.dat[q,2]+1)] <- "S"
          projdat[S.return.dat[q,1], (S.return.dat[q,2]+2):dim(projdat)[2]] <- 0} # end if
      } # end q loop
      if (length(P.return) > 0) {
        for (p in 1:length(P.return)) {
          projdat[P.return.dat[p,1], (P.return.dat[p,2]+1)] <- "P"
          projdat[P.return.dat[p,1], (P.return.dat[p,2]+2):dim(projdat)[2]] <- 0
        } # end p loop
      } # end if
    
    # do any newly pre-susceptible individuals post-recovery become re-susceptible?
      if (length(P.return) > 0) {
        P.new.dur <- ceiling(runif(length(P.return), disparams[2], disparams[3]))
        for (p in 1:length(P.return)) {
          if ((P.return.dat[p,2]+1+P.new.dur[p]) <= dim(projdat)[2]) {
            projdat[P.return.dat[p,1], (P.return.dat[p,2]+1):(P.return.dat[p,2]+1+P.new.dur[p])] <- "S"
            projdat[P.return.dat[p,1], (P.return.dat[p,2]+2+P.new.dur[p]):dim(projdat)[2]] <- 0} # end if
          if ((P.return.dat[p,2]+1+P.new.dur[p]) > dim(projdat)[2]) {
            projdat[P.return.dat[p,1], (P.return.dat[p,2]+1):dim(projdat)[2]] <- "P"} # end if
        } # end p loop
      } # end if
    
      # do any newly susceptible individuals post-recovery become re-exposed?
      projdat2 <- projdat[endR.before[,1],]
    } # end if
    
    if (length(R.end.before) == 0) {
      projdat2 <- projdat}
    
    ## now reproject through days with zeros
    itdiv <- n/100
    for (s in 1:dim(projdat2)[1]) {
      dat.tmp <- projdat2[s, ]
      ind.id <- as.numeric(rownames(dat.tmp))
      start.day <- which(dat.tmp == 0)[1]
      if (is.na(start.day) == F) {
        start.state <- as.character(dat.tmp[(start.day - 1)])
        d.len <- dim(projdat)[2] - start.day + 1
        
        for (d in 1:d.len) {
          today <- start.day + d - 1
          nS.onday <- length(which((projdat[,(today-1)]) == "S")) # number of susceptibles in population on the day before
          nI.onday <- length(which((projdat[,(today-1)]) == "I")) # number of infecteds in population on the day before
          nSexp.onday <- ceiling(ceiling(nI.onday * disparams[4]) * (nS.onday/n))
          pr.this.ind.E <- (nSexp.onday/nS.onday * disparams[5]) + (nS.onday * disparams[6] * disparams[7]) # probability this individual is exposed
          state.this.ind <- ifelse((rbinom(1,1,pr.this.ind.E)) == 1, "E", "S") # is this individual exposed?  
          if (state.this.ind == "S") {
            projdat[ind.id, today] <- "S"} # end if
          
          if (state.this.ind == "E") {
            E1.dur1 <- ceiling(runif(1, disparams[8], disparams[9]))
            E1.dur <- ifelse(E1.dur1 > (days+1-today), (days+1-today), E1.dur1)
            E1.dur.vec <- rep("E", E1.dur)
            E1.dur.vec[length(E1.dur.vec)] <- "I"
            if (length(E1.dur.vec) < (days+1-today)) {
              projdat[ind.id, (today+1):(E1.dur+today)] <- E1.dur.vec} # end if
            if (length(E1.dur.vec) == 0) {
              projdat[ind.id, dim(projdat)[2]] <- "E"} # end if
          } # end if
        } # end d loop
      } # end if
      if (s %% itdiv==0) print(paste("processed ", s, " of ", n, " individuals", sep=""))
    } # end s loop
  } # end if
  
  # accumulate individual hold-over characteristics for next epidemiologial projection
  indholdoverdat <- initdat
  indholdoverdat$permP <- 0
  indholdoverdat[perm.P,dim(indholdoverdat)[2]] <- 1 # permanently pre-susceptible
  indholdoverdat$permI <- 0
  indholdoverdat[perm.I,dim(indholdoverdat)[2]] <- 1 # permanently infectious
  indholdoverdat$days.remI <- 0
  indholdoverdat[nI2,dim(indholdoverdat)[2]] <- I.days.rem # days remaining infectious for next projection interval set
  if (nR.perm > 0) {
    indholdoverdat$perm.R <- 0
    indholdoverdat$perm.R[R.perm.ind] <- 1 # permanently resistant
    indholdoverdat$days.remR <- 0
    indholdoverdat[R.remain.nonzero[,1], dim(indholdoverdat)[2]] <- R.remain.nonzero[,2] # days remaining recovered for next projection interval set
  } # end if
  indholdoverdat$end.dis <- projdat[,dim(projdat)[2]]
  
  dis.freq.final <- table(projdat[,dim(projdat)[2]])/n
  outlist <- list(projdat,indholdoverdat,dis.freq.final)
  print("final proportions per disease state")
  print(outlist[[3]]) # final disease frequencies at end of total projection interval

  ## timeline plots
  par(mfrow=c(2,3))
  # proportion pre-susceptible
  proj.tmp1 <- ifelse(outlist[[1]] == "P", 1, 0)
  TnS <- apply(proj.tmp1, 2, sum) / n
  plot(1:(days+1), TnS, type="l", xlab="days", ylab="proportion pre-susceptible")
  
  # proportion susceptible
  proj.tmp1 <- ifelse(outlist[[1]] == "S", 1, 0)
  TnS <- apply(proj.tmp1, 2, sum) / n
  plot(1:(days+1), TnS, type="l", xlab="days", ylab="proportion susceptible")
  
  # proportion susceptible
  proj.tmp2 <- ifelse(outlist[[1]] == "E", 1, 0)
  TnE <- apply(proj.tmp2, 2, sum) / n
  plot(1:(days+1), TnE, type="l", xlab="days", ylab="proportion exposed")
  
  # proportion infected
  proj.tmp3 <- ifelse(outlist[[1]] == "I", 1, 0)
  TnI <- apply(proj.tmp3, 2, sum) / n
  plot(1:(days+1), TnI, type="l", xlab="days", ylab="proportion infected")
  
  # proportion recovered
  proj.tmp4 <- ifelse(outlist[[1]] == "R", 1, 0)
  TnR <- apply(proj.tmp4, 2, sum) / n
  plot(1:(days+1), TnR, type="l", xlab="days", ylab="proportion recovered") # end if
  par(mfrow=c(1,1))
  
  return(outlist)
  
} # end OutbreakR.init function