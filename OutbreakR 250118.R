# Emerging Case: OutbreakR: an R library for epidemiological projections of wildlife and human diseases
# 
# Project
# Program Outbreak (www.vortex10.org/Outbreak.aspx) is an individual-based model software package (open access) that simulates disease dynamics
# using transitions among susceptible, exposed, infectious and recovered individuals. Outbreak provides several options for modes of transmission
# (random contact within populations, spatially based transmission, contact with environmental sources of disease) and provides options for
# management through vaccination or culling (written by Robert Lacy, JP Pollak, PS Miller, L Hungerford, and P Bright).
# 
# While Outbreak is ideal for modelling the dynamics of diseases in small populations, it could be limited for applications involving
# larger (>> 1000s) populations because of the detailed, individual-to-individual contacts invoked in the individual-based model
# at each transition interval. We therefore propose to write an R Package (cran.r-project.org) library that includes a large component
# of the main Outbreak methodology, but calculates disease-state transitions based on a weighted population of infected and susceptible
# individuals at each time step.
# 
# The main idea is to hybridise the individual- and population-based algorithms such that the disease dynamics of larger populations
# can be modelled. We propose to start with a detailed matrix of individuals in a population (rows), and their characteristics 
# (age, sex, condition, etc.) and initial disease state (pre-susceptible, susceptible, exposed, infected, recovered) at the beginning
# of the projection. For each interval, the disease state of each individual is assessed as the product of the probabilities of contact,
# transmission, etc., weighted by the number of individuals in the corresponding disease state (e.g., susceptibles can transition
# to ‘exposed’ as a product of the contact probability, the transmission probability, and the number of infected individuals in the entire population).
# 
# At the end of the desired number of disease-transition intervals (e.g., daily for 1 year), the frequency table of individuals
# in each disease category can then be parsed to other applications in R (e.g., matrix population models) or for external applications
# (e.g., RAMAS/Metapop). Over time, the addition of spatial elements, environmental correlates, and other complexities can be added 
# to the base functions in the library.
# 
# Meetings & Activities
# RCN Synthesis Meeting, White Oak Conservation, Yulee, Florida, USA (15-19 January 2018)
# 
# Participants: 
# Robert Lacy (Chicago Zoological Society); Phil Miller (Species Survival Commission: Species Planning Specialist Group); JP Pollak 
# (University of Cornell), Ana Davidson (Colorado State University); Nadia Ali (University of Chicago); Loren Cassin Sackett 
# (University of South Florida); Taylor Callicrate (Chicago Zoological Society; Species Conservation Toolkit Initiative)


#####################################################################################################
#####################################################################################################
# create dummy input matrix
pop.size <- 10000
age.classes <- c("juv","sub","adu")
age.prop <- c(0.50, 0.30, 0.20)
life.span <- 17 # in years
sex.classes <- c("F","M")
sex.ratio <- 0.50

disease.states <- c("P", "S", "E", "I",  "R")
init.dis.prop <- c(0.05, 0.9, 0.0, 0.05, 0.0)

# How many days to project?
days.proj <- 365 # number of days to project disease transition states into the future

# randomise initial individual characteristics
# disease state 
dis.state.dat <- as.data.frame(t(rmultinom(pop.size, size = 1, prob = init.dis.prop)))
dis.state.vec <- rep(0, pop.size)
for (i in 1:pop.size) {
  sub.dis <- which(dis.state.dat[i,] == 1)
  dis.state.vec[i] <- disease.states[sub.dis]
}

# age 
age.dat <- as.data.frame(t(rmultinom(pop.size, size = 1, prob = age.prop)))
age.vec <- rep(0, pop.size)
for (i in 1:pop.size) {
  sub.age <- which(age.dat[i,] == 1)
  age.vec[i] <- age.classes[sub.age]
}

# sex
sex.vec <- ifelse(rbinom(pop.size,1,sex.ratio) == 1, "F", "M")

# set up init data frame
ind.vec <- 1:pop.size
init.dat <- data.frame(ind.vec, sex.vec, age.vec, dis.state.vec)
colnames(init.dat) <- c("ind","sex","age","dis")

# disease parameters
encounter.rate <- 5 # individuals per day
prob.never.S <- 0.05 # probability pre-susceptible never becomes susceptible (permanent immunity)
min.tim.b.ind.S <- 1 # minimum time before an individual becomes susceptible
max.tim.b.ind.S <- 5 # maximum time before an individual becomes susceptible
contact.rate <- 0.1 # proportion of the population an exposed/infected individual encounters/day
daily.trans.pr <- 0.27 # daily transmission probability
out.dis.pr.trans <- 0.0027 # outside disease source probabilty of transmission
min.lat.inc.E <- 4 # minimum (days) incubation (latent) period of infection for exposed individuals
max.lat.inc.E <- 20 # maximum (days) incubation (latent) period of infection for exposed individuals
prop.I.remain <- 0.425 # proportion of infected individuals remaining so indefinitely
min.tim.rem.I <- 60 # minimum number of days an animal can remain infectious
max.tim.rem.I <- 17*365 # maximum number of days an animal can remain infectious
pr.R <- 0.01 # probability of recovery and subsequent resistance after the minimum time being infectious 
pr.ret.S <- 0.9 # probability of returning to the susceptible state after recovery (assuming those not returning to susceptible are pre-susceptible)
pr.R.perm.imm <- 0.01 # proportion of recovered individuals acquiring permanent immunity
min.tim.res <- 20 # minimum number of days resistant
max.tim.res <- 60 # maximum number of days resistant


## transmission state projection matrix
proj.dat <- as.data.frame(matrix(0,nrow=pop.size,ncol=days.proj+1))
colnames(proj.dat)[1] <- "ST"
colnames(proj.dat)[2:(days.proj+1)] <- paste("D",seq(1:days.proj),sep="")
proj.dat[,1] <- init.dat$dis

# project
# pre-susceptibles remaining so indefinitely
nP <- which(proj.dat$ST == "P")
lnP <- length(nP)
n.perm.P <- round(lnP*prob.never.S,0)
perm.P <- nP[sort(round(runif(n.perm.P, min=1, max=lnP),0))]
if(length(perm.P) > 0) {
  proj.dat[perm.P, (2:(days.proj+1))] <- rep("P",days.proj)}

# time P->S
nP2 <- which(proj.dat$ST == "P")
lnP2 <- length(nP2)
P.max.tim <- ifelse(max.tim.b.ind.S > days.proj, days.proj, max.tim.b.ind.S)
P.dur.tim <- round(runif(lnP2,min.tim.b.ind.S, P.max.tim),0)
P.dur.mat <- as.data.frame(matrix(0,nrow=lnP2,ncol=P.max.tim))
for (i in 1:lnP2) {
  P.dur.vec <- c(rep("P",(P.dur.tim[i]-1)), "S")
  if (length(P.dur.vec) == 0) {
    P.dur.vec <- "S"}
  P.dur.mat[i,1:(length(P.dur.vec))] <- P.dur.vec
}
proj.dat[nP2, (2:(P.max.tim+1))] <- P.dur.mat

# remaining I indefinitely
nI <- which(proj.dat$ST == "I")
lnI <- length(nI)
if (round(prop.I.remain*lnI, 0) > 0) {
    perm.I <- sort(sample(nI,round(prop.I.remain*lnI, 0), replace=F))
    proj.dat[perm.I, (2:(days.proj+1))] <- rep("I",days.proj)
}

# S->E
nS <- which(proj.dat$ST == "S")
lnS <- length(nS)
nSexp <- round(round(lnI*contact.rate, 0) * lnS/pop.size, 0) # number of susceptibles exposed
nE.int <- round(nSexp * daily.trans.pr, 0) # new exposed from internal exposure
nE.ext <- round(lnS * out.dis.pr.trans, 0) # new exposed from outside source
nE.new <- nE.int + nE.ext # total number of new exposed
if (nE.new > 0) {
  E.new <- sort(sample(nS, nE.new, replace=F))
  proj.dat$D1[E.new] <- "E"} # newly exposed on DAY 1

E1.dur1 <- round(runif(nE.new, min.lat.inc.E, max.lat.inc.E), 0)
E1.dur <- ifelse(E1.dur1 > days.proj, days.proj, E1.dur1)
E1.dur.arr <- sapply(E1.dur, rep, x="E")
for (i in 1:nE.new) {
  E1.dur.arr[[i]][E1.dur[i]] <- "I"
  proj.dat[E.new[i], (3:(E1.dur[i]+2))] <- E1.dur.arr[[i]]
}

# assign remaining individuals on DAY1 as status on ST
D1.0 <- which(proj.dat$D1 == 0)
proj.dat[D1.0, 2] <- as.character(proj.dat[D1.0, 1])


# time remaining infectious
nI2 <- which(proj.dat$D1 == "I")
lnI2 <- length(nI2)
I.max.tim <- ifelse(max.tim.rem.I > days.proj, days.proj, max.tim.rem.I)
I.dur.tim1 <- round(runif(lnI2, min.tim.rem.I, max.tim.rem.I),0)
I.dur.tim <- ifelse(I.dur.tim1 > days.proj, (days.proj), I.dur.tim1)
I.days.rem1 <- I.dur.tim1 - (days.proj - 1)
I.days.rem <- ifelse(I.days.rem1 < 0, 0, I.days.rem1)
I.dur.mat <- as.data.frame(matrix(0,nrow=lnI2,ncol=I.max.tim))

for (i in 1:lnI2) {
  I.dur.vec <- rep("I",I.dur.tim[i])
  if (length(I.dur.vec) == 0) {
    I.dur.vec <- "I"}
  I.dur.mat[i,1:(length(I.dur.vec))] <- I.dur.vec
}
proj.dat[nI2, (2:(I.max.tim+1))] <- I.dur.mat



## project through all days
for (d in 3:(dim(proj.dat)[2])) {

  nS <- which(proj.dat[,(d-1)] == "S" & proj.dat[,d] == 0)
  lnS <- length(nS)
  nSexp <- round(round(lnI*contact.rate, 0) * lnS/pop.size, 0) # number of susceptibles exposed
  nE.int <- round(nSexp * daily.trans.pr, 0) # new exposed from internal exposure
  nE.ext <- round(lnS * out.dis.pr.trans, 0) # new exposed from outside source
  nE.new <- nE.int + nE.ext # total number of new exposed
  if (nE.new > 0) {
    E.new <- sort(sample(nS, nE.new, replace=F))
    proj.dat[E.new, d] <- "E"} # newly exposed on DAY d
  
  if (nE.new > 0) {
    E1.dur1 <- round(runif(nE.new, min.lat.inc.E, max.lat.inc.E), 0)
    E1.dur <- ifelse(E1.dur1 > (days.proj+1-d), (days.proj+1-d), E1.dur1)
    for (i in 1:nE.new) {
      E1.dur.vec <- rep("E", E1.dur[i])
      E1.dur.vec[length(E1.dur.vec)] <- "I"
      if (length(E1.dur.vec) < (days.proj+1-d)) {
        proj.dat[E.new[i], ((d+1):(E1.dur[i]+d))] <- E1.dur.vec} # end if
      if (length(E1.dur.vec) == 0) {
        proj.dat[E.new[i], dim(proj.dat)[2]] <- "E"} # end if
    } # i loop
  }  # end if
  
  # assign remaining individuals on DAY d as status on d-1
  D.0 <- which(proj.dat[d] == 0)
  if (d <= dim(proj.dat)[2]) {
    proj.dat[D.0, d] <- as.character(proj.dat[D.0, (d-1)])} # end if
    
} # end d loop



# resistance/recovery
# start recovery
proj.tmpI <- as.data.frame(ifelse(proj.dat == "I", 1, 0))
I.sub <- (ifelse(apply(proj.tmpI, 1, sum) > 0, 1, 0))
I.tot <- ind.vec[(which(I.sub == 1))]
nI.tot <- length(I.tot) # total number of individuals infected
nRpt <- round(nI.tot*pr.R, 0) # number of potential recovereds
Rpt <- sort(sample(I.tot, nRpt, replace=F)) # identify potential recovereds
nI.days <- I.1st <- R.min.start <- R.start <- R.start2 <- rep(0,nRpt)
for (r in 1:nRpt) {
  nI.days[r] <- length(which(proj.dat[Rpt[r], 2:(dim(proj.dat)[2])] == "I")) # number of days infected over total projection interval
  I.1st[r] <- which(proj.dat[Rpt[r], 1:(dim(proj.dat)[2])] == "I")[1] # day when first infected
  R.min.start[r] <- ifelse((I.1st[r] + min.tim.rem.I) > dim(proj.dat)[2], dim(proj.dat)[2], (I.1st[r] + min.tim.rem.I)) # minimum day recovery could start
  R.start[r] <- R.min.start[r] + round(runif(1, R.min.start[r], dim(proj.dat)[2]), 0) # actual day of recovery
  R.start2[r] <- ifelse(R.start[r] > dim(proj.dat)[2], dim(proj.dat)[2], R.start[r])
  proj.dat[Rpt[r], R.start2[r]] <- "R"} # end r loop
R.start.dat <- data.frame(Rpt, R.start2)

# fill R days after start R day
R.dur <- round(runif(nRpt, min.tim.res, max.tim.res), 0) # duration of recovered time between min and max
R.end <- R.start.dat[,2] + R.dur # end day of recovery
R.remain <- ifelse((R.end - dim(proj.dat)[2]) < 0, 0, (R.end - dim(proj.dat)[2])) # days remaining infected beyond final projection day (holdover for next run)
R.fin <- ifelse(R.end > dim(proj.dat)[2], dim(proj.dat)[2], R.end)
for (s in 1:nRpt) {
  proj.dat[Rpt[s], (R.start.dat[s,2]:R.fin[s])] <- "R"
  if (R.fin[s] == dim(proj.dat)[2]) {
    proj.dat[Rpt[s], dim(proj.dat)[2]] <- proj.dat[Rpt[s], (dim(proj.dat)[2]-1)]
  }
} # end s loop


# permanently resistant
R.remain.dat <- data.frame(R.start.dat[,1], R.remain)
colnames(R.remain.dat) <- c("recovered","R.remain")
if (length(which(R.remain == 0)) > 0) {
  R.remain.nonzero <- R.remain.dat[-which(R.remain == 0), ]} # end if
if (length(which(R.remain == 0)) == 0) {
  R.remain.nonzero <- R.remain.dat} # end if
nR.perm <- round(dim(R.remain.nonzero)[1]*pr.R.perm.imm, 0)
R.perm <- sort(sample(1:dim(R.remain.nonzero)[1], nR.perm, dim(R.remain.nonzero)[1]))
R.remain.nonzero[R.perm,2] <- life.span*365 # set time remaining to full lifespan (in days)

# returning to susceptible/pre-susceptible
R.fin.dat <- data.frame(R.start.dat[,1], R.fin)
colnames(R.fin.dat) <- c("recovered", "endR")
R.end.before <- which(R.fin.dat[,2] < dim(proj.dat)[2])
if (length(R.end.before > 0)) {
  endR.before <- R.fin.dat[R.end.before,] # which ones recovered before end of entire projection interval?
  S.return <- sort(sample(endR.before[,1], round(dim(endR.before)[1]*pr.ret.S, 0))) # which individuals become susceptible?
  S.return.dat <- endR.before[endR.before[,1] %in% S.return,]
  P.return <- endR.before[!(endR.before[,1] %in% S.return),1] # which become pre-susceptible?
  P.return.dat <- endR.before[!(endR.before[,1] %in% S.return),]
  for (q in 1:length(S.return)) {
    proj.dat[S.return.dat[q,1], (S.return.dat[q,2]+1)] <- "S"
    proj.dat[S.return.dat[q,1], (S.return.dat[q,2]+2):dim(proj.dat)[2]] <- 0
  } # end q loop
  for (p in 1:length(P.return)) {
    proj.dat[P.return.dat[p,1], (P.return.dat[p,2]+1)] <- "P"
    proj.dat[P.return.dat[p,1], (P.return.dat[p,2]+2):dim(proj.dat)[2]] <- 0
  } # end p loop
} # end if


# do any newly pre-susceptible individuals post-recovery become re-susceptible?
P.new.dur <- round(runif(length(P.return), min.tim.b.ind.S, max.tim.b.ind.S), 0)
for (p in 1:length(P.return)) {
  if ((P.return.dat[p,2]+1+P.new.dur[p]) <= dim(proj.dat)[2]) {
    proj.dat[P.return.dat[p,1], (P.return.dat[p,2]+1):(P.return.dat[p,2]+1+P.new.dur[p])] <- "S"
    proj.dat[P.return.dat[p,1], (P.return.dat[p,2]+2+P.new.dur[p]):dim(proj.dat)[2]] <- 0} # end if
  if ((P.return.dat[p,2]+1+P.new.dur[p]) > dim(proj.dat)[2]) {
    proj.dat[P.return.dat[p,1], (P.return.dat[p,2]+1):dim(proj.dat)[2]] <- "P"} # end if
} # end p loop

# do any newly susceptible individuals post-recovery become re-exposed?
proj.dat2 <- proj.dat[endR.before[,1],]

## now reproject through days with zeros
for (s in 1:dim(proj.dat2)[1]) {
  dat.tmp <- proj.dat2[s, ]
  ind.id <- as.numeric(rownames(dat.tmp))
  start.day <- which(dat.tmp == 0)[1]
  if (is.na(start.day) == F) {
    start.state <- as.character(dat.tmp[(start.day - 1)])
    d.len <- dim(proj.dat)[2] - start.day + 1
    
    for (d in 1:d.len) {
      today <- start.day + d - 1
      nS.onday <- length(which((proj.dat[,(today-1)]) == "S")) # number of susceptibles in population on the day before
      nI.onday <- length(which((proj.dat[,(today-1)]) == "I")) # number of infecteds in population on the day before
      nSexp.onday <- round(round((nI.onday * contact.rate), 0) * (nS.onday/pop.size), 0)
      pr.this.ind.E <- nSexp.onday/nS.onday * (daily.trans.pr+out.dis.pr.trans) # probability this individual is exposed
      state.this.ind <- ifelse((rbinom(1,1,pr.this.ind.E)) == 1, "E", "S") # is this individual exposed?  
      if (state.this.ind == "S") {
        proj.dat[ind.id, today] <- "S"} # end if
      
      if (state.this.ind == "E") {
        E1.dur1 <- round(runif(1, min.lat.inc.E, max.lat.inc.E), 0)
        E1.dur <- ifelse(E1.dur1 > (days.proj+1-today), (days.proj+1-today), E1.dur1)
        E1.dur.vec <- rep("E", E1.dur)
        E1.dur.vec[length(E1.dur.vec)] <- "I"
          if (length(E1.dur.vec) < (days.proj+1-today)) {
            proj.dat[ind.id, (today+1):(E1.dur+today)] <- E1.dur.vec} # end if
          if (length(E1.dur.vec) == 0) {
            proj.dat[ind.id, dim(proj.dat)[2]] <- "E"} # end if
      } # end if
    } # end d loop
  } # end if
} # end s loop


# accumulate individual hold-over characteristics for next epidemiologial projection
ind.holdover.dat <- init.dat
ind.holdover.dat$permP <- 0
ind.holdover.dat[perm.P,dim(ind.holdover.dat)[2]] <- 1 # permanently pre-susceptible
ind.holdover.dat$permI <- 0
ind.holdover.dat[perm.I,dim(ind.holdover.dat)[2]] <- 1 # permanently infectious
ind.holdover.dat$days.remI <- 0
ind.holdover.dat[nI2,dim(ind.holdover.dat)[2]] <- I.days.rem # days remaining infectious for next projection interval set
ind.holdover.dat$perm.R <- 0
ind.holdover.dat$days.remR <- 0
ind.holdover.dat[R.remain.nonzero[,1], dim(ind.holdover.dat)[2]] <- R.remain.nonzero[,2] # days remaining recovered for next projection interval set
ind.holdover.dat$end.dis <- proj.dat[,dim(proj.dat)[2]]


## timeline plots
par(mfrow=c(2,2))
# proportion susceptible
proj.tmp1 <- ifelse(proj.dat == "S", 1, 0)
TnS <- apply(proj.tmp1, 2, sum) / pop.size
plot(1:(days.proj+1), TnS, type="l", xlab="days", ylab="proportion susceptible")

# proportion susceptible
proj.tmp2 <- ifelse(proj.dat == "E", 1, 0)
TnE <- apply(proj.tmp2, 2, sum) / pop.size
plot(1:(days.proj+1), TnE, type="l", xlab="days", ylab="proportion exposed")

# proportion infected
proj.tmp3 <- ifelse(proj.dat == "I", 1, 0)
TnI <- apply(proj.tmp3, 2, sum) / pop.size
plot(1:(days.proj+1), TnI, type="l", xlab="days", ylab="proportion infected")

# proportion recovered
proj.tmp4 <- ifelse(proj.dat == "R", 1, 0)
TnR <- apply(proj.tmp4, 2, sum) / pop.size
plot(1:(days.proj+1), TnR, type="l", xlab="days", ylab="proportion recovered")
par(mfrow=c(1,1))
