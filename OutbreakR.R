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
# (University of South Florida); Taylor Callicrate (Chicago Zoological Society; Species Conservation Toolkit Initiative); Barry
# Brook (University of Tasmania)

rm(list=ls(all=TRUE))

# source functions
setwd("~/Documents/Papers/Disease/OutbreakR/analysis/functions")
source("age.disL.r") 
source("disParams.r") 
source("indChar.r")
source("OutbreakR.init.r")

#####################################################################################################
#####################################################################################################
# set initial parameters
pop.size <- 10000
#age.classes <- c("juv","sub","adu") # can also be a vector of ages
age.classes <- 0:17
#age.prop <- c(0.50, 0.30, 0.20) # can also be a stable st(age) distribution
age.prop <- c(0.098905556,0.08489517,0.07179324,0.05978244,0.04395387,0.031761248,0.02265974,0.017012779,0.012820156,0.008809558,0.006616935,0.0044243117,0.0036864416,0.0025847272,0.001665039,7.4535067E-4,5.5293506E-4,5.1192205E-4)
life.span <- 17 # in years
sex.classes <- c("F","M")
sex.ratio <- 1 # proportion female

disease.states <- c("P", "S", "E", "I",  "R")
init.dis.prop <- c(0.35, 0.5, 0.05, 0.1, 0.0) # proportion of population in each disease state; ensure sums to 1
sum(init.dis.prop)
library(crayon)
if (sum(init.dis.prop) != 1) {
  cat(red("** initial disease-state proportions do not sum to one **\n"))
}

# set age-specific likelihoods of disease state (e.g., all juveniles "P"; randomise rest)
age.disL.mat <- age.disL(age.classes, disease.states)
age.disL.mat[1,2:5] <- 0 # under-yearlings can only be P

# How many days to project?
days.proj <- 365 # number of days to project disease transition states into the future

# disease parameters
encounter.rate <- 5 # individuals per day
prob.never.S <- 0.001 # probability pre-susceptible never becomes susceptible (permanent immunity)
min.tim.b.ind.S <- 1 # minimum time before an individual becomes susceptible
max.tim.b.ind.S <- 2 # maximum time before an individual becomes susceptible
contact.rate <- 0.001 # proportion of the population an exposed/infected individual encounters/day
daily.trans.pr <- 0.00027 # daily transmission probability
out.dis.contact.rate <- 0.012 # outside source contact rate
out.dis.pr.trans <- 0.00027 # outside disease source probabilty of transmission
min.lat.inc.E <- 150 # minimum (days) incubation (latent) period of infection for exposed individuals
max.lat.inc.E <- 7*365 # maximum (days) incubation (latent) period of infection for exposed individuals
prop.I.remain <- 0.425 # proportion of infected individuals remaining so indefinitely
min.tim.rem.I <- 60 # minimum number of days an animal can remain infectious
max.tim.rem.I <- 17*365 # maximum number of days an animal can remain infectious
pr.R <- 0.001 # probability of recovery and subsequent resistance after the minimum time being infectious 
pr.ret.S <- 0.8 # probability of returning to the susceptible state after recovery (assuming those not returning to susceptible are pre-susceptible)
pr.R.perm.imm <- 0.01 # proportion of recovered individuals acquiring permanent immunity
min.tim.res <- 30 # minimum number of days resistant
max.tim.res <- 5*365 # maximum number of days resistant

# accumulate disease parameters
displist <- disParams(prob.never.S,min.tim.b.ind.S,max.tim.b.ind.S,contact.rate,daily.trans.pr,out.dis.contact.rate,out.dis.pr.trans,
                      min.lat.inc.E,max.lat.inc.E,prop.I.remain,min.tim.rem.I,max.tim.rem.I,pr.R,pr.ret.S,pr.R.perm.imm,min.tim.res,max.tim.res,life.span)

# randomise initial individual characteristics
init.dat <- indChar(a=age.prop, as=age.classes, s=sex.ratio, d=init.dis.prop, ds=disease.states, disL=age.disL.mat, n=pop.size)

# project first outbreak run
OutbreakR.list <- OutbreakR.init(initdat=init.dat, disparams=displist, days=days.proj, n=pop.size)


#################################################
# simple demography set by age/sex/disease state
surv.vec <- c(0.3,0.5,0.8,0.95,rep(0.97,(18-4))) # age-specific survival vector
pbreed.vec <- c(0,0,0,0.2,0.5,0.8,rep(0.95,(18.6))) # age-specific proportion breeding vector
fert.vec <- c(0,0,0,rep(1,(18-3))) # female fertilty (daughters

surv.dis.mod <- c(1,1,1,0.8,1) # disease state modifier of survival vector (P, S, E, I, R)
fert.dis.mod <- c(1,1,1,0.9,1) # disease state modifier of fertility vector (P, S, E, I, R)
dis.mod.dat <- data.frame(surv.dis.mod, fert.dis.mod)
colnames(dis.mod.dat) <- c("survmod", "fertmod")
rownames(dis.mod.dat) <- disease.states

# set new matrix of individuals taken from initial OutbreakR projection
for (j in 1:length(age.classes)) {
  age.tmp <- subset(OutbreakR.list[[1]], age==age.classes[j])
  
  for (i in 1:dim(OutbreakR.list[[1]])[1]) {
    age.tmp[i,]$end.dis
    
  } # end i loop
} # end j loop
