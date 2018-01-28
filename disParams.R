# nS = probability pre-susceptible never becomes susceptible (permanent immunity)
# minbS = minimum time before an individual becomes susceptible
# maxbS = maximum time before an individual becomes susceptible
# cr = proportion of the population an exposed/infected individual encounters/day
# trp = daily transmission probability
# ocr = outside source contact rate
# otrp = outside disease source probabilty of transmission
# minlat = minimum (days) incubation (latent) period of infection for exposed individuals
# maxlat = maximum (days) incubation (latent) period of infection for exposed individuals
# Irem = proportion of infected individuals remaining so indefinitely
# minIrem = minimum number of days an animal can remain infectious
# maxIrem = maximum number of days an animal can remain infectious
# pR = probability of recovery and subsequent resistance after the minimum time being infectious 
# prS = probability of returning to the susceptible state after recovery (assuming those not returning to susceptible are pre-susceptible)
# pRp = proportion of recovered individuals acquiring permanent immunity
# minR = 30 # minimum number of days resistant
# maxR = 5*365 # maximum number of days resistant
# lspan = species lifespan in years

disParams <- function(nS, minbS, maxbS, cr, trp, ocr, otrp, minlat, maxlat, Irem, minIrem, maxIrem, pR, prS, pRp, minR, maxR, lspan) {
  disParams.vec <- c(nS, minbS, maxbS, cr, trp, ocr, otrp, minlat, maxlat, Irem, minIrem, maxIrem, pR, prS, pRp, minR, maxR, lspan)
  return(disParams.vec)
} # end disParams function
