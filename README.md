# OutbreakR
Code development for hybrid Outbreak R package

Emerging Case: OutbreakR: an R library for epidemiological projections of wildlife and human diseases

Project: 
Program Outbreak (www.vortex10.org/Outbreak.aspx) is an individual-based model software package (open access) that simulates disease dynamics using transitions among susceptible, exposed, infectious and recovered individuals. Outbreak provides several options for modes of transmission (random contact within populations, spatially based transmission, contact with environmental sources of disease) and provides options for management through vaccination or culling (written by Robert Lacy, JP Pollak, PS Miller, L Hungerford, and P Bright).

While Outbreak is ideal for modelling the dynamics of diseases in small populations, it could be limited for applications involving larger (>> 1000s) populations because of the detailed, individual-to-individual contacts invoked in the individual-based model at each transition interval. We therefore propose to write an R Package (cran.r-project.org) library that includes a large component of the main Outbreak methodology, but calculates disease-state transitions based on a weighted population of infected and susceptible individuals at each time step.

The main idea is to hybridise the individual- and population-based algorithms such that the disease dynamics of larger populations can be modelled. We propose to start with a detailed matrix of individuals in a population (rows), and their characteristics (age, sex, condition, etc.) and initial disease state (pre-susceptible, susceptible, exposed, infected, recovered) at the beginning of the projection. For each interval, the disease state of each individual is assessed as the product of the probabilities of contact, transmission, etc., weighted by the number of individuals in the corresponding disease state (e.g., susceptibles can transition to ‘exposed’ as a product of the contact probability, the transmission probability, and the number of infected individuals in the entire population).

At the end of the desired number of disease-transition intervals (e.g., daily for 1 year), the frequency table of individuals in each disease category can then be parsed to other applications in R (e.g., matrix population models) or for external applications (e.g., RAMAS/Metapop). Over time, the addition of spatial elements, environmental correlates, and other complexities can be added to the base functions in the library.

Meetings & Activities: 
RCN Synthesis Meeting, White Oak Conservation, Yulee, Florida, USA (15-19 January 2018)

Participants & Partners: 
Corey Bradshaw (Flinders University), 
Robert Lacy (Chicago Zoological Society), 
Phil Miller (Species Survival Commission: Species Planning Specialist Group), 
JP Pollak (University of Cornell), 
Ana Davidson (Colorado State University), 
Nadia Ali (University of Chicago), 
Loren Cassin Sackett (University of South Florida), 
Taylor Callicrate (Chicago Zoological Society; Species Conservation Toolkit Initiative), 
Nadia Ali (University of Chicago), 
Barry Brook (University of Tasmania)


Corey J. A. Bradshaw, 
Flinders University, 
January 2018
