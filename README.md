# GrSRs

To explore how temperature and population dynamic features modulate the GrSR, we set up a mathematical model of host-viral interactions. 
Briefly, we assumed that each step of the infection cycle is temperature dependent, all virus species are ecologically equivalent and that the size of the virus set for a species follows a Poisson distribution. 
The temperature thresholds for heat and cold inductions were set as 37℃ and 4℃, respectively. 

The 'modeling-v3.R' file includes main functions for modeling and result visualization.  

The 'modeling_funcs' directory contains key functions for the simulations of host traits with given parameters.  
beta.R        Burst size, particles  
delta.R       Induction rate  
growth.R      Growth rate  
growth_stress.R      Temperature stress on growth  
immunity.R        Immunity rate  
Induction2.R        Induction rate of lysogeny across temperature  
theta.R       Temperature niche of host cells.  
