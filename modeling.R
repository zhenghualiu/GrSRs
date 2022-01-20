#0 Environment setting####
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#0.1 Packages####
library(PBSddesolve)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
#0.2 Source code####
Rfiles <- dir("./modeling_funcs/")
res <- lapply(Rfiles,function(x){
  tmp <- paste("./modeling_funcs/",x,sep="")
  source(tmp)
})

#1. model####
mod <- function(Time, State, Pars) {
  with(as.list(c(Pars)), {
    if (Time < tau)
      lag <- initial
    else
      lag <- pastvalue(Time-tau)
    #Col index
    B_.id <- grep("B_",names(initial))
    Bi.id <- grep("Bi",names(initial))
    P.id <- grep("P",names(initial))
    N.id <- grep("N",names(initial))
    Bb.id <- grep("Bb",names(initial))
    # Number of prophage in host
    n <- seq(1,Mn)
    #State
    B_ <- State[B_.id]
    Bi <- State[Bi.id]
    P <- State[P.id]
    N <- State[N.id]
    Bb <- State[Bb.id]
    #growth rate
    GR = gr(N,Vm,C)
    GR1 = GR * theta1(Tem,T0,sd0)
    GR2 = GR * theta1(Tem,T0,sd0)
    GRb = gr(N,Vmb,C)
    GRb1 = GRb * theta1(Tem,Tb,sd0)

    dBb <- GRb1 * Bb #B+ growth
    
    dB_ = 
      GR1 * B_ - #B- growth
      deltaC  * B_ * P + #Lysis, lysogenization
      deltaC *  
#      dpois(0,AMn + P/N) * B_ * P #Escape infection
      dpois(0,(AMn + P/N)*B_/N) * B_ * P
    
    dBi = GR2 * Bi + #Lysogens Bi growth
      deltaC * alpha * 
      dpois(n,(AMn + P/N)*B_/N) * B_ * P - #Lysogenization
      cauchy2(Tem,cauchy0,sc,T2,T3) * Bi #Induction
    
    dP = sum(deltaC *  (1 - alpha) * 
               dpois(n,(AMn + P/N)*B_/N) *
               beta0 * Bi * P,na.rm = T) + #Lysis
      cauchy2(Tem,cauchy0,sc,T2,T3) * sum( beta0 * lag[Bi.id] ,na.rm = T) - #Induction
   #   sum(deltaC * n * dpois(n, AMn + P/N) * B_ * P) #Adsorption
      deltaC * N * P
    
    dN = dB_ + #B-
      sum(dBi, na.rm = T ) + #Lysogens
      dBb # B+
    return(c(B_=dB_,Bi=dBi,P=dP,Bb=dBb,N=dN))
  })
}

#2. Define initial values and parameters####
Mn = 10 # Maximal number of prophage in host
B_ = 4*10^6 #susceptible cells, cells/mL
Bb = 8*10^6 #nonsusceptible cells, cells/mL
Bi = rep(0,Mn) #lysogens, cells/mL
P = 10^5 #virus, particles/mL
N = B_ + sum(Bi) + Bb # total cell population density, cells/mL

yinit <- c(B_ = B_, Bi = Bi, P = P, Bb = Bb, N = N)

parms <- list(
  initial = yinit, #Initial state
  tau = 1, # Latent time, h
  AMn = 2, # Averagely minnimal number of viruses in host
  alpha = 0.5, #Probability of being lysogens
  C = 2 * 10^8, # environment capacity, cells/mL
  Vm = 0.5, # maximal growth rate of susceptible cells, hr-1
  Vmb = 0.5, # maximal growth rate of nonsusceptible cells, hr-1
  Tem = 273.15 + 40, # Temperature at time t, K
  T0 = 273.15 + 30, # Optimal growth temperature of susceptible cells, K
  Tb = 273.15 + 40, # Optimal growth temperature of nonsusceptible cells, K
  sd0 = 15, # Niche breadth of temperature, K
  deltaC = 4*10^(-8), #Adsorption rate, mL/h
  cauchy0 = 5 * 10^(-5), # Induction rate constant, mL/h
  T2 = 273.15 + 37, # Temprature switch constant for heat lysis, K
  T3 = 273.15 + 4, # Temprature switch constant for cold lysis, K
  sc = 4, # constant parameter
  beta0 = 5*10^2 # Burst size, particles
)

Ngen=1000 # Number of time steps

#3. Qualitative analysis####
#Thermophilic, Mesophilic, Psychrophilc 
logGr <- seq(-2,1.2,0.1)
Vm <- 10^logGr
Thermophilic <- data.frame(Tem=37,Vm=Vm,T0=45,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Mesoophilic <- data.frame(Tem=30,Vm=Vm,T0=30,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Psychrophilc  <- data.frame(Tem=4,Vm=Vm,T0=20,Vmb=0.15,B_=4*10^6,Bb=8*10^6)

parms_set <- rbind.data.frame(Thermophilic,Mesoophilic)
parms_set <- rbind.data.frame(parms_set,Psychrophilc)

parms_set_g <- apply(parms_set[,1:6],1,function(x){
  B_ = x[5] #uninfected bacteria, cells/mL
  Bi = rep(0,10) #lysogens, cells/mL
  P = 10^5 #phage particles, particles/mL
  Bb = x[6] #background 
  N = B_ + sum(Bi) + Bb # population size
  yinit <- c(B_ = B_, Bi = Bi, P = P, Bb = Bb, N = N)
  parms$yinit <- yinit
  parms$Tem <- 273.15 + x[1]
  parms$T0 <- 273.15 + x[3]
  parms$Tb <- parms$Tem
  parms$Vm <- x[2]
  parms$Vmb <- x[4]
  out.tmp <- dde(y=yinit,times=seq(0,Ngen,1),func=mod,parms=parms)
  Bi.id <- grep("Bi",colnames(out.tmp))
  bi <- out.tmp[-1,Bi.id]
  bi2 <- t(apply(bi,1,function(x){x/sum(x)}))
  bi3 <- t(apply(bi2,1,function(x){x*seq(1,Mn)}))
  g_ <- rowSums(bi3)
  g <- sum(as.numeric(out.tmp[Ngen,Bi.id]/sum(out.tmp[Ngen,Bi.id])*seq(1,Mn)))
  return(list(g=g,g_=g_))
})
parms_set$g <- unlist(lapply(parms_set_g,
                             function(x){return(x$g)}))

save(parms_set,file="parms_set.Rdata")
load("parms_set.Rdata")
mytheme <- theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black",size=12),
        strip.background = element_blank(),
        axis.title = element_text(size=15,color="black"),
        strip.text = element_text(color="black",size=12)
  )
parms_set$Type <- "Mesophiles"
parms_set$Type <- 
  ifelse(parms_set$T0==45,"Thermophiles",
         parms_set$Type)
parms_set$Type <- 
  ifelse(parms_set$T0==20,"Psychrophiles",
         parms_set$Type)


parms_set$Type <-
  factor(parms_set$Type,levels = rev(
    c("Psychrophiles","Mesophiles","Thermophiles")))

parms_set %>% 
  ggplot(aes(x=log10(Vm),1/g,colour=Type))+
  geom_vline(xintercept = log10(0.15),
             linetype="dashed",
             colour="grey",
             size=1
             )+
#  geom_vline(xintercept = log10(1/(0.15/stress1(3,30,15))),
#             linetype="dashed",
#             colour="blue",
#             size=1
#  )+
  geom_line(size=0.8)+mytheme+
  labs(y="Virus specificity",
       x=expression(Log[10]~"host growth rate ("*hr^-1*")"))+
  guides(colour=
           guide_legend(title=""))+
  theme(legend.position = c(0.2,0.3)
        )+
  scale_colour_manual(values =  brewer.pal(5,"Set1")[c(1,2,5 )],
                      labels = c(
                        "Thermophiles (OGT = 45°C, Env = 37°C)",
                        "Mesophiles (OGT = 30°C, Env = 30°C)",
                        "Psychrophiles (OGT = 20°C Env = 4°C)" )
                                 ) +
  theme(legend.text = element_text(size=10),
        legend.background = element_rect(colour="white",color="black")
 )-> Fig.4A

#4. Temperature-dependent function####
#4.1 lytic switch ####
Temperature <- seq(-10,50,0.1)
Cauchy <- sapply(Temperature,cauchy2,1,4,37,4)
df <- data.frame(Tem=Temperature,Cau=Cauchy)
df %>% ggplot(aes(x=Tem,y=Cau)) + 
  geom_line()+mytheme+
  labs(x= "Temperature (°C)", y="Induction rate") +
  geom_vline(xintercept = 4,linetype="dashed",color="red")+
  geom_vline(xintercept = 37,linetype="dashed",color="red")+
  geom_text(x=7,y=0.5,label="4 °C")  +
  geom_text(x=33,y=0.5,label="37 °C")  -> p_induction_switch
ggsave("p_induction_switch.tiff",
       plot=p_induction_switch,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in")

#4.2 Without temperarture-dependent lytic switch ####
parms_set3 <- parms_set[,1:6]
#modify_parms1$g
parms_set_g2 <- apply(parms_set3,1,function(x){
  B_ = x[5] #uninfected bacteria, cells/mL
  Bi = rep(0,10) #lysogens, cells/mL
  P = 10^5 #phage particles, particles/mL
  Bb = x[6] #background 
  N = B_ + sum(Bi) + Bb # population size
  yinit <- c(B_ = B_, Bi = Bi, P = P, Bb = Bb, N = N)
  parms$yinit <- yinit
  parms$Tem <- 273.15 + x[1]
  parms$T0 <- 273.15 + x[3]
  parms$Tb <- parms$Tem
  parms$Vm <- x[2]
  parms$Vmb <- x[4]
  # Inactivating temperarture-independent lytic switch#
  parms$T2 <- 999
  parms$T3 <- -999
  #----------------------#
  out.tmp <- dde(y=yinit,times=seq(0,Ngen,1),func=mod,parms=parms)
  Bi.id <- grep("Bi",colnames(out.tmp))
  bi <- out.tmp[-1,Bi.id]
  bi2 <- t(apply(bi,1,function(x){x/sum(x)}))
  bi3 <- t(apply(bi2,1,function(x){x*seq(1,Mn)}))
  g_ <- rowSums(bi3)
  g <- sum(as.numeric(out.tmp[Ngen,Bi.id]/sum(out.tmp[Ngen,Bi.id])*seq(1,Mn)))
  return(list(g=g,g_=g_))
})
parms_set3$g <- unlist(lapply(parms_set_g2,
                             function(x){return(x$g)}))

save(parms_set3,file="parms_set3.Rdata")
load("parms_set3.Rdata")

parms_set3$Type <- "Mesophiles"
parms_set3$Type <- 
  ifelse(parms_set3$T0==45,"Thermophiles",
         parms_set3$Type)
parms_set3$Type <- 
  ifelse(parms_set3$T0==20,"Psychrophiles",
         parms_set3$Type)


parms_set3$Type <-
  factor(parms_set3$Type,levels = rev(
    c("Psychrophiles","Mesophiles","Thermophiles")))
parms_set3$LS <- "N"
parms_set$LS <- "Y"

rbind.data.frame(parms_set,
                 parms_set3) -> Fig.S4.plot
Fig.S4.plot$LS <- 
  factor(Fig.S4.plot$LS,levels = c("Y","N"))
Fig.S4.plot %>% 
  ggplot(aes(x=log10(Vm),1/g,colour=Type,linetype=LS))+
  geom_vline(xintercept = log10(0.15),
     #        linetype="dashed",
             colour="grey",
             size=0.8
  )+
  geom_line(size=0.8,alpha=8)+mytheme+
  labs(y="Specialization",
       x=expression(Log[10]~"host growth rate ("*hr^-1*")"))+
  guides(colour=F,
         linetype=F)+
  theme(legend.position = "bottom")+
  scale_colour_manual(values = brewer.pal(5,"Set1")[c(1,2,5 )],
                      labels = c(
                        "Thermophiles (OGT = 45°C, Env = 37°C)",
                        "Mesophiles (OGT = 30°C, Env = 30°C)",
                        "Psychrophiles (OGT = 20°C, Env = 4°C)"     )
  ) +
  theme(legend.text = element_text(size=10),
        legend.background = element_blank()
        #panel.grid.major  = element_line(size=0.3,colour="grey")
        )+
  facet_wrap(~Type)-> Fig.4SA


ggsave("p_induction_switch_Ly.tiff",
       plot=Fig.4SA,device = "tiff",compression="lzw",
       width = 9,height = 4,units = "in")

ggsave("p_induction_switch_Ly.pdf",
       plot=Fig.4SA,device = "pdf",
       width = 9,height = 4,scale = 1.1)
