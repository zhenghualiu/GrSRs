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

    dBb <- GRb1 * Bb #Bacterial growth
    
    dB_ = 
      GR1 * B_ - #Bacterial growth
      delta(Tem,deltaC) * (1-I(Tem,Ic,T1,sd1)) * B_ * P + #Lysis, lysogenization
      delta(Tem,deltaC) * (1-I(Tem,Ic,T1,sd1)) * 
      dpois(0,(AMn + P/N) * B_/N) * B_ * P #Escape infection
    
    dBi = GR2 * Bi + #Lysogeny growth
      delta(Tem,deltaC) * (1-I(Tem,Ic,T1,sd1)) * alpha * 
      dpois(n,(AMn + P/N) * B_/N) * B_ * P - #Lysogenization
      cauchy2(Tem,cauchy0,sc,T2,T3) * Bi #Induction
    
    dP = sum(delta(Tem,deltaC) * (1-I(Tem,Ic,T1,sd1)) *  (1 - alpha) * 
               dpois(n,(AMn + lag[P.id]/lag[N.id])  * lag[B_.id]/lag[N.id]) *
               beta(n,beta0) * lag[Bi.id] * lag[P.id],na.rm = T) + #Lysis
      cauchy2(Tem,cauchy0,sc,T2,T3) * sum( beta(n,beta0) * Bi ,na.rm = T) - #Induction
      delta(Tem,deltaC) * N * P #Adsorption
    
    dN = GR1 * B_ + #Bacterial growth
      sum(GR2 * Bi, na.rm = T ) + #Lysogeny growth
      dBb # Background bacteria
    
    return(c(B_=dB_,Bi=dBi,P=dP,Bb=dBb,N=dN))
  })
}

#2. Define initial values and parameters####
Mn = 100 # Maximal number of prophage in host
B_ = 4*10^6 #uninfected bacteria, cells/mL
Bb = 1*10^7 #background 
Bi = rep(0,100) #lysogeny, cells/mL
P = 10^5 #phages, phages/mL
N = B_ + sum(Bi) + Bb # population size

yinit <- c(B_ = B_, Bi = Bi, P = P, Bb = Bb, N = N)

parms <- list(
  initial = yinit, #Initial state
  tau = 1, # Latent time, h
  Mn = 100, # Maximal number of prophage in host
  AMn = 2, # Averagely minnimal number of phage in host
  alpha = 0.5, #Probability of lysogeny
  C = 2 * 10^8, # environment capacity, cells/mL
  Vm = 0.06, # maximal growth rate, h-1
  Vmb = 0.1, # maximal growth rate of background, h-1
  Mm = 0.1, # minimal mortality rate, h-1
  Tem = 273.15 + 40, # Temperature at time t, K
  T0 = 273.15 + 30, # Optimal growth temperature, K
  Tb = 273.15 + 40, # Optimal growth temperature of background, K
  sd0 = 15, # Niche breadth of temperature, K
  deltaC = 2.3*10^(-8), #Adsorption rate, mL/h
  Ic = 1, #Immunity rate
  T1 = 273.15 + 37, # Optimal temperature for immunity activity, K
  sd1 = 2, # Temperature breadth of immunity activity, K
  K0 = 100, # Immunity capacity, phages
  cauchy0 = 2 * 10^(-5), # Induction rate constant, mL/h
  T2 = 273.15 + 38, # Temprature switch constant for heat lysis, K
  T3 = 273.15 + 3, # Temprature switch constant for cold lysis, K
  sc = 4,
  beta0 = 8*10^2 # Burst size, particles
)

#3. Qualitative analysis####
library(ggplot2)
Ngen=1000 # Generations
#3.1 Temperature####
#Thermophilic, Mesophilic, Psychrophilc 
logDT <- seq(-1,2,0.1)
Vm <- 1/10^(logDT)
Thermophilic <- data.frame(Tem=45,Vm=Vm,T0=45,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Mesoophilic <- data.frame(Tem=30,Vm=Vm,T0=30,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Psychrophilc  <- data.frame(Tem=3,Vm=Vm,T0=30,Vmb=0.1,B_=4*10^6,Bb=8*10^6)
Psychrophilc2  <- data.frame(Tem=3,Vm=Vm,T0=30,Vmb=0.15,B_=4*10^6,Bb=8*10^6)

parms_set <- rbind.data.frame(Thermophilic,Mesoophilic)
parms_set <- rbind.data.frame(parms_set,Psychrophilc)
parms_set <- rbind.data.frame(parms_set,Psychrophilc2)


#modify_parms1$g
parms_set_g <- apply(parms_set,1,function(x){
#parms_set_g<- apply(parms_set,1,function(x){
  B_ = x[5] #uninfected bacteria, cells/mL
  Bi = rep(0,100) #lysogeny, cells/mL
  P = 10^5 #phages, phages/mL
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

#save(parms_set,file="parms_set.Rdata")
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
  ifelse(parms_set$Tem==45,"Thermophiles",
         parms_set$Type)
parms_set$Type <- 
  ifelse(parms_set$Tem==3,"Psychrotrophiles",
         parms_set$Type)


parms_set$Type <-
  factor(parms_set$Type,levels = rev(
    c("Psychrotrophiles","Mesophiles","Thermophiles")))
#Fig.4A####
parms_set %>% 
  filter(Vmb == 0.15) %>%
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
  labs(y=expression("Host specialization"),
       x=expression("Host growth rate (log10)" ))+
  guides(colour=
           guide_legend(title=""))+
  theme(legend.position = c(0.25,0.3))+
  scale_colour_manual(values =  brewer.pal(5,"Set1")[c(1,2,5 )],
                      labels = c(
                        "Thermophiles (OGT = 45℃, Env = 45℃)",
                        "Mesophiles (OGT = 30℃, Env = 30℃)",
                        "Psychrotrophiles (OGT = 30℃, Env = 3℃)"                   )
                                 ) +
  theme(legend.text = element_text(size=10),
        legend.background = element_blank())-> Fig.4A

#3.2 rates of sensitive/non-sensitive population density ####
logDT <- seq(-1,2,0.1)
Vm <- 1/10^(logDT)
Thermophilic <- data.frame(Tem=45,Vm=Vm,T0=45,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Mesoophilic <- data.frame(Tem=30,Vm=Vm,T0=30,Vmb=0.15,B_=4*10^6,Bb=8*10^6)
Psychrophilc  <- data.frame(Tem=3,Vm=Vm,T0=30,Vmb=0.1,B_=4*10^6,Bb=8*10^6)
Psychrophilc2.2  <- data.frame(Tem=3,Vm=Vm,T0=30,Vmb=0.15,B_=4*10^6,Bb=8*10^6)

tmp1 <- Thermophilic[,1:5]
tmp2 <- Mesoophilic[,1:5]
tmp3 <- Psychrophilc[,1:5]
tmp4 <- Psychrophilc2.2[,1:5]
Thermophilic2 <- NULL
Mesoophilic2 <- NULL
Psychrophilc2 <- NULL
Psychrophilc2.2.2 <- NULL
s_ns_ratio <- c(1/10^(seq(-1,5,0.2)))
for(i in s_ns_ratio)
{
  tmp1$Bb <- tmp1$B_*i
  tmp2$Bb <- tmp2$B_*i
  tmp3$Bb <- tmp3$B_*i
  tmp4$Bb <- tmp4$B_*i
  if(is.null(Thermophilic2)) {
    Thermophilic2 <- tmp1
    Mesoophilic2 <- tmp2
    Psychrophilc2 <- tmp3
    Psychrophilc2.2.2 <- tmp4
  }else{
    Thermophilic2 <- rbind.data.frame(Thermophilic2,tmp1)
    Mesoophilic2 <- rbind.data.frame(Mesoophilic2,tmp2)
    Psychrophilc2 <- rbind.data.frame(Psychrophilc2,tmp3)
    Psychrophilc2.2.2 <- rbind.data.frame(Psychrophilc2.2.2,tmp4)
    }
}


parms_set2 <- rbind.data.frame(
  Thermophilic2,Mesoophilic2,Psychrophilc2,
  Psychrophilc2.2.2
)

parms_set2_g<- apply(parms_set2,1,function(x){
  B_ = x[5] #uninfected bacteria, cells/mL
  Bi = rep(0,100) #lysogeny, cells/mL
  P = 10^5 #phages, phages/mL
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

parms_set2$g <- unlist(lapply(parms_set2_g,function(x){return(x$g)}))
save(parms_set2,file="parms_set2.Rdata")
load("parms_set2.Rdata")

parms_set2$group <- log10(parms_set2$B_/parms_set2$Bb)
parms_set2$group <- parms_set2$group
parms_set2$Type <- "Mesophiles"
parms_set2$Type <-
  ifelse(parms_set2$T0==45,"Thermophiles",parms_set2$Type)
parms_set2$Type <-
  ifelse(parms_set2$Tem==3,"Psychrotrophiles",parms_set2$Type)

parms_set2$Type <- 
  factor(parms_set2$Type,
         levels = c("Thermophiles","Mesophiles","Psychrotrophiles"))

load("parms_set2.Rdata")
#Fig.4B####
parms_set2 %>% 
#  filter(Type=="Psychrotrophilies") %>%
#  filter(Type=="Mesophilies") %>%
#  filter(Type=="Thermophiles") %>%
  filter(group >= 0.5 &
           Vmb == 0.15
           ) %>% 
  ggplot(aes(x=log10(1/Vm),1/g,
             colour=group,
             group=group))+
  geom_line(size=0.8)+mytheme+
  labs(y=expression("Host specialization"),
       x=expression("Host growth rate (log10)" ))+
#  theme(legend.position = c(0.2,0.7))+
  theme(legend.text = element_text(size=10),
        legend.background = element_blank(),
        legend.position = "bottom"
        )+
  facet_wrap(~Type,scales = "free_y")+
  guides(colour=guide_colorbar(title=expression(
    B^"-"/B^"+"~ratio)))+
  scale_colour_gradient2(high="#CC3333",low="#0066CC",#mid="grey",
                         midpoint = 3,
                         breaks=seq(1,5),
                         labels=c(expression("10"^"1"),
                                  expression("10"^"2"),
                                  expression("10"^"3"),
                                  expression("10"^"4"),
                                  expression("10"^"5")
                                  )) -> Fig.4B

cowplot::plot_grid(Fig.4A,Fig.4B,
                   labels = c('A', 'B'),nrow=2) -> Fig.4
ggsave("../04.Figure/02.Publish/Fig.4.tiff",
       plot=Fig.4,device = "tiff",compression="lzw",
       width = 6.6*1.1,height = 4.8*1.5*1.1,units = "in")


#4. Temperature-dependent function####
#4.1 lytic switch ####
Temperature <- seq(-10,50,0.1)
Cauchy <- sapply(Temperature,cauchy2,1,4,38,3)
df <- data.frame(Tem=Temperature,Cau=Cauchy)
df %>% ggplot(aes(x=Tem,y=Cau)) + 
  geom_line()+mytheme+
  labs(x= "Temperature (℃)", y="Induction rate") +
  geom_vline(xintercept = 4,linetype="dashed",color="red")+
  geom_vline(xintercept = 37,linetype="dashed",color="red")+
  geom_text(x=7,y=0.5,label="4 ℃")  +
  geom_text(x=33,y=0.5,label="37 ℃")  -> p_induction_switch
ggsave("../04.Figure/02.Publish/p_induction_switch.tiff",
       plot=p_induction_switch,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in")

#4.2 Temperarture-dependent lytic switch ####

#modify_parms1$g
parms_set_g <- apply(parms_set,1,function(x){
  #parms_set_g<- apply(parms_set,1,function(x){
  B_ = x[5] #uninfected bacteria, cells/mL
  Bi = rep(0,100) #lysogeny, cells/mL
  P = 10^5 #phages, phages/mL
  Bb = x[6] #background 
  N = B_ + sum(Bi) + Bb # population size
  yinit <- c(B_ = B_, Bi = Bi, P = P, Bb = Bb, N = N)
  parms$yinit <- yinit
  parms$Tem <- 273.15 + x[1]
  parms$T0 <- 273.15 + x[3]
  parms$Tb <- parms$Tem
  parms$Vm <- x[2]
  parms$Vmb <- x[4]
  #Setting Temperarture-independent lytic switch#
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
parms_set$g <- unlist(lapply(parms_set_g,
                             function(x){return(x$g)}))
parms_set3 <- parms_set

#save(parms_set3,file="parms_set3.Rdata")
load("parms_set3.Rdata")

parms_set3$Type <- "Mesophiles"
parms_set3$Type <- 
  ifelse(parms_set3$Tem==45,"Thermophiles",
         parms_set3$Type)
parms_set3$Type <- 
  ifelse(parms_set3$Tem==3,"Psychrotrophiles",
         parms_set3$Type)


parms_set3$Type <-
  factor(parms_set3$Type,levels = rev(
    c("Psychrotrophiles","Mesophiles","Thermophiles")))
parms_set3 %>% 
  filter(Vmb == 0.15) %>%
  ggplot(aes(x=log10(Vm),1/g,colour=Type))+
  geom_line(size=0.8)+mytheme+
  labs(x=expression(" Host growth rate (log10)"),y="Host specialization")+
  guides(colour=
           guide_legend(title=""))+
  theme(legend.position = c(0.38,0.45))+
  scale_colour_manual(values = brewer.pal(5,"Set1")[c(1,2,5 )],
                      labels = c(
                        "Thermophiles (OGT = 45℃, Env = 45℃)",
                        "Mesophiles (OGT = 30℃, Env = 30℃)",
                        "Psychrotrophiles (OGT = 30℃, Env = 3℃)"                   )
  ) +
  theme(legend.text = element_text(size=10),
        legend.background = element_blank()) -> Fig.4SA

ggsave("../04.Figure/02.Publish/p_induction_switch_no_Ly.tiff",
       plot=Fig.4SA,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in")
