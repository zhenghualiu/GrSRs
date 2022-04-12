#rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(reshape2)
library(reshape)
library(ggplot2)
library(ggpubr)
library(stringr)
library(plyr)
library(rlist)
library(RColorBrewer)
library(vegan)
library(GGally)
library(lme4)
source("../../../1.script/R/reshaping_object/VecListBind.R")
source("../../../1.script/R/longDFfuncs/cal_lm_p.R")
source("../../../1.script/R/longDFfuncs/cal_cor.R")
source("../../../1.script/R/longDFfuncs/compare_mean.R")
source("../../../1.script/R/reshaping_object/windows_split.R")
source("./mantel_genomic_trait_net.R")

load("../../02.Data/02.Rdata/0.host_metadata.Rdata")
load("../../02.Data/02.Rdata/0.virus_metadata.Rdata")
load("../../02.Data/02.Rdata/0.Vpc_funcs.Rdata")
load("../../02.Data/02.Rdata/1.host_pc_M_muti.Rdata")
load("../../02.Data/02.Rdata/1.host_vc_M_muti.Rdata")
load("../../02.Data/02.Rdata/2.dIndex.Rdata")
load("../../02.Data/02.Rdata/2.fast_res.Rdata")
load("../../02.Data/02.Rdata/2.iso_mod_members_vc_R2.Rdata")

load("../../02.Data/02.Rdata/0.vpc_file.Rdata")

iso_host_metadata$DoublingTime <-
  iso_host_metadata$DoublingTime/24
#0.ReName####
#0.1.host####
{
    {
  Old_VarNames_continuous_iso <- c(
    #Genome Traits
    "Assembly.Accession","Size_Mb","Genes","GC","num_5S","num_16S","num_23S",
    "num_sigma_factor","OGT", "CUBHE","ConsistencyHE","CPB","DoublingTime",
    "num_prophage","num_VC","num_VPC","num_VP","Num Temperate",
    "Num Virulent",
    #Phylogeny
    #Environment
    "latitude","longitude"
      )
    }#iso
    {
      Old_VarNames_continuous_bin <- c(
        #Genome Traits
        "bin_id","Size","num_CDS","GC","num_5S","num_16S","num_23S",
        "num_sigma_factor","OGT", "CUBHE","ConsistencyHE","CPB","DoublingTime",
        "num_prophage","num_VC","num_VPC","num_VP", "Num Temperate",
        "Num Virulent",
        #Phylogeny
        #Environment
        "latitude","longitude"
      )
    }#bin
    {
    New_VarNames_continuous <- c(
      #Genome Traits
      "GenomeID","Size","CDS","GC","Num 5S","Num 16S","Num 23S",
      "Num sigma factor","OGT","CUBHE","ConsistencyHE","CPB","DoublingTime",
      "Num prophage","Num VC","Num VPC","Num VP", "Num Temperate",
      "Num Virulent",
      #Phylogeny
      #Environment
      "Latitude","Longitude")
      }#New variable names
  iso_host_metadata_cont <-
    iso_host_metadata[,Old_VarNames_continuous_iso]
  bin_host_metadata_cont <-
    bin_host_metadata[,Old_VarNames_continuous_bin]
  colnames(iso_host_metadata_cont) <- New_VarNames_continuous
  colnames(bin_host_metadata_cont) <- New_VarNames_continuous
  iso_host_metadata_cont$Type="Iso"
  bin_host_metadata_cont$Type="MAG"
}#continuous variable
{
  all_vars <- colnames(iso_host_metadata)
    {
      Old_VarNames_discrete_iso <- c(
        "Assembly.Accession",
        #Phylogeny
        "kingdom","phylum","class","order","family","genus","species",
        #Environment
        "empo_3.5"
      )
      #Phenotype
      Old_VarNames_phenotype_iso <- 
        c("Assembly.Accession","cell_shape","motility","sporulation",
      "temperature_range","salinity","oxygen_requirement","gram_stain")
      #AntiCRISPR
      Old_VarNames_AntiCrispr <- 
        c("Assembly.Accession",all_vars[grepl("AntiCRISPR",all_vars)])
    }#iso
    {
      Old_VarNames_discrete_bin <- c(
        "bin_id",
        #Phylogeny
        "domain","phylum","class","order","family","genus","species",
        #Environment
        "empo_3.5"
      )
    }#bin
    {
      New_VarNames_discrete <- c(
        "GenomeID",
        #Phylogeny
        "Kingdom","Phylum","Class","Order","Family","Genus","Species",
        #Environment
        "Habitat")
      New_VarNames_phenotype <- 
        c("GenomeID","Cell shape","Motility","Sporulation",
          "Temperature range","Salinity","Oxygen requirement","Gram strain")
      New_VarNames_AntiCrispr <- 
        c("GenomeID",all_vars[grepl("AntiCRISPR",all_vars)])
    }#New variable names
  iso_host_metadata_disc <-
    iso_host_metadata[,Old_VarNames_discrete_iso]
  bin_host_metadata_disc <-
    bin_host_metadata[,Old_VarNames_discrete_bin]
  iso_host_metadata_phen <- 
    iso_host_metadata[,Old_VarNames_phenotype_iso]
  iso_host_metadata_antiCrispr <- 
    iso_host_metadata[,Old_VarNames_AntiCrispr]
  colnames(iso_host_metadata_disc) <- New_VarNames_discrete
  colnames(bin_host_metadata_disc) <- New_VarNames_discrete
  colnames(iso_host_metadata_phen) <- New_VarNames_phenotype
  colnames(iso_host_metadata_antiCrispr) <- New_VarNames_AntiCrispr
  iso_host_metadata_antiCrispr
  iso_host_metadata_disc$Type="Iso"
  bin_host_metadata_disc$Type="MAG"
}#Discrete variable
#0.2.virus####
{
  {
    Old_VarNames_continuous_iso_vc <- 
      c("VC", "Mean.Size.VC","Mean.GC.VC")
  }#iso
  {
    Old_VarNames_continuous_bin_vc <- 
       c("VC", "Mean.Size.VC","Mean.GC.VC")
  }#bin
  New_VarNames_continuous_vc <- c("VC", "Mean.Size","Mean.GC")
  iso_vc_metadata_cont <-
    iso_virus_meta2[,Old_VarNames_continuous_iso_vc]
  bin_vc_metadata_cont <-
    bin_virus_meta2[,Old_VarNames_continuous_bin_vc]
  colnames(iso_vc_metadata_cont) <- New_VarNames_continuous_vc
  colnames(bin_vc_metadata_cont) <- New_VarNames_continuous_vc
  iso_vc_metadata_cont$Type="Iso"
  bin_vc_metadata_cont$Type="MAG"
  iso_vc_metadata_cont <- 
    iso_vc_metadata_cont[!duplicated(iso_vc_metadata_cont),]
  bin_vc_metadata_cont <- 
    bin_vc_metadata_cont[!duplicated(bin_vc_metadata_cont),]
}#continuous variable
{
  {
    Old_VarNames_discrete_iso_vc <- 
      c("VC", "Order","Family","Genus","LifeStyle" )
  }#iso
  {
    Old_VarNames_discrete_bin_vc <- 
      c("VC", "Order","Family","Genus","LifeStyle" )
  }#bin
  iso_vc_metadata_disc <-
    iso_virus_meta2[,Old_VarNames_discrete_iso_vc]
  bin_vc_metadata_disc <-
    bin_virus_meta2[,Old_VarNames_discrete_bin_vc]
  iso_vc_metadata_disc$Type="Iso"
  bin_vc_metadata_disc$Type="MAG"
  iso_vc_metadata_disc <- 
    iso_vc_metadata_disc[!duplicated(iso_vc_metadata_disc),]
  bin_vc_metadata_disc <- 
    bin_vc_metadata_disc[!duplicated(bin_vc_metadata_disc),]
}#Discrete variable

#1.Merge Data####
#1.1.dIndex for host####
{
   {
    iso_host_pc_d_host <-
       dIndex$host_pc_d_host[grep("iso",names(dIndex$host_pc_d_host))]
    iso_host_vc_d_host <-
      dIndex$host_vc_d_host[grep("iso",names(dIndex$host_vc_d_host))]
    iso_host_pc_d_pc <-
      dIndex$host_pc_d_host[grep("iso",names(dIndex$host_pc_d_pc))]
    iso_host_vc_d_vc <-
      dIndex$host_vc_d_host[grep("iso",names(dIndex$host_vc_d_vc))]
    #list 2 dataframe
    iso_d_host <- VecListBind(c(iso_host_pc_d_host,iso_host_vc_d_host))
    iso_d_pc <- VecListBind(iso_host_pc_d_pc)
    iso_d_vc <- VecListBind(iso_host_vc_d_vc)
    #rename col
    colnames(iso_d_host) <- str_replace(colnames(iso_d_host),"iso_","")
    colnames(iso_d_pc) <- str_replace(colnames(iso_d_pc),"iso_","")
    colnames(iso_d_vc) <- str_replace(colnames(iso_d_vc),"iso_","")
    #
    iso_host_stat_cont <- merge.data.frame(
      iso_host_metadata_cont,iso_d_host,by.x="GenomeID",by.y = 0)
    iso_host_stat_disc <- merge.data.frame(
      iso_host_metadata_disc,iso_d_host,by.x="GenomeID",by.y = 0)
    iso_host_stat_phen <- merge.data.frame(
      iso_host_metadata_phen,iso_d_host,by.x = "GenomeID",by.y = 0)
    iso_host_stat_all <- merge.data.frame(
      iso_host_metadata,iso_d_host,by.x = "Assembly.Accession",by.y = 0)
  }#iso 
   {
    bin_host_pc_d_host <-
      dIndex$host_pc_d_host[grep("bin",names(dIndex$host_pc_d_host))]
    bin_host_vc_d_host <-
      dIndex$host_vc_d_host[grep("bin",names(dIndex$host_vc_d_host))]
    bin_host_pc_d_pc <-
      dIndex$host_pc_d_host[grep("bin",names(dIndex$host_pc_d_pc))]
    bin_host_vc_d_vc <-
      dIndex$host_vc_d_host[grep("bin",names(dIndex$host_vc_d_vc))]
    #list 2 dataframe
    bin_d_host <- VecListBind(c(bin_host_pc_d_host,bin_host_vc_d_host))
    bin_d_pc <- VecListBind(bin_host_pc_d_pc)
    bin_d_vc <- VecListBind(bin_host_vc_d_vc)
    #rename col
    colnames(bin_d_host) <- str_replace(colnames(bin_d_host),"bin_","")
    colnames(bin_d_pc) <- str_replace(colnames(bin_d_pc),"bin_","")
    colnames(bin_d_vc) <- str_replace(colnames(bin_d_vc),"bin_","")
    #
    bin_host_stat_cont <- merge.data.frame(
      bin_host_metadata_cont,bin_d_host,by.x="GenomeID",by.y = 0)
    bin_host_stat_disc <- merge.data.frame(
      bin_host_metadata_disc,bin_d_host,by.x="GenomeID",by.y = 0)
    
  }#bin
}
#1.2 fast parameters####
{
{
    fast_higer <- lapply(fast_res,function(x){
    x$higher$ID <- rownames(x$higher)
    return(x$higher)})
  fast_lower <- lapply(fast_res,function(x){
    x$lower$ID <- rownames(x$lower)
    return(x$lower)})
  fast_host_merge <- ldply(fast_lower)
  colnames(fast_host_merge)[1] <- "Net"
  fast_host_merge$Net <- str_replace(fast_host_merge$Net,"iso_","")
  fast_host_merge$Net <- str_replace(fast_host_merge$Net,"bin_","")
}#host
  fast_higer_merge <- ldply(fast_higer)
}

Target_vars <- c(colnames(iso_d_host))

#2.Reshaping data to long data####
#2.1.host####
{
  {
host_stat_cont <- rbind.data.frame(iso_host_stat_cont,bin_host_stat_cont)
host_stat_disc <- rbind.data.frame(iso_host_stat_disc,bin_host_stat_disc)
host_stat_cont$Latitude <- abs(host_stat_cont$Latitude)
host_stat_cont$Longitude <- abs(host_stat_cont$Longitude)
{
host_stat_cont.melt <- 
  melt.data.frame(host_stat_cont,measure.vars = New_VarNames_continuous[-1],
                  variable_name = "Cont_Vars")
host_stat_cont.melt$value <- as.numeric(host_stat_cont.melt$value)
host_stat_cont.melt %>%
  filter(!(Cont_Vars == "DoublingTime" & value > 50)) -> host_stat_cont.melt
colnames(host_stat_cont.melt)[
  which(colnames(host_stat_cont.melt) == "value")] <- "Cont_Var_value"
host_stat_cont.melt <- 
  melt.data.frame(host_stat_cont.melt,
                  measure.vars = Target_vars,
                  variable_name = "Target_Vars")
host_stat_cont.melt$value <- as.numeric(host_stat_cont.melt$value)
colnames(host_stat_cont.melt)[
  which(colnames(host_stat_cont.melt) == "value")] <- "Target_Var_value"
}#d'
{
  host_metadata_cont <- 
    rbind(iso_host_metadata_cont,
          bin_host_metadata_cont)
  fast_host_metadata <- 
    merge.data.frame(fast_host_merge,host_metadata_cont,by.x = "ID",
                     by.y="GenomeID")
  fast_host_metadata.melt <- 
    melt.data.frame(fast_host_metadata,
                    measure.vars = New_VarNames_continuous[-1],
                    variable_name = "Cont_Vars")
  colnames(fast_host_metadata.melt)[
    colnames(fast_host_metadata.melt) == "value"] <- "Cont_Var_value"
  fast_host_metadata.melt$Cont_Var_value <- 
    as.numeric(fast_host_metadata.melt$Cont_Var_value)
}#fast parameters
} #host_stat_cont
  {
  iso_host_stat_phen.melt <- 
    melt.data.frame(iso_host_stat_phen,
                    measure.vars = New_VarNames_phenotype[-1],
                    variable_name = "Phenotype")
  colnames(iso_host_stat_phen.melt)[
    which(colnames(iso_host_stat_phen.melt) == "value")] <- "Phenotype_value"
  iso_host_stat_phen.melt <- 
    melt.data.frame(
      iso_host_stat_phen.melt,
      measure.vars = Target_vars,
      variable_name = "Target_Vars")
  colnames(iso_host_stat_phen.melt)[
    which(colnames(iso_host_stat_phen.melt) == "value")] <- "Target_Vars_value"
  iso_host_stat_phen.melt %>% filter(!is.na(Phenotype_value)) %>%
    group_by(Phenotype,Target_Vars) %>%
    do(data.frame(GenomeID=.$GenomeID,
                  Phenotype=.$Phenotype,
                  Phenotype_value=.$Phenotype_value,
                  Target_Vars=.$Target_Vars,
                  Target_Vars_value=.$Target_Vars_value,
                  P=compare_mean(.$Target_Vars_value,.$Phenotype_value)
                  )) -> iso_host_stat_phen.melt
  iso_host_stat_phen.melt$sig <- 
    ifelse(iso_host_stat_phen.melt$P < 0.05,"p<0.05","ns")
}#phenotype_d
}
#2.2.viral cluster####
{
  vc_metadata_cont <- rbind.data.frame(iso_vc_metadata_cont,bin_vc_metadata_cont)
  vc_metadata_disc <- rbind.data.frame(iso_vc_metadata_disc,bin_vc_metadata_disc)
}

#3.ggplot Theme####
{
  #Mytheme####
  Col=brewer.pal(3,"Set1")
  mytheme <- theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(color="black",size=12),
          strip.background = element_blank(),
          axis.title = element_text(size=15,color="black"),
          strip.text = element_text(color="black",size=15)
    )
  
  mytheme2 <- theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(color="black",size=12+3),
          strip.background = element_blank(),
          axis.title = element_text(size=15+3,color="black"),
          strip.text = element_text(color="black",size=15+3)
    )
}
#3.Variable searches####
{
  host_stat_cont.melt_cor <- host_stat_cont.melt %>%
    filter(Cont_Var_value != 0) %>%
    group_by(Type,Cont_Vars,Target_Vars) %>%
    do(cor_res=cal_cor(.$Target_Var_value,.$Cont_Var_value))
  host_stat_cont.melt_cor <- 
    cbind.data.frame(host_stat_cont.melt_cor,
                     list.rbind(host_stat_cont.melt_cor$cor_res))
  host_stat_cont.melt_cor %>% 
    filter(!is.na(Cor)) -> host_stat_cont.melt_cor
  host_stat_cont.melt_cor$sig <-
    ifelse(host_stat_cont.melt_cor$P<0.05,
           round(host_stat_cont.melt_cor$Cor,digits = 2),"")
  host_stat_cont.melt_cor %>% filter(Type == "Iso" & P < 0.05) %>% 
    select(Cont_Vars) %>% unique() -> sig_cont_var
  
host_stat_cont.melt_cor %>%
  ggplot(aes(x=Cont_Vars,y=Target_Vars,fill=Cor,label=sig))+
  geom_tile(color="black")+
  geom_text(size=2)+
  mytheme+
  scale_fill_gradient2(high="red",low="blue")+
  scale_x_discrete(expand = c(0,0,0,0))+
  scale_y_discrete(expand = c(0,0,0,0))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  labs(y="Specialisation d' ",x="")+
  facet_wrap(~Type,nrow=2) -> p_cont_d
ggsave("../../04.Figure/01.Raw/p_cont_d.tiff",
       plot=p_cont_d,device = "tiff",compression="lzw",
       width = 6.2*1.6,height = 4.8*1.6,units = "in")
{
fast_host_metadata.melt_cor <- fast_host_metadata.melt %>%
  filter(Cont_Var_value != 0) %>%
  group_by(Type,Cont_Vars,Net) %>%
  do(cor_res=cal_cor(.$PDI,.$Cont_Var_value))
fast_host_metadata.melt_cor <- 
  cbind.data.frame(fast_host_metadata.melt_cor,
                   list.rbind(fast_host_metadata.melt_cor$cor_res))
fast_host_metadata.melt_cor %>% 
  filter(!is.na(Cor)) -> fast_host_metadata.melt_cor
fast_host_metadata.melt_cor$sig <-
  ifelse(fast_host_metadata.melt_cor$P<0.05,
         round(fast_host_metadata.melt_cor$Cor,digits = 2),"")
#plot
fast_host_metadata.melt_cor %>%
  ggplot(aes(x=Cont_Vars,y=Net,fill=Cor,label=sig))+
  geom_tile(color="black")+
  geom_text(size=2)+
  mytheme+
  scale_fill_gradient2(high="red",low="blue")+
  scale_x_discrete(expand = c(0,0,0,0))+
  scale_y_discrete(expand = c(0,0,0,0))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  labs(y="PDI",x="")+
  facet_wrap(~Type,nrow=2) -> p_cont_pdi
ggsave("../../04.Figure/01.Raw/p_cont_pdi.tiff",
       plot=p_cont_pdi,device = "tiff",compression="lzw",
       width = 6.2*1.6,height = 4.8*1.6,units = "in")
}#PDI
{
  fast_host_metadata.melt_cor <- fast_host_metadata.melt %>%
    filter(Cont_Var_value != 0) %>%
    group_by(Type,Cont_Vars,Net) %>%
    do(cor_res=cal_cor(.$ND,.$Cont_Var_value))
  fast_host_metadata.melt_cor <- 
    cbind.data.frame(fast_host_metadata.melt_cor,
                     list.rbind(fast_host_metadata.melt_cor$cor_res))
  fast_host_metadata.melt_cor %>% 
    filter(!is.na(Cor)) -> fast_host_metadata.melt_cor
  fast_host_metadata.melt_cor$sig <-
    ifelse(fast_host_metadata.melt_cor$P<0.05,
           round(fast_host_metadata.melt_cor$Cor,digits = 2),"")
  #plot
  fast_host_metadata.melt_cor %>%
    ggplot(aes(x=Cont_Vars,y=Net,fill=Cor,label=sig))+
    geom_tile(color="black")+
    geom_text(size=2)+
    mytheme+
    scale_fill_gradient2(high="red",low="blue")+
    scale_x_discrete(expand = c(0,0,0,0))+
    scale_y_discrete(expand = c(0,0,0,0))+
    theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
    labs(y="Normalized degree",x="")+
    facet_wrap(~Type,nrow=2) -> p_cont_nd
  ggsave("../../04.Figure/01.Raw/p_cont_nd.tiff",
         plot=p_cont_nd,device = "tiff",compression="lzw",
         width = 6.2*1.6,height = 4.8*1.6,units = "in")
}#ND
{
  fast_host_metadata.melt_cor <- fast_host_metadata.melt %>%
    filter(Cont_Var_value != 0) %>%
    group_by(Type,Cont_Vars,Net) %>%
    do(cor_res=cal_cor(.$Strength,.$Cont_Var_value))
  fast_host_metadata.melt_cor <- 
    cbind.data.frame(fast_host_metadata.melt_cor,
                     list.rbind(fast_host_metadata.melt_cor$cor_res))
  fast_host_metadata.melt_cor %>% 
    filter(!is.na(Cor)) -> fast_host_metadata.melt_cor
  fast_host_metadata.melt_cor$sig <-
    ifelse(fast_host_metadata.melt_cor$P<0.05,
           round(fast_host_metadata.melt_cor$Cor,digits = 2),"")
  #plot
  fast_host_metadata.melt_cor %>%
    ggplot(aes(x=Cont_Vars,y=Net,fill=Cor,label=sig))+
    geom_tile(color="black")+
    geom_text(size=2)+
    mytheme+
    scale_fill_gradient2(high="red",low="blue")+
    scale_x_discrete(expand = c(0,0,0,0))+
    scale_y_discrete(expand = c(0,0,0,0))+
    theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
    labs(y="Species Strength",x="")+
    facet_wrap(~Type,nrow=2) -> p_cont_Strength
  ggsave("../../04.Figure/01.Raw/p_cont_Strength.tiff",
         plot=p_cont_Strength,device = "tiff",compression="lzw",
         width = 6.2*1.6,height = 4.8*1.6,units = "in")
}#Strength
{
  fast_host_metadata.melt_cor <- fast_host_metadata.melt %>%
    filter(Cont_Var_value != 0) %>%
    group_by(Type,Cont_Vars,Net) %>%
    do(cor_res=cal_cor(.$Nestedrank,.$Cont_Var_value))
  fast_host_metadata.melt_cor <- 
    cbind.data.frame(fast_host_metadata.melt_cor,
                     list.rbind(fast_host_metadata.melt_cor$cor_res))
  fast_host_metadata.melt_cor %>% 
    filter(!is.na(Cor)) -> fast_host_metadata.melt_cor
  fast_host_metadata.melt_cor$sig <-
    ifelse(fast_host_metadata.melt_cor$P<0.05,
           round(fast_host_metadata.melt_cor$Cor,digits = 2),"")
  #plot
  fast_host_metadata.melt_cor %>%
    ggplot(aes(x=Cont_Vars,y=Net,fill=Cor,label=sig))+
    geom_tile(color="black")+
    geom_text(size=2)+
    mytheme+
    scale_fill_gradient2(high="red",low="blue")+
    scale_x_discrete(expand = c(0,0,0,0))+
    scale_y_discrete(expand = c(0,0,0,0))+
    theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
    labs(y="Nestrank",x="")+
    facet_wrap(~Type,nrow=2) -> p_cont_nestrank
  ggsave("../../04.Figure/01.Raw/p_cont_nestrank.tiff",
         plot=p_cont_nestrank,device = "tiff",compression="lzw",
         width = 6.2*1.6,height = 4.8*1.6,units = "in")
}#Nestrank

}#cont_d heatmap
{
  Target_vars_use <- c("host_pc_R2","host_pc_R3",
                       "host_vc_R2", "host_vc_R3")
  host_stat_cont.melt_lm <- host_stat_cont.melt %>% 
    filter(Cont_Vars != "CPB" & 
           Cont_Var_value != 0   
          ) %>% 
    group_by(Type,Cont_Vars,Target_Vars) %>%
    do(
      data.frame(GenomeID=.$GenomeID,
                 Type=.$Type,
                 Cont_Vars=.$Cont_Vars,
                 Cont_Var_value=.$Cont_Var_value,
                 Target_Vars=.$Target_Vars,
                 Target_Var_value=.$Target_Var_value,
                 lm_P=cal_lm_p(.$Target_Var_value,
                               log(.$Cont_Var_value))[2]
       ))
  host_stat_cont.melt_lm$lm_sig <- 
    ifelse(host_stat_cont.melt_lm$lm_P <0.05,"sig","ns")

  host_stat_cont.melt_lm$lm_sig <- 
    factor(host_stat_cont.melt_lm$lm_sig,levels = c("sig","ns"))
  host_stat_cont.melt_lm %>% 
    filter(Cont_Vars %in% as.vector(sig_cont_var[,1])[-18] &
             Target_Vars %in% Target_vars_use) %>%
  ggplot(aes(x=log(Cont_Var_value),y=Target_Var_value))+
  geom_point(aes(colour=Type),alpha=0.3)+
  geom_smooth(aes(colour=Type,linetype=lm_sig),method="lm",se=F)+
  facet_grid(Target_Vars~Cont_Vars,scales = "free")+
  mytheme +
  scale_colour_manual(values=Col)+
  labs(y="Specialisation d' ",x="")+
  guides(colour=guide_legend(title = "")) -> p_sig_cont_d
  
ggsave("../../04.Figure/01.Raw/p_sig_cont_d.tiff",
       plot=p_sig_cont_d,device = "tiff",compression="lzw",
       width = 6.2*2*2,height = 4.8*1.3*2,units = "in")

}#cont_d point
#Host###
#4.1 OGT####
{
  {
    temperatures <- seq(15,55,5)
    iso_ogt_cut1 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% filter(OGT >= i & 
                                             !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut1 <- rbind.data.frame(iso_ogt_cut1,tmp)
    }
    iso_ogt_cut1 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(cor_res = cal_cor(.$host_vc_R2,
                           log(.$DoublingTime))) -> iso_ogt_cut1_cor
    iso_ogt_cut1_cor <- 
      cbind.data.frame(iso_ogt_cut1_cor,
                       list.rbind(iso_ogt_cut1_cor$cor_res))
    iso_ogt_cut1_cor$sig <- ifelse(iso_ogt_cut1_cor$P < 0.05,"sig","ns")
    iso_ogt_cut1_cor %>% 
      ggplot(aes(x=OGT_cut,y=Cor))+
      geom_point(aes(size=Cor,colour=sig))+
      geom_errorbar(
        aes(ymin = Cor-StdErr, ymax = Cor+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      scale_x_continuous(breaks = temperatures,
                         labels = c("[15,)","[20,)","[25,)","[30,)","[35,)",
                                    "[40,)","[45,)","[50,)","[55,)"
                         ))+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
      labs(x="OGT", 
           y="Pearson's r") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Cor+StdErr+0.1,label=N)) +
      guides(size=F,colour=F)  -> p_slop_ogt_d_1
    #Slope1
    temperatures <- seq(15,55,1)
    iso_ogt_cut1 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% filter(OGT >= i & 
                                             !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut1 <- rbind.data.frame(iso_ogt_cut1,tmp)
    }
    iso_ogt_cut1 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(lm_res = cal_lm_p(.$host_vc_R2,
                           log10(1/.$DoublingTime))) -> iso_ogt_cut1_lm
    iso_ogt_cut1_lm <- 
      cbind.data.frame(iso_ogt_cut1_lm,
                       list.rbind(iso_ogt_cut1_lm$lm_res))
    iso_ogt_cut1_lm$sig <- ifelse(iso_ogt_cut1_lm$P < 0.05,"sig","ns")
    iso_ogt_cut1_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr),
        position = position_dodge(.5), width = 0 
      )+
      geom_point(aes(size=Adj_R,fill=sig),color="black",pch=21)+
      #mytheme+
      mytheme2+
      theme(#axis.text.x = element_text(hjust=1,vjust=1),
        legend.position = c(0.4,0.4),
        legend.text = element_text(size=15),
        #legend.background = element_blank(),
        legend.background = element_rect(colour="white",color="black"),
        legend.title = element_text(size=15))+
      labs(y=expression("Slope"),
           x=expression("Optimal growth temperature (℃)")) +
      scale_fill_manual(values=c("white","black"))+
      scale_y_continuous(limits = c(-0.4,0.07))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.03,label=N),
                angle=90) +
      guides(fill=F,
             size=guide_legend(title = "Adjust\nR-square"))  -> p_slop_ogt_d_1.2
    #Fig.2C####
    Fig.2C <- p_slop_ogt_d_1.2
    iso_ogt_cut1_lm_fw <- iso_ogt_cut1_lm
    #Slop2
    temperatures <- seq(19,55,1)
    iso_ogt_cut1 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% filter(OGT <= i & 
                                             !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut1 <- rbind.data.frame(iso_ogt_cut1,tmp)
    }
    iso_ogt_cut1 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(lm_res = cal_lm_p(.$host_vc_R2,
                           log10(1/.$DoublingTime))) -> iso_ogt_cut1_lm
    iso_ogt_cut1_lm <- 
      cbind.data.frame(iso_ogt_cut1_lm,
                       list.rbind(iso_ogt_cut1_lm$lm_res))
    iso_ogt_cut1_lm$sig <- ifelse(iso_ogt_cut1_lm$P < 0.05,"sig","ns")
    iso_ogt_cut1_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr),
        position = position_dodge(.5), width = 0.1 
      )+
      geom_point(aes(size=Adj_R,fill=sig),pch=21)+
      #mytheme+
      mytheme2+
      theme(#axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        legend.position = c(0.8,0.4),
        legend.title = element_text(size=15),
        legend.background = element_rect(colour="white",color="black"),
        legend.text = element_text(size=15))+
      labs(y=expression("Slope"),
           x=expression("Optimal growth temperature (℃)")) +
      scale_fill_manual(values=c("white","black"))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.03,label=N),angle=90) +
      scale_y_continuous(limits = c(-0.4,0.07))+
      guides(fill=F,
             size=guide_legend(title = "Adjust\nR-square"))  -> p_slop_ogt_d_1.3
    #Fig.2D####
    Fig.2D <- p_slop_ogt_d_1.3
    iso_ogt_cut1_lm_bw <- iso_ogt_cut1_lm
    cowplot::plot_grid(p_slop_ogt_d_1.2,p_slop_ogt_d_1.3,nrow=2,
                       labels = c("A","B")) -> p_slpp_ogt_d.12_3
    ggsave("../../04.Figure/02.Publish/p_slpp_ogt_d.12_3.tiff",
           device = "tiff",plot = p_slpp_ogt_d.12_3,
           width = 6.6,height = 4.8*2,units = "in",dpi=300,
           compression = "lzw")
    #merge
    iso_ogt_cut1_lm_fw$type="Forward"
    iso_ogt_cut1_lm_bw$type="Backward"
    iso_ogt_cut1_lm <- rbind.data.frame(iso_ogt_cut1_lm_fw,
                                        iso_ogt_cut1_lm_bw)
    iso_ogt_cut1_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_point(aes(size=Adj_R,colour=sig))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1)
      )+
      labs(x="OGT", 
           y="Slope") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.01,label=N),
                angle=90,size=3) +
      facet_wrap(~type,nrow=2)+
      guides(colour=F,
             size=guide_legend(title = "Adjust\nR-square")) -> p_slop_bf
    
    ggsave("../../04.Figure/02.Publish/p_slop_bf.tiff",
           plot=p_slop_bf,device = "tiff",compression="lzw",
           width = 6.8,height = 4.8*2,units = "in")
    
    
    #Slop3
    temperatures <- seq(30,40,1)
    iso_ogt_cut1 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% 
        filter(OGT >= 30 & OGT <=40) %>%
        filter(OGT >= i & !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut1 <- rbind.data.frame(iso_ogt_cut1,tmp)
    }
    iso_ogt_cut1 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(lm_res = cal_lm_p(.$host_vc_R2,
                           log(.$DoublingTime))) -> iso_ogt_cut1_lm
    iso_ogt_cut1_lm <- 
      cbind.data.frame(iso_ogt_cut1_lm,
                       list.rbind(iso_ogt_cut1_lm$lm_res))
    iso_ogt_cut1_lm$sig <- ifelse(iso_ogt_cut1_lm$P < 0.05,"sig","ns")
    iso_ogt_cut1_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_point(aes(size=Adj_R,colour=sig))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
            legend.position = c(0.8,0.3))+
      labs(x="OGT", 
           y="Slope") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.01,label=N),angle=90) +
      guides(colour=F,
             size=guide_legend(title = "Adjust\nR-square")) +
      scale_x_continuous(limits = c(30,40))-> p_slop_ogt_d_1.4
    iso_ogt_cut1_lm_mfw <- iso_ogt_cut1_lm
    #Slop4
    temperatures <- seq(30,40,1)
    iso_ogt_cut1 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% 
        filter(OGT >= 30 & OGT <=40) %>%
        filter(OGT <= i & !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut1 <- rbind.data.frame(iso_ogt_cut1,tmp)
    }
    iso_ogt_cut1 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(lm_res = cal_lm_p(.$host_vc_R2,
                           log(.$DoublingTime))) -> iso_ogt_cut1_lm
    iso_ogt_cut1_lm <- 
      cbind.data.frame(iso_ogt_cut1_lm,
                       list.rbind(iso_ogt_cut1_lm$lm_res))
    iso_ogt_cut1_lm$sig <- ifelse(iso_ogt_cut1_lm$P < 0.05,"sig","ns")
    iso_ogt_cut1_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_point(aes(size=Adj_R,colour=sig))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
            legend.position = c(0.8,0.3))+
      labs(x="OGT", 
           y="Slope") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.01,label=N),angle=90) +
      guides(colour=F,
             size=guide_legend(title = "Adjust\nR-square")) +
      scale_x_continuous(limits = c(30,40))-> p_slop_ogt_d_1.5
    iso_ogt_cut1_lm_mbw <- iso_ogt_cut1_lm
    #merge
    iso_ogt_cut1_lm_mbw$type="Backward"
    iso_ogt_cut1_lm_mfw$type="Forward"
    iso_ogt_cut1_lm_m <- rbind.data.frame(iso_ogt_cut1_lm_mbw,
                                          iso_ogt_cut1_lm_mfw)
    iso_ogt_cut1_lm_m %>%
      ggplot(aes(x=OGT_cut,y=Coef))+
      geom_point(aes(size=Adj_R,colour=sig))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1)
      )+
      labs(x="OGT", 
           y="Slope") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.01,label=N),
                angle=90,size=3) +
      facet_wrap(~type)+
      guides(colour=F,
             size=guide_legend(title = "Adjust\nR-square")) +
      scale_x_continuous(breaks = seq(30,40,2))-> p_slop_mbf
    
    ggsave("../../04.Figure/02.Publish/p_slop_mbf.tiff",
           plot=p_slop_mbf,device = "tiff",compression="lzw",
           width = 6.2*2,height = 4.8*1.3,units = "in")
  }# method 1
  {
    temperatures <- seq(10,50,10)
    iso_ogt_cut2 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% 
        filter(OGT >= i & OGT < i+10 &
                 !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut2 <- rbind.data.frame(iso_ogt_cut2,tmp)
    }
    iso_ogt_cut2 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(cor_res = cal_cor(.$host_vc_R2,
                           log(.$DoublingTime))) -> iso_ogt_cut2_cor
    iso_ogt_cut2_cor <- 
      cbind.data.frame(iso_ogt_cut2_cor,
                       list.rbind(iso_ogt_cut2_cor$cor_res))
    iso_ogt_cut2_cor$sig <- ifelse(iso_ogt_cut2_cor$P < 0.05,"sig","ns")
    
    iso_ogt_cut2_cor %>% 
      ggplot(aes(x=OGT_cut,y=Cor)) +
      geom_point(aes(size=Cor,colour=sig))+
      geom_errorbar(
        aes(ymin = Cor-StdErr, ymax = Cor+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      scale_x_continuous(breaks = temperatures,
                         labels = c("[10,20)","[20,30)","[30,40)","[40,50)",
                                    "[50,60)"))+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
      labs(x="OGT", 
           y="Pearson's r") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Cor+StdErr+0.1,label=N)) +
      guides(size=F,colour=F)  -> p_slop_ogt_d_2
    
    #Slop
    temperatures <- seq(15,50,10)
    iso_ogt_cut2 <- NULL
    for(i in temperatures)
    {
      tmp <- iso_host_stat_cont %>% 
        filter(OGT >= i & OGT < i+10 &
                 !is.na(host_vc_R2))
      tmp$OGT_cut = i
      iso_ogt_cut2 <- rbind.data.frame(iso_ogt_cut2,tmp)
    }
    iso_ogt_cut2 %>% filter(!is.na(host_vc_R2)) %>% group_by(OGT_cut) %>%
      do(lm_res = cal_lm_p(.$host_vc_R2,
                           log(.$DoublingTime))) -> iso_ogt_cut2_lm
    iso_ogt_cut2_lm <- 
      cbind.data.frame(iso_ogt_cut2_lm,
                       list.rbind(iso_ogt_cut2_lm$lm_res))
    iso_ogt_cut2_lm$sig <- ifelse(iso_ogt_cut2_lm$P < 0.05,"sig","ns")
    
    iso_ogt_cut2_lm %>% 
      ggplot(aes(x=OGT_cut,y=Coef)) +
      geom_point(aes(size=Adj_R,colour=sig))+
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr,colour=sig),
        position = position_dodge(.5), width = 0.1 
      )+mytheme+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
      labs(x="OGT", 
           y="Slope") + 
      scale_colour_manual(values=c("grey",Col[1]))+
      geom_text(aes(x=OGT_cut,y=Coef+StdErr+0.1,label=N)) +
      guides(colour=F)  -> p_slop_ogt_d_2.2
  }# method 2
  {
    #10 window####
    window_size <- 10
    temperatures <- seq(ceiling(min(iso_host_stat_cont$OGT))-1,
                        ceiling(max(iso_host_stat_cont$OGT))-window_size,1)
    OGT_interval <- list()
    for(i in 1:length(temperatures))
    {
        iso_host_stat_cont %>% 
        filter(OGT >= temperatures[i] &
                 OGT < temperatures[i] + window_size) -> OGT_interval_tmp
      OGT_interval_tmp$OGT_interval_start <- temperatures[i]
      OGT_interval[[i]] <-  OGT_interval_tmp
    }
    names(OGT_interval) <- paste("OGT:",temperatures,
                                 "-",temperatures+window_size,sep="")
    OGT_interval_long <- ldply(OGT_interval,.id="OGT_interval")
    OGT_interval_long %>% 
      group_by(OGT_interval) %>%
      do(lm_res=cal_lm_p(.$host_vc_R2,
                         log10(1/.$DoublingTime))) -> OGT_interval_lm
    OGT_interval_lm <-
      cbind.data.frame(OGT_interval_lm[,1],ldply(OGT_interval_lm$lm_res))
    OGT_interval_lm$sig <- ifelse(OGT_interval_lm$P < 0.05,"sig","ns")
    OGT_interval_lm$OGT <- 
      as.numeric(str_extract(OGT_interval_lm$OGT_interval,"[0-9]+"))
    OGT_interval_lm %>%
      filter(OGT <= 60) %>%
      ggplot(aes(x=OGT,y=Coef)) +
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr),
        position = position_dodge(.5), width = 0 
      )+
      geom_point(aes(size=Adj_R,fill=sig),color="black",pch=21)+
      mytheme2+
      theme(#axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        legend.position = c(0.3,0.3),
        legend.title = element_text(size=15),
        legend.background = element_rect(colour="white",color="black"),
        legend.text = element_text(size=15),
        plot.title=element_text(hjust=0.5,size=18)
        )+
      labs(y=expression("Slope"),
           x=expression("Optimal growth temperature (°C)"),
           title="10°C window"
           )+
      scale_fill_manual(values=c("white","black"))+
      geom_text(aes(x=OGT,y=Coef+StdErr+0.03,label=N),angle=90,size=3.5) +
      scale_y_continuous(limits = c(-0.4,0.07))+
      scale_x_continuous(limits = c(10,60))+
      guides(fill=F,
             size=guide_legend(title = "Adjusted\nR-square"))  -> OGT_interval_lm.10
    #20 window####
    window_size <- 20
    temperatures <- seq(ceiling(min(iso_host_stat_cont$OGT))-1,
                        ceiling(max(iso_host_stat_cont$OGT))-window_size,1)
    OGT_interval <- list()
    for(i in 1:length(temperatures))
    {
      iso_host_stat_cont %>% 
        filter(OGT >= temperatures[i] &
                 OGT < temperatures[i] + window_size) -> OGT_interval_tmp
      OGT_interval_tmp$OGT_interval_start <- temperatures[i]
      OGT_interval[[i]] <-  OGT_interval_tmp
    }
    names(OGT_interval) <- paste("OGT:",temperatures,
                                 "-",temperatures+window_size,sep="")
    OGT_interval_long <- ldply(OGT_interval,.id="OGT_interval")
    OGT_interval_long %>% 
      group_by(OGT_interval) %>%
      do(lm_res=cal_lm_p(.$host_vc_R2,
                         log10(1/.$DoublingTime))) -> OGT_interval_lm
    OGT_interval_lm <-
      cbind.data.frame(OGT_interval_lm[,1],ldply(OGT_interval_lm$lm_res))
    OGT_interval_lm$sig <- ifelse(OGT_interval_lm$P < 0.05,"sig","ns")
    OGT_interval_lm$OGT <- 
      as.numeric(str_extract(OGT_interval_lm$OGT_interval,"[0-9]+"))
    OGT_interval_lm %>%
      filter(OGT <= 60) %>%
      ggplot(aes(x=OGT,y=Coef)) +
      geom_errorbar(
        aes(ymin = Coef-StdErr, ymax = Coef+StdErr),
        position = position_dodge(.5), width = 0 
      )+
      geom_point(aes(size=Adj_R,fill=sig),color="black",pch=21)+
      mytheme2+
      theme(#axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        legend.position = c(0.3,0.3),
        legend.title = element_text(size=15),
        legend.background = element_rect(colour="white",color="black"),
        legend.text = element_text(size=15),
        plot.title=element_text(hjust=0.5,size=18)
      )+
      labs(y=expression("Slope"),
           x=expression("Optimal growth temperature (°C)"),
           title="20°C window"
      )+
      scale_fill_manual(values=c("white","black"))+
      geom_text(aes(x=OGT,y=Coef+StdErr+0.03,label=N),angle=90,size=3.5) +
      scale_y_continuous(limits = c(-0.4,0.07))+
      scale_x_continuous(limits = c(10,60))+
      guides(fill=F,
             size=guide_legend(title = "Adjusted\nR-square"))  -> OGT_interval_lm.20
    
  }# method 3
  FIG.OGT_DT_d <- cowplot::plot_grid(p_slop_ogt_d_1,
                                     p_slop_ogt_d_2,nrow=2,labels = c("A","B")
  )
  ggsave("../../04.Figure/02.Publish/FIG.OGT_DT_d.tiff",device = "tiff",plot = FIG.OGT_DT_d,
         width = 6.6,height = 4.8*2,units = "in",dpi=300,
         compression = "lzw")
  tmp <- iso_host_stat_cont %>% filter(!is.na(host_vc_R2))
  tmp %>% ggplot()+
    geom_histogram(aes(OGT),bins=12,fill="grey",color="black")+
    mytheme+scale_y_continuous(expand = c(0,0,0,100)) +
    labs(y="The numer of genomes") -> p_hist_OGT
  ggsave("../../04.Figure/02.Publish/Fig_hist_OGT.tiff",
         device = "tiff",plot = p_hist_OGT,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")
  
  tmp <- iso_host_stat_cont
  tmp$OGT_cut <- ifelse(tmp$OGT >=40,"≥ 40","Others")
  tmp$OGT_cut <- factor(tmp$OGT_cut,levels = c("≥ 40","Others"))
  tmp %>% 
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2,
               colour=OGT_cut))+
    geom_point(alpha=0.5)+
    geom_smooth(method="lm",se=F)+
    mytheme+
    scale_colour_manual(values=c(Col[1],"#CCCCCC"))+
    labs(y="d'")+
    theme(legend.position = c(0.8,0.3),
          legend.background = element_blank(),
          legend.text = element_text(size=10))+
    guides(colour=guide_legend(title = "OGT")) -> Fig_p_OGT_DT_d
  ggsave("../../04.Figure/02.Publish/Fig_p_OGT_DT_d.tiff",
         device = "tiff",plot = Fig_p_OGT_DT_d,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")
}#Temperature range

{
  
  iso_host_stat_phen.melt %>% filter(!is.na(Phenotype_value)) %>%
    ggplot(aes(x=Phenotype_value,y=Target_Vars_value))+
    geom_boxplot(aes(fill=sig))+
    mytheme+
    facet_grid(Target_Vars~Phenotype,scales = "free")
}#phenotype_d boxplot
{
  Net_Use <- c("host_vc_R2","host_vc_R3")
  {
    slop_df <- iso_host_metadata_cont <-
      iso_host_metadata[,c(Old_VarNames_continuous_iso,
                           Old_VarNames_discrete_iso[-1],
                           Old_VarNames_phenotype_iso[-1]
                           )]
    colnames(slop_df) <- 
      c(New_VarNames_continuous,
        New_VarNames_discrete[-1],
        New_VarNames_phenotype[-1])
    slop_df <- merge.data.frame(
      slop_df,iso_d_host,by.x = "GenomeID",by.y = 0)
    slop_df <- melt.data.frame(
      slop_df,measure.vars = Target_vars,variable_name = "Net")
    colnames(slop_df)[colnames(slop_df)=="value"] <- "d"
    slop_df <- melt.data.frame(
      slop_df,measure.vars = c("Kingdom","Phylum","Class","Order","Family",
                               "Habitat","Genus",New_VarNames_phenotype[-1]),
      variable_name = "Discrete_Var")
    colnames(slop_df)[colnames(slop_df)=="value"] <- "Discrete_Var_value"
    slop_df$n <- 1
    slop_df %>% 
      filter(!is.na(d) & !is.na(Discrete_Var_value)) %>%
      group_by(Net,Discrete_Var,Discrete_Var_value) %>%
      do(data.frame(
        d=.$d,
        DoublingTime=.$DoublingTime,
        N=sum(.$n),
        OGT=.$OGT
      )) -> slop_df2
    slop_df2 %>% 
      filter(N > 8 & log10(DoublingTime) <= 2) %>% 
      group_by(Net,Discrete_Var,Discrete_Var_value) %>%
      do(lm_res=cal_lm_p(.$d,log10(1/.$DoublingTime))) -> slop_df3
    slop_df3 <- cbind.data.frame(
      slop_df3,list.rbind(slop_df3$lm_res))
    slop_df3$sig <- ifelse(slop_df3$P<0.05,"sig","ns")
    {
    slop_df3 %>% 
      filter(Net %in% Net_Use &
        Discrete_Var %in% c("Kingdom","Phylum","Class"
                            ) &
          Adj_R > 0
               ) %>%
      ggplot(aes(x=Discrete_Var_value,y=Coef))+
      geom_point(aes(colour=sig,size=Adj_R))+
      geom_errorbar( 
        aes(colour=sig,ymin = Coef+StdErr, ymax = Coef-StdErr),
        position = position_dodge(.5), width = 0.1 
      )+
      geom_hline(yintercept = 0,linetype="dashed")+
      facet_grid(Discrete_Var~Net,scales = "free_y")+
      coord_flip()+mytheme+
      labs(x="",y="Slop(d'~log(Doubling time))") -> p_slop_df3.1
    ggsave("../../04.Figure/01.Raw/p_slop_df3.taxnomy.tiff",
           plot=p_slop_df3.1,device = "tiff",compression="lzw",
           width = 6.2*2,height = 4.8*2,units = "in")
    
    slop_df3 %>% 
      filter(Discrete_Var == "Phylum" & 
               Net == "host_vc_R2" ) -> Phylum_use
    slop_df3 %>% 
      filter( Net=="host_vc_R2" & sig == "sig") -> all.sig
    slop_df2 %>%
      filter(
        Discrete_Var_value %in% Phylum_use$Discrete_Var_value &
          Discrete_Var == "Phylum" &
          Net=="host_vc_R2") ->slop_df2_phylum
    slop_df2_phylum$sig <- "ns"
    slop_df2_phylum$sig <- 
      ifelse(slop_df2_phylum$Discrete_Var_value %in% 
               all.sig$Discrete_Var_value,
             "sig","ns")
    
    slop_df2_phylum %>%
      filter(Discrete_Var_value %in% c(
        "Actinobacteria","Bacteroidetes",#"delta/epsilon subdivisions",
        #"Euryarchaeota",
        "Firmicutes","Proteobacteria"
        
      )) ->  slop_df2_phylum_abund
    slop_df2_phylum_abund$Discrete_Var_value <-
      as.character.factor(slop_df2_phylum_abund$Discrete_Var_value)
    slop_df2_phylum_abund$Discrete_Var_value[
      slop_df2_phylum_abund$Discrete_Var_value == "delta/epsilon subdivisions"
    ] <- "Delta/epsilon subdivisions"
    slop_df2_phylum_abund$Discrete_Var_value <-
      factor(slop_df2_phylum_abund$Discrete_Var_value,
             levels=c("Proteobacteria","Actinobacteria","Firmicutes",
                      "Bacteroidetes"))
    slop_df2_phylum_abund$sig <-
      factor(slop_df2_phylum_abund$sig,
             levels = c("sig","ns"))
    #slop_df2_phylum_abund %>%
    slop_df2_phylum$sig <- 
      factor(slop_df2_phylum$sig,levels=c("sig","ns"))
      slop_df2_phylum %>%
        filter( (log10(DoublingTime) < 1 &
                   Discrete_Var_value != "Deinococcus-Thermus"
                   )|
                 ( log10(DoublingTime) < 1 & 
                     Discrete_Var_value == "Deinococcus-Thermus" ) &
                 OGT > 30
                 ) %>%
    ggplot(aes(x=log10(1/DoublingTime),y=d))+
      geom_point(fill="grey",alpha=0.5)+
      facet_wrap(~Discrete_Var_value,scales = "free")+
      #scale_colour_gradient2(high="red",low="blue",midpoint = 45)+
      mytheme+
      geom_smooth(aes(linetype=sig),
                  method="lm",se=F,colour="#FF0033")+
        labs(y=expression("Specialization"*~italic("d'")),
             x=expression(Log[10]~"host growth rate (doublings/day)"))+
      guides(linetype=F) -> p_phylum_d_dt
    ggsave("../../04.Figure/01.Raw/p_phylum_d_dt.tiff",
           plot= p_phylum_d_dt,device = "tiff",compression="lzw",
           width = 6.8*1.4,height = 4.8,units = "in")
    ggsave("../../04.Figure/01.Raw/p_phylum_d_dt.pdf",
           plot=p_phylum_d_dt,
           width = 6.8*1.4,height = 4.8,units = "in")
    
    slop_df3 %>% 
      filter(sig=="sig" & Discrete_Var == "Class" & 
               Net == "host_vc_R2") -> class.sig
    slop_df2 %>%
      filter(Discrete_Var_value %in% 
               class.sig$Discrete_Var_value
             & Discrete_Var == "Class" &
               Net=="host_vc_R2") %>%
      ggplot(aes(x=log10(DoublingTime),y=d))+
      geom_point(aes(colour=OGT))+
      facet_wrap(~Discrete_Var_value,scales = "free")+
      scale_colour_gradient2(high="red",low="blue",midpoint = 45)+
      mytheme+
      geom_smooth(method="lm",se=F)+
      labs(y=expression(Specialisation~italic("d'")),
           x=expression(log[10]*italic(DT)))+
      guides(colour=guide_colorbar(title="OGT (℃)"))+
      theme(strip.text = element_text(face ="italic")) -> p_class_d_dt
    ggsave("../../04.Figure/01.Raw/p_class_d_dt.tiff",
           plot=p_class_d_dt,device = "tiff",compression="lzw",
           width = 6.8*1.2,height = 4.8*1.2,units = "in")
    
    slop_df3 %>% 
      filter(sig=="sig" & Discrete_Var == "Order" & 
               Net == "host_vc_R2") -> Order.sig
    slop_df2 %>%
      filter(Discrete_Var_value %in% 
               Order.sig$Discrete_Var_value 
             & Discrete_Var == "Order"
             & Net=="host_vc_R2"
             & log10(DoublingTime) < 2) %>%
      ggplot(aes(x=log10(1/DoublingTime),y=d))+
      geom_point(colour="black",alpha=0.5)+
      facet_wrap(~Discrete_Var_value,scales = "free")+
      #scale_colour_gradient2(high="red",low="blue",midpoint = 45)+
      mytheme+
      geom_smooth(method="lm",se=F,colour="#FF0033")+
      labs(y=expression(Specialization~italic("d'")),
           x=expression(log[10]*italic(G)))+
      theme(strip.text = element_text(face ="italic"))+
      guides(colour=guide_colourbar(title = "OGT (℃)"))-> p_Order_d_dt
    ggsave("../../04.Figure/01.Raw/p_Order_d_dt.tiff",
           plot=p_Order_d_dt,device = "tiff",compression="lzw",
           width = 6.8*1.3,height = 4.8*1.3,units = "in")
    
    slop_df3 %>% 
      filter(sig=="sig" & Discrete_Var == "Family" & 
               Net == "host_vc_R2") -> Family.sig
    slop_df2 %>%
      filter(Discrete_Var_value %in% 
               Family.sig$Discrete_Var_value
             & Discrete_Var == "Family"
             & Net=="host_vc_R2" 
             & Discrete_Var_value != ""
             ) %>%
      filter(Discrete_Var_value != "Alicyclobacillaceae" 
             #| 
            #   (Discrete_Var_value == "Alicyclobacillaceae" & OGT > 43)
             ) %>%
      ggplot(aes(x=log10(1/DoublingTime),y=d))+
      geom_point(fill="grey",alpha=0.5)+
      facet_wrap(~Discrete_Var_value,scales = "free")+
      #scale_colour_gradient2(high="red",low="blue",midpoint = 45)+
      mytheme+
      geom_smooth(method="lm",se=F,colour="#FF0033")+
      labs(y=expression("Host specialization"*~italic("d'")),
           x=expression("Host growth rate (log10)" ))+
      #theme(strip.text = element_text(face ="italic")) +
      guides(linetype=F,
             colour=guide_colorbar(title="OGT (℃)"))->
      p_Family_d_dt
    ggsave("../../04.Figure/01.Raw/p_Family_d_dt.tiff",
           plot=p_Family_d_dt,device = "tiff",compression="lzw",
           width = 6.8*1.3,height = 4.8*1.3,units = "in")
    
    slop_df3 %>% 
      filter(sig=="sig" & Discrete_Var == "Genus" & 
               Net == "host_vc_R2") -> Genus.sig
    slop_df2 %>%
      dplyr::filter(Discrete_Var_value %in% 
               Genus.sig$Discrete_Var_value
             & Discrete_Var == "Genus"
             & Net=="host_vc_R2" 
             & Discrete_Var_value != ""
      ) %>%
      dplyr::filter(Discrete_Var_value != "Alicyclobacillaceae" 
             #| 
             #   (Discrete_Var_value == "Alicyclobacillaceae" & OGT > 43)
      ) %>%
      ggplot(aes(x=log10(1/DoublingTime),y=d))+
      geom_point(fill="grey",alpha=0.5)+
      facet_wrap(~Discrete_Var_value,scales = "free")+
      #scale_colour_gradient2(high="red",low="blue",midpoint = 45)+
      mytheme+
      geom_smooth(method="lm",se=F,colour="#FF0033")+
      labs(y=expression("Virus specificity"~italic("d'")),
           x=expression(Log[10]~"host growth rate (doublings/day)"))+
      #theme(strip.text = element_text(face ="italic")) +
      guides(linetype=F,
             colour=guide_colorbar(title="OGT (℃)"))->
      p_Genus_d_dt
    ggsave("../../04.Figure/01.Raw/p_Genus_d_dt.tiff",
           plot=p_Genus_d_dt,device = "tiff",compression="lzw",
           width = 6.8,height = 4.8/1.5,units = "in")
    
  }#taxonomy 
    {
      slop_df3 %>% 
        filter(Net %in% c("host_pc_R2","host_pc_R3","host_pc_R4",
                          "host_vc_R2","host_vc_R3","host_vc_R4"
        ) &
          Discrete_Var %in% c("Cell shape","Motility",
                              "Sporulation","Temperature range","Salinity",
                              "Oxygen requirement","Gram strain"
          ) &
          Adj_R > 0
        ) %>%
        ggplot(aes(x=Discrete_Var_value,y=Coef))+
        geom_point(aes(colour=sig,size=Adj_R))+
        geom_errorbar( 
          aes(colour=sig,ymin = Coef+StdErr, ymax = Coef-StdErr),
          position = position_dodge(.5), width = 0.1 
        )+
        geom_hline(yintercept = 0,linetype="dashed")+
        facet_grid(Discrete_Var~Net,scales = "free_y")+
        coord_flip()+mytheme+
        labs(x="",y="Slop(d'~log(Doubling time))") -> p_slop_df3.2
      ggsave("../../04.Figure/01.Raw/p_slop_df3.phenotype.tiff",
             plot=p_slop_df3.2,device = "tiff",compression="lzw",
             width = 6.2*2.3,height = 4.8*2.3,units = "in")
    }#phenotype
    {
      slop_df3 %>% 
        filter(Net %in% c(
                          "host_vc_R2"
        ) &
          Discrete_Var %in% c("Habitat")
        #&
       #   Adj_R > 0
        ) -> slop_df3.3_tmp
      slop_df3.3_tmp$Discrete_Var_value <-
        str_replace(slop_df3.3_tmp$Discrete_Var_value,"_"," ")
      slop_df3.3_tmp$Adj_P <-
        p.adjust(slop_df3.3_tmp$P,
                 method = "BH",
                 n=length(slop_df3.3_tmp$P))
      slop_df3.3_tmp$sig <- 
        ifelse(slop_df3.3_tmp$Adj_P < 0.05,"sig","ns")
      slop_df3.3_tmp %>%
        filter(Discrete_Var_value!="Unknown") %>%
        ggplot(aes(x=Discrete_Var_value,y=Coef))+
        geom_errorbar( 
          aes(ymin = Coef+StdErr, ymax = Coef-StdErr),
          position = position_dodge(.5), width = 0
        )+
        geom_point(aes(fill=sig,size=Adj_R),pch=21)+
        geom_hline(yintercept = 0,linetype="dashed")+
        #facet_grid(Discrete_Var~Net,scales = "free_y")+
        coord_flip()+mytheme+
        labs(x="",
             y="Slope") +
        guides(fill=F,
               size=guide_legend(title="Adjust\nR-square")) +
        scale_fill_manual(values=c("white","black"))-> p_slop_df3.3
      ggsave("../../04.Figure/01.Raw/p_slop_df3.habitat.tiff",
             plot=p_slop_df3.3,device = "tiff",compression="lzw",
             width = 6.6,height = 4.8,units = "in")
      
    }#habitat
    {
      slop_df3 %>% filter(Net == "host_vc_R2" & sig=="sig") %>% 
        dplyr::select(Discrete_Var_value) %>% 
        unique() -> slop_sig
      slop_sig <- as.character.factor(slop_sig$Discrete_Var_value)
      dvv_levels <-  levels(slop_df2$Discrete_Var_value)
      dvv <-  slop_df2$Discrete_Var_value
      dvv_levels <- str_replace_all(dvv_levels,"_"," ")
      dvv <- str_replace_all(dvv,"_"," ")
      dvv <- factor(dvv,levels = dvv_levels)
      slop_df2$Discrete_Var_value <- dvv
      slop_sig2 <- str_replace_all(slop_sig,"_"," ")
      slop_df2 %>% filter(
        Discrete_Var %in% c("Habitat","Temperature range",
                            "Sporulation",
                            "Cell shape","Motility") &
          Net %in% c("host_vc_R2") &
          Discrete_Var_value %in%  slop_sig2 &
          Discrete_Var_value %in% 
          c("Thermal terrestrial","Sludge","Saline water","Saline sediment",
           "Plant","Nonsaline soil","Rod","Nonmotile","Motile","Thermophilic",
           "Mesophilic"
            )
      ) %>%
        ggplot(aes(x=log10(DoublingTime),y=d))+
        geom_point(aes(colour=OGT),alpha=0.6)+
        geom_smooth(method="lm",se=F,colour="red")+
        facet_wrap(~Discrete_Var_value,
                   scales = "free") +
        mytheme +
        labs(y="Specialisation d'")+
        scale_colour_gradient2(high="red",low="blue",
                               midpoint=45) -> p_point_smooth
      ggsave("../../04.Figure/01.Raw/p_point_smooth.tiff",
             plot=p_point_smooth,device = "tiff",compression="lzw",
             width = 6.2*1.5,height = 4.8*1.5,units = "in")
      
      #n.s
      
      slop_df2 %>% filter(
        Discrete_Var %in% c("Habitat","Temperature range",
                            "Sporulation",
                            "Cell shape","Motility") &
          Net %in% c("host_vc_R2") &
          !(Discrete_Var_value %in%  slop_sig2 &
              Discrete_Var_value %in% 
              c("Thermal terrestrial","Sludge","Saline water","Saline sediment",
                "Plant","Nonsaline soil","Rod","Nonmotile","Motile","Thermophilic",
                "Mesophilic","Unknown"
              )) &  !Discrete_Var_value %in% c("Fungus")
      ) %>%
        ggplot(aes(x=log10(DoublingTime),y=d))+
        geom_point(aes(colour=OGT),alpha=0.6)+
        geom_smooth(method="lm",se=F,colour="red",linetype="dashed")+
        facet_wrap(~Discrete_Var_value,
                   scales = "free") +
        mytheme +
        labs(y="Specialisation d'")+
        scale_colour_gradient2(high="red",low="blue",
                               midpoint=45) -> p_point_smooth.ns
      ggsave("../../04.Figure/01.Raw/p_point_smooth.ns.tiff",
             plot=p_point_smooth.ns,device = "tiff",
             compression="lzw",
             width = 6.2*1.7,height = 4.8*1.5,units = "in")
      
#habitat###
      slop_df2 %>% filter(
        Discrete_Var %in% c("Habitat") &
          Net %in% c("host_vc_R2") &
          Discrete_Var_value %in%  slop_sig2
      ) -> habitat_sig  
      habitat_sig$sig <- "sig"
      slop_df2 %>% filter(
        Discrete_Var %in% c("Habitat") &
          Net %in% c("host_vc_R2") &
          !(Discrete_Var_value %in% slop_sig2)
      ) -> habitat_ns  
      habitat_ns$sig <- "ns"
      habitat_sig <- rbind.data.frame(habitat_sig,habitat_ns)
      habitat_sig$Discrete_Var_value <-
        as.character.factor(habitat_sig$Discrete_Var_value)
      habitat_sig %>% 
        filter(!Discrete_Var_value %in% 
                 c("Unknown","Fungus","Thermal marine")
               ) -> habitat_sig
      habitat_sig$sig <- factor(
        habitat_sig$sig,levels=c("sig","ns"))
      habitat_sig$Discrete_Var_value <- factor(
        habitat_sig$Discrete_Var_value,levels = c(
          "Nonsaline sediment","Saline sediment",
          "Nonsaline soil","Saline soil","Nonsaline water",
          "Saline water","Nonsaline surface","Saline surface",
          "Thermal terrestrial","Animal","Plant",
          "Sludge","Wastewater","Mine area","Bioreactor"))
      habitat_sig$sig[habitat_sig$Discrete_Var_value 
                      == "Nonsaline surface"] <- "ns"
      habitat_sig$sig[habitat_sig$Discrete_Var_value 
                      == "Bioreactor"] <- "sig"
      habitat_sig  %>%
        ggplot(aes(x=log10(DoublingTime),y=d))+
        geom_point(aes(colour=OGT),alpha=0.6)+
        geom_smooth(aes(linetype=sig),method="lm",
                    se=F,colour="red")+
        facet_wrap(~Discrete_Var_value,ncol=5,
                   scales ="free_x") +
        mytheme +
        labs(y="Specialisation d'",
             x=expression(log[10]*DT))+
        scale_colour_gradient2(high="red",low="blue",
                               midpoint=45) +
        guides(linetype=F,
               colour=guide_colorbar(title="OGT(℃)")) + 
        scale_x_continuous(
          expand = c(0.1,0.1,0,0.3))-> p_point_habitat
      ggsave("../../04.Figure/01.Raw/p_point_habitat.tiff",
             plot=p_point_habitat,device = "tiff",
             compression="lzw",
             width = 6.6*1.5,height = 4.8*1.1,units = "in")
      Fig.1 <- p_point_habitat 
      
      slop_df2 %>% filter(
        Discrete_Var %in% c(
                            "Temperature range"
                            ) &
          Net %in% c("host_vc_R2") &
          Discrete_Var_value %in%  slop_sig2
      ) -> genotype_sig  
      genotype_sig$sig <- "sig"
      slop_df2 %>% filter(
        Discrete_Var %in% c(
                            "Temperature range"
        ) &
          Net %in% c("host_vc_R2") &
          !Discrete_Var_value %in%  slop_sig2
      ) -> genotype_ns  
      genotype_ns$sig <- "ns"
      genotype_sig <- rbind.data.frame(genotype_sig,genotype_ns)
      genotype_sig$sig <- factor(
        genotype_sig$sig,levels = c("sig","ns")
      )
      genotype_sig$Discrete_Var_value <- 
        factor(genotype_sig$Discrete_Var_value,
               levels = c("Thermophilic","Mesophilic","Psychrophilic") )
      genotype_sig %>%
        ggplot(aes(x=log10(1/DoublingTime),y=d,
                   group=Discrete_Var_value))+
        geom_point(#color="black",
                   aes(colour=Discrete_Var_value),alpha=0.6)+
        geom_smooth(aes(linetype=sig,
                        colour=Discrete_Var_value),method="lm",
                    se=F)+
        #facet_wrap(~Discrete_Var_value,ncol=3,
        #           scales = "free") +
        #mytheme +
        mytheme2 +
        labs(y=expression("Specialization"~italic("d'")),
             x=expression(Log[10]~"host growth rate (doublings/day)"))+
        scale_colour_manual(values=brewer.pal(5,"Set1")[c(1,2,5 )],
                            labels=c("Thermophiles","Mesophiles",
                                     "Psychrophiles")) +
        guides(linetype=F,
               colour=guide_legend(title="")) +
        theme(
          legend.position = c(0.83,0.2),
          #legend.background = element_blank(),
          legend.background = element_rect(colour="white",color="black"),
          legend.text = element_text(size=15)
        )-> p_point_genotype_tem
      ggsave("../../04.Figure/01.Raw/p_point_genotype_tem.tiff",
             plot=p_point_genotype_tem,device = "tiff",
             compression="lzw",
             width = 6.6,height = 4.8,units = "in")
      #Fig.2A####
      Fig.2A <- p_point_genotype_tem
      
      
      slop_df2 %>% filter(
        Discrete_Var %in% c("Motility","Sporulation",
                            "Salinity","Cell shape",
                            "Oxygen requirement","Gram strain"
        ) &
          Net %in% c("host_vc_R2") &
          Discrete_Var_value %in%  slop_sig2
      ) -> genotype_sig
      genotype_sig$sig <- "sig"
      slop_df2 %>% filter(
        Discrete_Var %in% c("Motility","Sporulation",
                            "Salinity","Cell shape",
                            "Oxygen requirement","Gram strain"
        ) &
          Net %in% c("host_vc_R2") &
          !Discrete_Var_value %in%  slop_sig2
      ) -> genotype_ns  
      genotype_ns$sig <- "ns"
      genotype_sig <- rbind.data.frame(genotype_sig,genotype_ns)
      genotype_sig$sig <- factor(
        genotype_sig$sig,levels = c("sig","ns")
      )
      genotype_sig$Discrete_Var_value <- 
        as.character.factor(genotype_sig$Discrete_Var_value)
      genotype_sig$Discrete_Var_value <-
        ifelse(genotype_sig$Discrete_Var_value == "NonHalophilic",
               "Nonhalophilic",genotype_sig$Discrete_Var_value)
      genotype_sig$Discrete_Var_value <- 
        factor(genotype_sig$Discrete_Var_value,
               levels=c("Motile","Nonmotile",
                        "Sporulating","Nonsporulating",
                        "Halophilic","Nonhalophilic",
                        "Aerobe","Anaerobe",
                        "Positive","Negative",
                        "Rod","Coccus","Filament",
                        "Spirilla","Vibrio"
                        ))
      genotype_sig %>% 
        filter( !Discrete_Var_value %in%
                  c("Spirilla","Vibrio","Rod","Coccus",
                    "Filament")) %>%
        ggplot(aes(x=log10(DoublingTime),y=d))+
        geom_point(aes(colour=OGT),alpha=0.6)+
        geom_smooth(aes(linetype=sig),method="lm",
                    se=F,colour="red")+
        facet_wrap(~Discrete_Var_value,ncol=5) +
        mytheme +
        labs(y="Specialisation d'",
             x=expression(log[10]*DT))+
        scale_colour_gradient2(high="red",low="blue",
                               midpoint=45) +
        guides(linetype=F,
               colour=guide_colorbar(title="OGT(℃)")) +
        theme(strip.text = element_text(size=10)) +
        scale_x_continuous(expand=c(0.1,0.1,0.1,0.1)
                           )-> p_point_genotype
      ggsave("../../04.Figure/01.Raw/p_point_genotype.tiff",
             plot=p_point_genotype,device = "tiff",
             compression="lzw",
             width = 6.6*1.5,height = 4.8,units = "in")
      Fig.2 <- p_point_genotype
      #Figure 1####
      f1_df <-rbind.data.frame(habitat_sig#,
                               #genotype_sig
                               )
      f1_df <- f1_df[f1_df$Discrete_Var!="Cell shape",]
      save(f1_df,file="../../04.Figure/01.Raw/figure1_df.Rdata")
      load("../../04.Figure/01.Raw/figure1_df.Rdata")
      f1_df$log10DT <- log10(1/f1_df$DoublingTime)
      lmer.fit1 <- lmer(d~log10DT+
                         (1+log10DT|Discrete_Var_value),
                       data = f1_df,REML = F)
      lmer.fit2 <- lmer(d~(1+log10DT|Discrete_Var_value),
                       data = f1_df,REML = F)
      anova(lmer.fit1,lmer.fit2,test="LRT")
      pred.mm <- ggpredict(lmer.fit1) 
      f1_df$sig <-
        ifelse(f1_df$Discrete_Var_value %in% 
                 c("Sludge","Bioreactor"),"ns",f1_df$sig)
      
      f1_df %>%
      #  filter( Discrete_Var = "Cell shape") %>%
        filter(log10(DoublingTime) < 1) %>%
        ggplot(aes(x=log10(1/DoublingTime),y=d))+
        #geom_point(aes(colour=OGT),alpha=0.6)+
        geom_point(fill="grey",alpha=0.5)+
        geom_smooth(aes(linetype=sig),method="lm",
                    se=T,colour="#FF0033")+
        facet_wrap(~Discrete_Var_value,ncol=5) +
        mytheme +
        labs(y=expression("Specialization"*~italic("d'")),
             x=expression(Log[10]~"host growth rate (doublings/day)"))+
        #scale_colour_gradient2(high="red",low="blue",
        #                       midpoint=45) +
        guides(linetype=F) +
        theme(strip.text = element_text(size=13),
              axis.text = element_text(size=11)) +
        scale_x_continuous(expand = c(0.02,0.02)) -> Figure1
      ggsave("../../04.Figure/01.Raw/Figure1.tiff",
             plot=Figure1,device = "tiff",
             compression="lzw",
             width = 6.6*1.5,height = 4.8*1.3,units = "in")
      
      ggsave("../../04.Figure/01.Raw/Figure1.pdf",
             plot=Figure1,device = "pdf",
             #compression="lzw",
             width = 6.6*1.5,height = 4.8*1.3)
  }#sig var point smooth
  }#reshaping data
  {
    iso_host_stat_cont
  }
}#iso: slop(d'~log(DoublingTime))

save(iso_host_stat_cont,file="../../02.Data/02.Rdata/iso_host_stat_cont.Rdata")

save(iso_host_stat_all,
     file="../../04.Figure/01.Raw/iso_host_stat_all.Rdata")
#4.2 AntiCrispr####
{
  iso_host_all_antiCrispr.melt <- 
    melt.data.frame(iso_host_stat_all,
                    measure.vars = New_VarNames_AntiCrispr[-1],
                    variable_name = "AntiCrisprType")
  colnames(iso_host_all_antiCrispr.melt)[
    colnames(iso_host_all_antiCrispr.melt)=="value"] <- "Num AntiCrispr"
  iso_host_all_antiCrispr.melt$AntiCrisprType <-
    str_replace(iso_host_all_antiCrispr.melt$AntiCrisprType,"AntiCRISPR","")
  iso_host_all_antiCrispr.melt$AntiCrisprType <-
    str_remove(iso_host_all_antiCrispr.melt$AntiCrisprType,"[|]")
  iso_host_all_antiCrispr.melt %>% ggplot(aes(
    x=`Num AntiCrispr`,y=host_vc_R2))+geom_point()+
    geom_smooth(method = "lm")+
    facet_wrap(~AntiCrisprType,scales = "free") +
    mytheme+
    labs(y="d'",x="Number of AntiCRISPR") -> p_antiCrispr
  ggsave("../../04.Figure/01.Raw/p_antiCrispr.tiff",
         plot=p_antiCrispr,device = "tiff",compression="lzw",
         width = 6.2*1.6,height = 4.8*1.6,units = "in")
  #Presence and absence
  AntiCrisprType <-
    iso_host_all_antiCrispr.melt$AntiCrisprType
  AntiCrisprType <- factor(AntiCrisprType)
  #iso_host_all_antiCrispr.melt2 <- iso_host_all_antiCrispr.melt
  iso_host_all_antiCrispr.melt$AntiCrisprType <-
    AntiCrisprType
  iso_host_all_antiCrispr.melt$AntiCrispr_Presence <- 
    ifelse(iso_host_all_antiCrispr.melt$`Num AntiCrispr` >0,"+","-")
  iso_host_all_antiCrispr.melt$AntiCrispr_Presence <- 
    ifelse(is.na(iso_host_all_antiCrispr.melt$AntiCrispr_Presence),"-",
           iso_host_all_antiCrispr.melt$AntiCrispr_Presence)
  iso_host_all_antiCrispr.melt %>%
    filter(log10(.$DoublingTime) < 2) %>%
    group_by(AntiCrisprType,AntiCrispr_Presence) %>%
    do(lm_res=cal_lm_p(.$host_vc_R2,log10(1/.$DoublingTime))) -> slop_df_anticrispr
  slop_df_anticrispr <- cbind.data.frame(
    slop_df_anticrispr,list.rbind(slop_df_anticrispr$lm_res))
  slop_df_anticrispr$Adj_P <- 
    p.adjust(slop_df_anticrispr$P, 
             method = "BH", 
             n = length(slop_df_anticrispr$P))
  slop_df_anticrispr$sig <- ifelse(slop_df_anticrispr$Adj_P<0.05,"sig","ns")
  slop_df_anticrispr$sig[is.na(slop_df_anticrispr$sig)] <- "ns"
  slop_df_anticrispr %>%
    filter(AntiCrispr_Presence == "+" &
    !AntiCrisprType %in% c("V-A,I-C","I-E,I-F") &
      Adj_R >=0) %>%
    ggplot(aes(x=AntiCrisprType,y=Coef,fill=sig))+
    geom_errorbar( 
      aes(ymin = Coef+StdErr, ymax = Coef-StdErr),
      position = position_dodge(.5), width = 0 
    )+
    geom_point(aes(size=Adj_R),colour="black",pch=21)+
    geom_hline(yintercept = 0,linetype="dashed")+
    coord_flip()+mytheme+
    scale_y_continuous(limits = c(-0.15,0.1))+
    labs(x="AntiCRISPR",y="Slope") +
    scale_fill_manual(values=c("white","black"))+
    guides(fill=F,
           size=guide_legend(title="Adjust\nR-square")
           )-> p_slop_anticrispr
  ggsave("../../04.Figure/01.Raw/p_slop_anticrispr.tiff",
         plot=p_slop_anticrispr,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
  iso_host_all_antiCrispr.melt %>% 
    filter(AntiCrispr_Presence != "-" &
            AntiCrisprType %in% c("I-C","II-A","II-C","VI-B")
             ) %>%
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2))+
    geom_point(aes(colour=OGT),alpha=0.6)+
    geom_smooth(method="lm",se=F,color="red")+
    facet_wrap(~AntiCrisprType,scales = "free_x") + mytheme +
    guides(colour=guide_colourbar(title = "OGT"))+
    labs(y="Specialisation d'") +
    scale_colour_gradient2(low="blue",high="red",
                           midpoint=45) -> p_d_GR_crispr
  ggsave("../../04.Figure/01.Raw/p_d_GR_crispr.tiff",
         plot=p_d_GR_crispr,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
  iso_host_all_antiCrispr.melt %>% 
    ggplot(aes(x=AntiCrispr_Presence,y=host_vc_R2))+
    geom_boxplot(alpha=0.3)+
    facet_wrap(~AntiCrisprType) + mytheme +
    labs(y="d'") -> p_d_boxplot_crispr
  ggsave("../../04.Figure/01.Raw/p_d_boxplot_crispr.tiff",
         plot=p_d_boxplot_crispr,device = "tiff",compression="lzw",
         width = 6.2*1.5,height = 4.8*1.5,units = "in")
  
  iso_host_all_antiCrispr.melt$AntiCrisprType <- 
    as.character.factor(iso_host_all_antiCrispr.melt$AntiCrisprType)
  
  iso_host_all_antiCrispr.melt$`Num AntiCrispr`[
    is.na(iso_host_all_antiCrispr.melt$`Num AntiCrispr`)] <- 0
  iso_host_all_antiCrispr.melt$AntiCrisprType <-
    ifelse(!iso_host_all_antiCrispr.melt$`Num AntiCrispr` >=1,
           "Absence",iso_host_all_antiCrispr.melt$AntiCrisprType)
  
  iso_host_all_antiCrispr.melt %>% 
    group_by(AntiCrisprType) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),
                       .$host_vc_R2)) -> iso_host_all_antiCrispr.cor
  iso_host_all_antiCrispr.cor <- iso_host_all_antiCrispr.cor[,-2]
  
  iso_host_all_antiCrispr.melt$AntiCrisprType2 <- "Others"
  iso_host_all_antiCrispr.melt$AntiCrisprType2 <-
    ifelse(grepl("I-C",iso_host_all_antiCrispr.melt$AntiCrisprType),
           "I-C",iso_host_all_antiCrispr.melt$AntiCrisprType2)
  iso_host_all_antiCrispr.melt$AntiCrisprType2 <-
    ifelse(grepl("VI-B",iso_host_all_antiCrispr.melt$AntiCrisprType),
           "VI-B",iso_host_all_antiCrispr.melt$AntiCrisprType2)
  iso_host_all_antiCrispr.melt %>% 
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2,colour=AntiCrisprType2))+
    geom_point()+geom_smooth(method="lm",se=F)+
    scale_colour_manual(values = rainbow(3))

  }#AntiCrispr
#4.3 Cold shock protein####
{
  iso_host_metadata_csp <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0
  )
  iso_host_metadata_csp$`Cold_shock_protein|7 & ScoF` <- 
    ifelse(!is.na(iso_host_metadata_csp$`Cold_shock_protein|7`) &
             !is.na(iso_host_metadata_csp$`Cold_shock_protein|ScoF`),
           1,NA)
  csp_id <- grep("Cold_shock_protein",colnames(iso_host_metadata_csp))
  iso_host_metadata_csp$num_csp <- 
    rowSums(iso_host_metadata_csp[,csp_id[1:16]],na.rm = T)
  
  iso_host_metadata_csp.melt <-
    melt.data.frame(iso_host_metadata_csp,
                    measure.vars = colnames(iso_host_metadata_csp)[csp_id],
                    variable_name = "csp_genes")
  iso_host_metadata_csp.melt$csp_genes <-
    str_remove(iso_host_metadata_csp.melt$csp_genes,
               "Cold_shock_protein[|]")
  
#  iso_host_metadata_csp.melt$csp_genes <- 
#    as.character.factor(iso_host_metadata_csp.melt$csp_genes)
  iso_host_metadata_csp.melt$csp_genes <-
    ifelse(iso_host_metadata_csp.melt$num_csp==0,"Absence",
           iso_host_metadata_csp.melt$csp_genes)
  iso_host_metadata_csp.melt <- 
    iso_host_metadata_csp.melt[
      !duplicated(
        iso_host_metadata_csp.melt[,c("Assembly.Accession","csp_genes")]),]
  iso_host_metadata_csp.melt$value <-
    ifelse(iso_host_metadata_csp.melt$num_csp==0,1,
           iso_host_metadata_csp.melt$value)
  iso_host_metadata_csp.melt %>% filter(!is.na(value) &
                                          !is.na(host_vc_R2)) %>% 
    group_by(csp_genes) %>% 
    do(cor_res = 
         cal_cor(log10(1/.$DoublingTime),
                 .$host_vc_R2)) -> iso_host_metadata_csp.melt_cor
  iso_host_metadata_csp.melt_cor <- 
    cbind.data.frame(iso_host_metadata_csp.melt_cor,
                     list.rbind(iso_host_metadata_csp.melt_cor$cor_res))
  iso_host_metadata_csp.melt_cor$Adj_P <-
    p.adjust(iso_host_metadata_csp.melt_cor$P,method="BH",
             n=length(iso_host_metadata_csp.melt_cor$P))
  
  iso_host_metadata_csp.melt_cor$sig <- 
    ifelse(iso_host_metadata_csp.melt_cor$Adj_P < 0.05,"sig","ns")
  csp_levels <- c("Absence","1","2","7","CapB","CspA","CspB","CspC","CspD",
                  "CspE","CspG","CspJ","CspLA","CspV","ScoF","YdfK",
                  "Unknown","7 & ScoF")
  iso_host_metadata_csp.melt_cor$csp_genes <-
    factor(iso_host_metadata_csp.melt_cor$csp_genes,
           levels = csp_levels)
  iso_host_metadata_csp.melt_cor %>%
    filter(csp_genes != "Unknown") -> tmp
  tmp %>%
    ggplot(aes(x=csp_genes,y=Cor))+
    geom_errorbar( 
      aes(ymin = Cor-StdErr, ymax = Cor+StdErr),
      position = position_dodge(.5), width = 0 
    )+
    geom_point(aes(fill=sig),fatten = abs(tmp$Cor)*6+1,pch=21)+
    geom_hline(yintercept = 0,linetype="dashed")+
    coord_flip()+
    #mytheme+
    mytheme2+
    scale_fill_manual(values = c("white","black"))+
    labs(x="CSP genes",
         y=expression("Pearson's"~italic(r)))+
    geom_text(aes(x=csp_genes,y=Cor+StdErr+0.05,label=N))+
    guides(size=F,fill=F) + 
    theme(axis.text.y = element_text(face="italic"),
          legend.text = element_text(size=15)
          ) + 
    scale_x_discrete(
      labels=c("Absence","csp1","csp2","csp7","capB","cspA","cspB","cspC","cspD",
               "cspE","cspG","cspJ","cspLA","cspV","scof","ydfK",
               "csp7 & scof")
    )-> p_cor_csp
  #Fig.2E####
  Fig.2E <- p_cor_csp
  ggsave("../../04.Figure/01.Raw/p_cor_csp.tiff",
         plot=p_cor_csp,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  iso_host_metadata_csp.melt_cor %>% 
    filter(sig=="sig") -> csp_sig
  
  iso_host_metadata_csp.melt$sig <- 
    ifelse(iso_host_metadata_csp.melt$csp_genes %in%
             as.character.factor(csp_sig$csp_genes),
           "sig","ns")
  iso_host_metadata_csp.melt$sig <- 
    factor(iso_host_metadata_csp.melt$sig,levels=c("sig","ns"))
  iso_host_metadata_csp.melt$csp_genes <- 
    factor(iso_host_metadata_csp.melt$csp_genes,
           levels = csp_levels )
  iso_host_metadata_csp.genes <- 
    as.character.factor(iso_host_metadata_csp.melt$csp_genes)
  iso_host_metadata_csp.levels <- 
    levels(iso_host_metadata_csp.melt$csp_genes)
  iso_host_metadata_csp.genes[
    grep("[0-9]",iso_host_metadata_csp.genes)] <-
    paste("CSP",iso_host_metadata_csp.genes[
      grep("[0-9]",iso_host_metadata_csp.genes)],sep="")
  iso_host_metadata_csp.levels[
    grep("[0-9]",iso_host_metadata_csp.levels)] <-
    paste("CSP",iso_host_metadata_csp.levels[
      grep("[0-9]",iso_host_metadata_csp.levels)],sep="")
  iso_host_metadata_csp.genes <- 
    factor(iso_host_metadata_csp.genes,
           levels=iso_host_metadata_csp.levels)
  iso_host_metadata_csp.melt2 <- iso_host_metadata_csp.melt
  iso_host_metadata_csp.melt2$csp_genes <- iso_host_metadata_csp.genes
  csp_factors <- levels(iso_host_metadata_csp.melt2$csp_genes)
  iso_host_metadata_csp.melt2$csp_genes <- 
    as.character.factor(iso_host_metadata_csp.melt2$csp_genes)
  
  iso_host_metadata_csp.melt2$csp_genes <- 
    factor(iso_host_metadata_csp.melt2$csp_genes,
           levels = csp_factors)
  
  csp_genes <- as.character.factor(
    unique(iso_host_metadata_csp.melt2$csp_genes))
  
  csp_genes.labs <- c("cspLA","Absence","cspA","cspB","cspD","cspC","Unknown",
                      "cspE","cspG","scoF","cspV","capB","cspJ","csp7","csp1",
                      "csp2","ydfK","csp7 & scofF")
  names(csp_genes.labs) <- csp_genes
  
  iso_host_metadata_csp.melt2 %>% 
    filter(!is.na(value) & log10(DoublingTime) < 1 &
             csp_genes != "Unknown"
             ) -> csp_gene.plot
  csp_genes <- as.character.factor(
    unique(csp_gene.plot$csp_genes))
  
  csp_genes.labs <- c("cspLA","Absence","cspA","cspB","cspD","cspC",
                      "cspE","cspG","scoF","cspV","capB","cspJ","csp7","csp1",
                      "csp2","ydfK","csp7 & scofF")
  names(csp_genes.labs) <- csp_genes
  csp_gene.plot$csp_genes <-
    factor(csp_gene.plot$csp_genes,
           levels = c("Absence",levels(csp_gene.plot$csp_genes)[-1]))
  
  csp_gene.plot %>%
    ggplot(aes(x=log10(1/DoublingTime),y=host_vc_R2))+
    geom_point(fill="grey",alpha=0.5)+
    geom_smooth(aes(linetype=sig),method="lm",se=F,colour="#FF0033")+
    facet_wrap(~csp_genes,scales = "free",
               labeller = labeller(csp_genes=csp_genes.labs)
               )+
    mytheme+
    #scale_colour_gradient2(high = "red",low="blue",midpoint = 45)+
    labs(y=expression("Specialization"*~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)"))+
    guides(linetype=F) +
    theme(strip.text = element_text(face="italic"))-> p_csp_point_smooth
  ggsave("../../04.Figure/01.Raw/p_csp_point_smooth.tiff",
         plot=p_csp_point_smooth,device = "tiff",compression="lzw",
         width = 6.6*2,height = 4.8*1.5,units = "in")
  ggsave("../../04.Figure/01.Raw/p_csp_point_smooth.pdf",
         plot=p_csp_point_smooth,device = "pdf",
         width = 6.6*1.5,height = 4.8*1.5,units = "in")
  
  iso_host_metadata_csp %>% 
    filter(!is.na(`Cold_shock_protein|7`) &
             !is.na(`Cold_shock_protein|ScoF`) &
             log10(DoublingTime) < 2
             ) %>%
    ggplot(aes(x=log10(1/DoublingTime),y=host_vc_R2))+
    geom_point(aes(colour=OGT))+
    geom_smooth(method = "lm",se=F)+
    scale_colour_gradient2(high = "red",low="blue",midpoint = 45)+
    labs(y="Specialisation d'",title = "CSP 7 & ScoF")+mytheme+
    theme(plot.title = element_text(hjust=0.5)) -> p_csp_7_scof
  ggsave("../../04.Figure/01.Raw/p_csp_7_scof.tiff",
         plot=p_csp_7_scof,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
             
  
  {
    {
      iso_host_metadata_csp %>% 
        filter(is.na(`Cold_shock_protein|7`) & 
                 is.na(`Cold_shock_protein|ScoF`) &
                 num_csp !=0) -> iso_host_metadata_csp.pos
      iso_host_metadata_csp %>% 
        filter(num_csp ==0) -> iso_host_metadata_csp.absence
      iso_host_metadata_csp.absence$csp_cut <- "Absence"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|7`)) -> iso_host_metadata_csp.7
      iso_host_metadata_csp.7$csp_cut <- "CSP 7"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|ScoF`)) -> iso_host_metadata_csp.ScoF
      iso_host_metadata_csp.ScoF$csp_cut <- "ScoF"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|7 & ScoF`)) -> iso_host_metadata_csp.dneg
      iso_host_metadata_csp.dneg$csp_cut <- "CSP 7 & ScoF"
    }#divide
    num_csp <- seq(1,9)
    iso_csp_cut1 <- NULL
    
    for(i in num_csp)
    {
      tmp <- iso_host_metadata_csp.pos %>% 
        filter(num_csp >= i &
                 !is.na(host_vc_R2))
      tmp$csp_cut = as.character(i)
      iso_csp_cut1 <- rbind.data.frame(iso_csp_cut1,tmp)
    }
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.absence)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.7)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.ScoF)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.dneg)
    
    iso_csp_cut1 %>% group_by(csp_cut) %>%
      do(cor_res=cal_cor(.$host_vc_R2,log(.$DoublingTime))) -> iso_csp_cut1_cor
    iso_csp_cut1_cor <- 
      cbind.data.frame(iso_csp_cut1_cor,
                       list.rbind(iso_csp_cut1_cor$cor_res))
    iso_csp_cut1_cor$sig <- 
      ifelse(iso_csp_cut1_cor$P < 0.05,"sig","ns")
    iso_csp_cut1_cor$csp_cut <-
      factor(iso_csp_cut1_cor$csp_cut,
             levels = c("CSP 7 & ScoF","CSP 7","ScoF","Absence",
                        as.character(seq(1,9))))
    iso_csp_cut1_cor %>%
    ggplot(aes(x=csp_cut,y=Cor))+
      geom_point(aes(colour=sig,size=abs(Cor)))+
      geom_errorbar( 
        aes(colour=sig,ymin = Cor-StdErr, ymax = Cor+StdErr),
        position = position_dodge(.5), width = 0.1 
      )+
      geom_hline(yintercept = 0,linetype="dashed")+
      coord_flip()+mytheme+
      scale_colour_manual(values = c("grey",Col[1]))+
      labs(x="")+
      geom_text(aes(x=csp_cut,y=Cor+StdErr+0.05,label=N))+
      guides(size=F,colour=F)+
      scale_x_discrete(labels=c(
        "CSP 7 & ScoF","CSP 7","ScoF","Absence",
        paste("[",seq(1,10),", )",sep=""))
      )+labs(x="Number of CSP",y="Pearson's r")+
      scale_y_continuous(expand = c(0,0.1))+
      theme(
        axis.title.y = 
          element_text(hjust = 0.7,vjust=-13)) -> p_num_csp_cor 
    ggsave("../../04.Figure/01.Raw/p_num_csp_cor.tiff",
           plot=p_num_csp_cor,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8,units = "in")
    iso_csp_cut1$csp_cut <- 
      factor(iso_csp_cut1$csp_cut,
             levels =  rev(c("CSP 7 & ScoF","CSP 7","ScoF","Absence",
                         as.character(seq(1,9))
      )))
    col_rain <- rainbow(13)
    col_rain[10] <- "black"
    iso_csp_cut1$sig <- 
      ifelse(iso_csp_cut1$csp_cut == "Absence","1","0")
    iso_csp_cut1  %>% 
      ggplot(aes(x=log10(DoublingTime),
                 y=host_vc_R2,
                 colour=csp_cut))+
      geom_point()+
      geom_smooth(aes(linetype=sig),method="lm",se=F)+
      scale_colour_manual(values = col_rain,
                          labels=c(paste("[",seq(1,9),", )",sep=""),
                                   "Absence","ScoF","CSP 7","CSP 7 & ScoF"
                                   ))+
      mytheme+guides(colour=guide_legend(title="",ncol=2),
                     linetype=F)+
      theme(legend.position = c(0.85,0.4),
            legend.background = element_blank(),
            legend.box.background = element_blank()) +
      labs(y=expression(Specialization~italic("d'")),
           x=expression(log[10]*italic(DT))
           ) -> p_point_smooth_num_csp
    #Fig.2F###
    #Fig.2F <- p_point_smooth_num_csp
    ggsave("../../04.Figure/01.Raw/p_point_smooth_num_csp.tiff",
           plot=p_point_smooth_num_csp,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8,units = "in")
    FIG.cps_DT_d1 <- cowplot::plot_grid(p_point_smooth_num_csp,
                                       p_num_csp_cor,nrow=2,
                                       labels = c("A","B"))
    ggsave("../../04.Figure/02.Publish/FIG.cps_DT_d1.tiff",
           plot=FIG.cps_DT_d1,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8*2,units = "in")
    
  }# method 1
  {
    {
      iso_host_metadata_csp %>% 
        filter(is.na(`Cold_shock_protein|7`) & 
                 is.na(`Cold_shock_protein|ScoF`) &
                 num_csp !=0) -> iso_host_metadata_csp.pos
      iso_host_metadata_csp %>% 
        filter(num_csp ==0) -> iso_host_metadata_csp.absence
      iso_host_metadata_csp.absence$csp_cut <- "Absence"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|7`)) -> iso_host_metadata_csp.7
      iso_host_metadata_csp.7$csp_cut <- "CSP 7"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|ScoF`)) -> iso_host_metadata_csp.ScoF
      iso_host_metadata_csp.ScoF$csp_cut <- "ScoF"
      iso_host_metadata_csp %>% 
        filter(!is.na(`Cold_shock_protein|7 & ScoF`)) -> iso_host_metadata_csp.dneg
      iso_host_metadata_csp.dneg$csp_cut <- "CSP 7 & ScoF"
    }#divide
    num_csp <- seq(1,10)
    iso_csp_cut1 <- NULL
    
    for(i in num_csp)
    {
      tmp <- iso_host_metadata_csp.pos %>% 
        filter(num_csp == i& 
                 !is.na(host_vc_R2))
      tmp$csp_cut = as.character(i)
      iso_csp_cut1 <- rbind.data.frame(iso_csp_cut1,tmp)
    }
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.absence)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.7)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.ScoF)
    iso_csp_cut1 <- 
      rbind.data.frame(iso_csp_cut1,iso_host_metadata_csp.dneg)
    
    iso_csp_cut1 %>% group_by(csp_cut) %>%
      do(cor_res=cal_cor(.$host_vc_R2,log(.$DoublingTime))) -> iso_csp_cut1_cor
    iso_csp_cut1_cor <- 
      cbind.data.frame(iso_csp_cut1_cor,
                       list.rbind(iso_csp_cut1_cor$cor_res))
    iso_csp_cut1_cor$sig <- 
      ifelse(iso_csp_cut1_cor$P < 0.05,"sig","ns")
    iso_csp_cut1_cor$csp_cut <-
      factor(iso_csp_cut1_cor$csp_cut,
             levels = c("CSP 7 & ScoF","CSP 7","ScoF","Absence",
                        as.character(seq(1,10))))
    iso_csp_cut1_cor %>% 
      filter(csp_cut %in% factor(seq(1,9))) %>%
      ggplot(aes(x=csp_cut,y=Cor))+
      geom_point(aes(colour=sig,size=abs(Cor)))+
      geom_errorbar( 
        aes(colour=sig,ymin = Cor-StdErr, ymax = Cor+StdErr),
        position = position_dodge(.5), width = 0.1 
      )+
      geom_hline(yintercept = 0,linetype="dashed")+
      coord_flip()+mytheme+
      scale_colour_manual(values = c("grey",Col[1]))+
      labs(x="")+
      geom_text(aes(x=csp_cut,y=Cor+StdErr+0.05,label=N))+
      guides(size=F,colour=F)+
      scale_x_discrete(labels=c(
        #"CSP 7 & ScoF","CSP 7","ScoF","Absence",
        seq(1,15))
      )+labs(x="Number of CSP",y="Pearson's r")+
      scale_y_continuous(expand = c(0,0.1)) -> p_num_csp_cor2
    ggsave("../../04.Figure/01.Raw/p_num_csp_cor2.tiff",
           plot=p_num_csp_cor2,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8,units = "in")
  }# method 2
  {
    iso_host_metadata_csp %>% 
      filter(num_csp != 0) -> iso_host_metadata_csp2
    iso_host_metadata_csp2 %>%
    select(c("Assembly.Accession",rev(csp_id)[-1])) -> 
      iso_host_metadata_csp_cluster_df
    rownames(iso_host_metadata_csp_cluster_df) <- 
      iso_host_metadata_csp_cluster_df$Assembly.Accession
    iso_host_metadata_csp_cluster_df <-  
      iso_host_metadata_csp_cluster_df[,-1]
    iso_host_metadata_csp_cluster_df[
      is.na(iso_host_metadata_csp_cluster_df)] <- 0
    dis="bray"
    hcl <- hclust(vegdist(iso_host_metadata_csp_cluster_df,method = dis))
    clusMember = cutree(hcl, 40)
    cl_count <- table(clusMember)
    
    iso_host_metadata_csp2$csp_cluster <- clusMember
    iso_host_metadata_csp2 %>% 
      filter(csp_cluster %in% which(cl_count>10) ) %>%
      group_by(csp_cluster) %>%
      do(cor_res=cal_cor(.$host_vc_R2,log(.$DoublingTime))) -> a
    cbind.data.frame(a,list.rbind(a$cor_res))[,-2] %>% arrange(Cor)
    iso_host_metadata_csp2 %>% filter(csp_cluster ==17) -> b
    b[is.na(b)] <- 0
    colMeans(b[,rev(csp_id)[-1]])
    iso_host_metadata_csp2 %>% 
#      filter(csp_cluster==9) %>%
      ggplot(aes(x=log(DoublingTime),
                 y=host_vc_R2,
                 colour=as.character(csp_cluster)))+
      geom_point()+geom_smooth(method="lm",se=F)+mytheme+
      scale_colour_manual(values = rainbow(40))
    
    {#mantel
    iso_host_metadata_csp2 %>% 
      filter(!is.na(host_vc_R2) & 
              # `Cold_shock_protein|CspC` >=3 &
               `Cold_shock_protein|CspC` >=4 ) -> iso_host_metadata_csp2.2
    a <- iso_host_metadata_csp2.2[,rev(csp_id)[-1]]
    a[is.na(a)] <- 0
    mantel(vegdist(a),
           vegdist(iso_host_metadata_csp2.2$host_vc_R2))
    #Mantel statistic r: 0.0232 
    #Significance: 0.001 
      }#mantel
  }#method 3
  {
    iso_host_metadata_csp2.melt <-
      melt.data.frame(iso_host_metadata_csp2,measure.vars = 
                        colnames(iso_host_metadata_csp2)[rev(csp_id)[-1]],
                      variable_name = "CSP_Genes")
    num_csp <- seq(1,10)
    iso_host_metadata_csp2.melt %>%
      filter(!is.na(value)& 
               !is.na(host_vc_R2)) -> iso_host_metadata_csp2.melt
    iso_host_metadata_csp2.melt2 <- NULL
    for(i in num_csp){
      tmp <- iso_host_metadata_csp2.melt %>% filter(value >=i)
      tmp$csp_cut <- i
      iso_host_metadata_csp2.melt2 <- 
        rbind.data.frame(iso_host_metadata_csp2.melt2,tmp)
    }
    iso_host_metadata_csp2.melt2 %>% 
      filter(!is.na(host_vc_R2)) %>%
      group_by(csp_cut,CSP_Genes) %>%
      do(cor_res=cal_cor(.$host_vc_R2,
                         log10(1/.$DoublingTime))) -> iso_host_metadata_csp2.melt2.Cor
    iso_host_metadata_csp2.melt2.Cor <- 
      cbind.data.frame(iso_host_metadata_csp2.melt2.Cor,
                       list.rbind(iso_host_metadata_csp2.melt2.Cor$cor_res))
    iso_host_metadata_csp2.melt2.Cor$CSP_Genes <-
      str_remove(iso_host_metadata_csp2.melt2.Cor$CSP_Genes,
               "Cold_shock_protein[|]")
    iso_host_metadata_csp2.melt2.Cor %>% 
      filter(CSP_Genes %in% c("CspC",
                              "CspD","CspLA") &
               N >= 10) %>% arrange(CSP_Genes) -> tmp
    tmp$CSP_Genes <- ifelse(tmp$CSP_Genes == "CspC","cspC",tmp$CSP_Genes)
    tmp$CSP_Genes <- ifelse(tmp$CSP_Genes == "CspD","cspD",tmp$CSP_Genes)
    tmp$CSP_Genes <- ifelse(tmp$CSP_Genes == "CspLA","cspLA",tmp$CSP_Genes)
    tmp$csp_cut <- as.numeric(tmp$csp_cut)
    tmp %>%
      ggplot(aes(x=as.factor(csp_cut),
                 y=Cor,color=CSP_Genes))+
      geom_pointrange(aes(ymin = Cor-StdErr, 
                          ymax = Cor+StdErr
                          ),fatten = abs(tmp$Cor)*6+1,
                      position=position_dodge(width=0.3))+
      #mytheme+
      mytheme2+
      labs(y=expression("Pearson's"~italic(r))
           ,x="Number of CSP genes")+
      guides(size=F,
             #colour=F,
             color=guide_legend(title = "")) + 
      theme(
        legend.position = c(0.15,0.28),
        #legend.background = element_blank(),
        legend.background = element_rect(colour="white",color="black"),
        legend.text = element_text(size=15,face="italic")
      )+
      scale_x_discrete(labels = 
                           c(paste("≥",seq(1,5),sep=""))) +
      scale_color_manual(values = c(
        "#CC3399","#660099","#FF6600")) +
      coord_flip() -> p_cor_num_csp
    ggsave("../../04.Figure/01.Raw/p_cor_num_csp.tiff",
           plot=p_cor_num_csp,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8,units = "in")
    
    p_num_csp_cor3 <- cowplot::plot_grid(
      p_cor_csp, p_cor_num_csp,labels = LETTERS[1:2],nrow=2
    )
    ggsave("../../04.Figure/02.Publish//p_num_csp_cor3.tiff",
           plot=p_num_csp_cor3,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8*2,units = "in")
    #Fig.2F####
    Fig.2F <- p_cor_num_csp
  }#method 
  iso_host_metadata_csp$OGT_cut <-
    ifelse(iso_host_metadata_csp$OGT >=40,"≥40","Others")
  iso_host_metadata_csp%>%
    filter(`Cold_shock_protein|CspB` >=3)->CspB
  iso_host_metadata_csp%>%
    filter(`Cold_shock_protein|CspC` >=3)->CspC
  iso_host_metadata_csp%>%
    filter(`Cold_shock_protein|CspD` >=3)->CspD
  iso_host_metadata_csp%>%
    filter(`Cold_shock_protein|CspLA` >=3)->CspLA
  iso_host_metadata_csp %>%
    filter(log10(DoublingTime) < 0.8) %>%
  ggplot(
         aes(x=log10(1/DoublingTime),
             y=host_vc_R2,colour=OGT_cut))+
    geom_point(alpha=0.3)+geom_smooth(method="lm",se=F)+
    #mytheme+
    mytheme2+
    labs(y=expression("Specialization"~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)"))+
    theme(legend.position = c(0.18,0.15),
          #legend.background = element_blank(),
          legend.background = element_rect(colour="white",color="black"),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15))+
    guides(colour=guide_legend(title = "OGT(°C)"))+
    scale_colour_manual(values=c("#E41A1C","#CCCCCC"))->
    Fig_p_OGT_DT_d2
  #  geom_smooth(data=CspC,aes(log10(DoublingTime),
  #                            y=host_vc_R2),
  #              colour="orange",method="lm",se=F)+
  #  geom_point(data=CspC,aes(log10(DoublingTime),
  #                           y=host_vc_R2),
  #             colour="orange")
  #geom_smooth(data=CspD,aes(log10(DoublingTime),
  #                          y=host_vc_R2),
  #            colour="blue",method="lm",se=F)+
  #  geom_point(data=CspD,aes(log10(DoublingTime),
  #                           y=host_vc_R2),
  #             colour="blue")+
  #  geom_smooth(data=CspB,aes(log10(DoublingTime),
  #                            y=host_vc_R2),
  #              colour="green",method="lm",se=F)+
  #  geom_point(data=CspB,aes(log10(DoublingTime),
  #                           y=host_vc_R2),
  #             colour="green")+
  #  geom_smooth(data=CspLA,aes(log10(DoublingTime),
  #                            y=host_vc_R2),
  #              colour="black",method="lm",se=F)+
  #  geom_point(data=CspLA,aes(log10(DoublingTime),
  #                          y=host_vc_R2),
  #            colour="black")+
  #  geom_text(x=-0.65,y=0.53,label="CspB",colour="green")+
  #  geom_text(x=-0.65,y=0.60,label="CspC",colour="orange")+
  #  geom_text(x=-0.65,y=0.30,label="CspD",colour="blue")+
  #  geom_text(x=-0.65,y=0.40,label="CspLA",colour="black")
   
  iso_host_metadata_csp %>%
    filter(log10(DoublingTime) < 1 &
             OGT_cut != "Others"
             ) %>%
    ggplot(
      aes(x=log10(1/DoublingTime),
          y=host_vc_R2,colour=OGT_cut))+
   # geom_point(alpha=0.3)+
    geom_smooth(method="lm",se=F)+mytheme+
    labs(y=expression("Spcialization"*~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)"))+
    theme(legend.position = c(0.88,0.15),
          legend.background = element_blank(),
          legend.text = element_text(size=10))+
    guides(colour=F)+
    scale_colour_manual(values=c("#E41A1C","#CCCCCC"))+
      geom_smooth(data=CspC,aes(log10(1/DoublingTime),
                                y=host_vc_R2),
                  colour="orange",method="lm",se=F)+
    geom_smooth(data=CspD,aes(log10(1/DoublingTime),
                              y=host_vc_R2),
                colour="blue",method="lm",se=F)+
    geom_smooth(data=CspB,aes(log10(1/DoublingTime),
                              y=host_vc_R2),
                colour="green",method="lm",se=F)+
    geom_smooth(data=CspLA,aes(log10(1/DoublingTime),
                             y=host_vc_R2),
                colour="black",method="lm",se=F)+
    geom_text(x=-0.42,y=0.7,label="CspB",colour="green")+
    geom_text(x=0.38,y=0.82,label="CspC",colour="orange")+
    geom_text(x=0.7,y=0.78,label="CspD",colour="blue")+
    geom_text(x=0.38,y=0.735,label="CspLA",colour="black")+
    geom_text(x=1.3,y=0.35,label="OGT ≥ 40°C",colour="red")->
    Fig.Psy_Ther
  ggsave("../../04.Figure/02.Publish/Fig.Psy_Ther.tiff",
         device = "tiff",plot = Fig.Psy_Ther,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")
  ggsave("../../04.Figure/02.Publish/Fig.Psy_Ther.pdf",
         device = "pdf",plot = Fig.Psy_Ther,
         width = 6.6,height = 4.8,units = "in")
   
   #Fig.2B####
  Fig.2B <- Fig_p_OGT_DT_d2
  ggsave("../../04.Figure/02.Publish/Fig_p_OGT_DT_d2.tiff",
         device = "tiff",plot = Fig_p_OGT_DT_d2,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")
  CspB$csp_type="CspB"
  CspC$csp_type="CspC"
  CspD$csp_type="CspD"
  CspLA$csp_type="CspLA"
  Csp_df <- rbind.data.frame(CspB,CspC)
  Csp_df <- rbind.data.frame(Csp_df,CspD)
  Csp_df <- rbind.data.frame(Csp_df,CspLA)
  Csp_df %>% ggplot()+
    geom_histogram(aes(x=OGT),fill="white",colour="black")+
    facet_wrap(~csp_type) +
    mytheme+
    scale_y_continuous(expand = c(0,0,0,5))+
    labs(y=" Number of genomes",x="OGT (℃)"
                   ) -> p_hist_ogt_csp
  ggsave("../../04.Figure/02.Publish/p_hist_ogt_csp.tiff",
         device = "tiff",plot = p_hist_ogt_csp,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")  
}#Cold shock protein

#Figure 2 ####
Fig.2 <- cowplot::plot_grid(OGT_interval_lm.10,OGT_interval_lm.20,
                            Fig.2B,Fig.2A,
                            Fig.2E,Fig.2F,
                            ncol=2,
                            labels = LETTERS[1:6],
                            label_size = 20)
ggsave("../../04.Figure/02.Publish/Fig.2.tiff",
       plot=Fig.2,device = "tiff",compression="lzw",
       width = 6.8*2.6,height = 4.8*2.6,units = "in")

ggsave("../../04.Figure/02.Publish/Fig.2.pdf",
       plot=Fig.2,device = "pdf",
       width = 6.8*2.5,height = 4.8*3,units = "in")

#4.4 Phage shock protein####
{
  iso_host_metadata_psp <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0
  )
  psp_id <- grep("Phage shock protein",colnames(iso_host_metadata_psp))
  iso_host_metadata_psp$num_psp <- 
    rowSums(iso_host_metadata_psp[,psp_id],na.rm=T)
  iso_host_metadata_psp %>% filter(num_psp>0)->iso_host_metadata_psp1
  iso_host_metadata_psp_cluster_df <- iso_host_metadata_psp1[,psp_id]
  rownames(iso_host_metadata_psp_cluster_df) <- 
    iso_host_metadata_psp1$Assembly.Accession
  iso_host_metadata_psp_cluster_df[
    is.na(iso_host_metadata_psp_cluster_df)] <- 0
  dis="bray"
  hcl <- hclust(vegdist(iso_host_metadata_psp_cluster_df,method = dis))
  clusMember = cutree(hcl, 10)
  cl_count <- table(clusMember)
  
  iso_host_metadata_psp1$psp_cluster <- clusMember
  iso_host_metadata_psp1 %>% 
    filter(!(is.na(host_vc_R2))) %>%
    group_by(psp_cluster) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),
                       .$host_vc_R2)) %>% as.data.frame()
  iso_host_metadata_psp1 %>% filter(psp_cluster == 2) -> a
  b <- a[,psp_id];b[is.na(b)] <- 0
  colMeans(b)
  
  iso_host_psp.melt <-
    melt.data.frame(iso_host_metadata_psp,
                    measure.vars = colnames(iso_host_metadata_psp)[psp_id],
                    variable_name = "PSP_Genes")
  iso_host_psp.melt$PSP_Genes <- 
    as.character.factor(iso_host_psp.melt$PSP_Genes)
  iso_host_psp.melt$PSP_Genes <-
    ifelse(iso_host_psp.melt$num_psp==0,
           "Absence",iso_host_psp.melt$PSP_Genes)
  iso_host_psp.melt %>%
    filter( !is.na(host_vc_R2) &
             (!is.na(value) & num_psp !=0) |
             (is.na(value) & num_psp ==0)) %>%
    group_by(PSP_Genes) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),
                       .$host_vc_R2)) -> iso_host_psp.melt_cor
  iso_host_psp.melt_cor <-
    cbind.data.frame(iso_host_psp.melt_cor,
                     list.rbind(iso_host_psp.melt_cor$cor_res))[,-2]
  iso_host_psp.melt_cor$PSP_Genes <- 
    str_remove(iso_host_psp.melt_cor$PSP_Genes,"Phage shock protein[|]")
  iso_host_psp.melt_cor$sig <-
    ifelse(iso_host_psp.melt_cor$P < 0.05,"sig","ns")
  iso_host_psp.melt_cor$PSP_Genes <-
    factor(iso_host_psp.melt_cor$PSP_Genes,
           levels = c("Absence","PspG","PspD","PspC","PspB","PspA"))
  iso_host_psp.melt_cor %>% 
    ggplot(aes(x=PSP_Genes,y=Cor))+
    geom_point(aes(colour=sig,size=Cor))+
    coord_flip()+mytheme+
    geom_errorbar(
      aes(ymin = Cor-StdErr, ymax = Cor+StdErr,colour=sig),
      position = position_dodge(.5), width = 0.1 
      )+
    scale_colour_manual(values = c("grey",Col[1]))+
    guides(size=F,colour=F)+
    labs(x="PSP genes",y="Pearson's r") -> p_psp_cor
  
  iso_host_psp.melt %>%
    filter( !is.na(host_vc_R2) &
              (!is.na(value) & num_psp !=0) |
              (is.na(value) & num_psp ==0)) -> iso_host_psp.melt
  
  iso_host_psp.melt$PSP_Genes <- 
    str_remove(iso_host_psp.melt$PSP_Genes,"Phage shock protein[|]")
  iso_host_psp.melt$sig <- "ns"
  iso_host_psp.melt$sig <- 
    ifelse(iso_host_psp.melt$PSP_Genes %in% c("Absence","PspA"),"sig","ns")
  iso_host_psp.melt$sig <- 
    factor(iso_host_psp.melt$sig,levels = c("sig","ns"))
  iso_host_psp.melt$PSP_Genes <- 
    factor(iso_host_psp.melt$PSP_Genes,
           levels = c("PspA","PspB","PspC","PspD","PspG","Absence"))
  iso_host_psp.melt %>% 
    ggplot(aes(x=log10(DoublingTime),y=host_vc_R2,colour=PSP_Genes))+
    geom_point()+mytheme+
    geom_smooth(aes(linetype=sig),method="lm",se=F)+
    scale_colour_manual(values = c(rainbow(5),"grey"))+
    guides(linetype=F,colour=guide_legend(title=""))+
    theme(legend.position = c(0.85,0.4),
          legend.background = element_blank())+
    labs(y=expression(Specialisation~italic("d'")),
         x=expression(log[10]*italic(DT))) -> p_psp_point
  ggsave("../../04.Figure/02.Publish/p_psp_point.tiff",
         device = "tiff",plot = p_psp_point,
         width = 6.6,height = 4.8,units = "in",dpi=300,
         compression = "lzw")
  
  
  FIG.psp_DT_d <- 
    cowplot::plot_grid(p_psp_point,p_psp_cor,
                       nrow=2,labels = c("A","B"))
  ggsave("../../04.Figure/02.Publish/FIG.psp_DT_d.tiff",
         device = "tiff",plot = FIG.psp_DT_d,
         width = 6.6,height = 4.8*2,units = "in",dpi=300,
         compression = "lzw")
  
  
  
  iso_host_metadata_psp %>% filter(num_psp == 0) -> psp_absence
  iso_host_metadata_psp %>% filter(num_psp != 0) -> psp_presence
  t.test(psp_absence$host_vc_R2,psp_presence$host_vc_R2)
}#Phage shock protein
#4.4 Acid and alkaline shock protein####
{
{
  iso_host_metadata_acid <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0
  )
  iso_host_metadata_acid$acid_presence <- 
    ifelse(!is.na(iso_host_metadata_acid$`Acid shock protein|Asp`),
           1,0)
  iso_host_metadata_acid %>% 
    filter(!is.na(host_vc_R2)) %>% group_by(acid_presence) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),.$host_vc_R2)) %>%
    as.data.frame()
  iso_host_metadata_acid$sig <-
    ifelse(iso_host_metadata_acid$acid_presence==0,"sig","ns")
  iso_host_metadata_acid$sig <-
    factor(iso_host_metadata_acid$sig,levels = c("sig","ns"))
  
  iso_host_metadata_acid %>%
    ggplot(aes(x=log10(DoublingTime),y=host_vc_R2,
               colour=as.factor(acid_presence)))+
    geom_point()+
    geom_smooth(aes(linetype=sig),method="lm",se=F)+
    mytheme+
    guides(linetype=F,
           colour=guide_legend(title="")
    )+
    theme(legend.position = c(0.8,0.4),
          legend.background = element_blank(),
          legend.text = element_text(size=12))+
    scale_colour_manual(values = c("grey",Col[1]),
                        labels=c("Absence","Acid shock protein"))+
    labs(y=expression(Specialisation~italic("d'")),
         x=expression(log[10]*italic(DT))) -> p_acid
}#Acid shock protein
{
  iso_host_metadata_alkaline <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0
  )
  iso_host_metadata_alkaline$alkaline_presence <- 
    ifelse(!is.na(iso_host_metadata_alkaline$`Alkaline shock protein|Asp23`),
           1,0)
  iso_host_metadata_alkaline %>% 
    filter(!is.na(host_vc_R2)) %>% group_by(alkaline_presence) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),.$host_vc_R2)) %>%
    as.data.frame()
  iso_host_metadata_alkaline$sig <-
    ifelse(iso_host_metadata_alkaline$alkaline_presence==0,"sig","ns")
  iso_host_metadata_alkaline$sig <-
    factor(iso_host_metadata_alkaline$sig,levels = c("sig","ns"))
  
  iso_host_metadata_alkaline %>%
    ggplot(aes(x=log10(DoublingTime),y=host_vc_R2,
               colour=as.factor(alkaline_presence)))+
    geom_point()+
    geom_smooth(aes(linetype=sig),method="lm",se=F)+
    mytheme+
    guides(linetype=F,
           colour=guide_legend(title="")
           )+
    theme(legend.position = c(0.8,0.4),
          legend.background = element_blank(),
          legend.text = element_text(size=12))+
    scale_colour_manual(values = c("grey",Col[1]),
                        labels=c("Absence","Alkaline shock protein"))+
    labs(y=expression(Specialisation~italic("d'")),
         x=expression(log[10]*italic(DT))) -> p_alkaline
  
}#Alkaline shock protein
  FIG.pH_DT_d <- 
    cowplot::plot_grid(p_acid,
                  p_alkaline,nrow=2,labels = c("A","B"))
  ggsave("../../04.Figure/02.Publish/FIG.pH_DT_d.tiff",
         plot=FIG.pH_DT_d,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8*2,units = "in")
}#pH
#4.5 Heat shock protein####
{
  iso_host_metadata_hsp <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0
  )
  iso_host_metadata_hsp %>% 
    filter(!is.na(host_vc_R2)) -> iso_host_metadata_hsp
  hsp_id <- grep("Heat shock protein",colnames(iso_host_metadata_hsp))
  iso_host_metadata_hsp$num_hsp <-
    rowSums(iso_host_metadata_hsp[,hsp_id],na.rm = T)
#  iso_host_metadata_hsp %>% 
#    filter(num_hsp >=1) -> iso_host_metadata_hsp
  iso_host_metadata_hsp.melt <-
    melt.data.frame(iso_host_metadata_hsp,
                    measure.vars = colnames(iso_host_metadata_hsp)[hsp_id],
                    variable_name = "HSP_genes")
  iso_host_metadata_hsp.melt %>%
    filter(!is.na(value)) %>%
    group_by(HSP_genes) %>%
    do(cor_res=cal_cor(log10(1/.$DoublingTime),
                       .$host_vc_R2)) %>% as.data.frame()
  iso_host_metadata_hsp.melt %>%
    filter(OGT >= 40) %>%
    ggplot(aes(x=value,y=OGT))+
    geom_point()+facet_wrap(~HSP_genes,scale="free")
  iso_host_metadata_hsp_com <- iso_host_metadata_hsp
  iso_host_metadata_hsp_com2 <- iso_host_metadata_hsp_com[
    !duplicated(iso_host_metadata_hsp$Assembly.Accession),
  ]
  rownames(iso_host_metadata_hsp_com2) <- 
    iso_host_metadata_hsp_com2$Assembly.Accession
  iso_host_metadata_hsp_com <- 
    iso_host_metadata_hsp_com[iso_host_metadata_hsp$num_hsp!=0,]
  iso_host_metadata_hsp_com_cl <- 
    iso_host_metadata_hsp_com[,hsp_id]
  iso_host_metadata_hsp_com_cl[is.na(iso_host_metadata_hsp_com_cl)] <-0
  dis="euclidean"
  hcl <- hclust(vegdist(iso_host_metadata_hsp_com_cl,method = dis,
                        binary =F
                        ))
  clusMember = cutree(hcl,30)
  cl_count <- table(clusMember)
  
  iso_host_metadata_hsp_com$hsp_cluster <- clusMember
  iso_host_metadata_hsp_com %>% group_by(hsp_cluster) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),
                       .$host_vc_R2)) -> iso_host_metadata_hsp.cor
  iso_host_metadata_hsp.cor <- 
    cbind.data.frame(iso_host_metadata_hsp.cor)
  
  colMeans(iso_host_metadata_hsp_com[
    iso_host_metadata_hsp_com$hsp_cluster==13,hsp_id],na.rm = T)
  iso_host_metadata_hsp_com %>%
    filter(hsp_cluster ==13) %>%
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2))+
    geom_point()+
    geom_smooth(method = "lm",se=F)
  
  iso_host_metadata_hsp_com %>% 
    filter(`Heat shock protein|HSP20` > 10 &
             `Heat shock protein|HSP40` > 100 &
             `Heat shock protein|HSP60` > 5 &
             `Heat shock protein|HSP70` > 15 &
             `Heat shock protein|HSP100` > 20
             )%>%
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2))+
    geom_point()+
    geom_smooth(method = "lm",se=F)
  
  iso_host_metadata_hsp.melt %>% filter(!is.na(value)) %>%
    group_by(HSP_genes) %>%
    do(cor_res=cal_cor(log(.$DoublingTime),.$host_vc_R2)) %>% 
    as.data.frame()
  
  iso_host_metadata_hsp %>%
    filter(`Heat shock protein|HSP60` >=20) %>%
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2))+
    geom_point()+
    geom_smooth(method = "lm",se=F)
  
  hsp_cor <- list()
  Seq_split <- list(
    seq(0,2,1),seq(0,2,1),seq(0,3,1),
    seq(0,2,1),seq(0,1,1),seq(0,4,1)#,seq(0,2,1)
  )
  hsp_id <- grep("Heat shock protein",
        colnames(iso_host_metadata_hsp))
  hsp_n <- colnames(iso_host_metadata_hsp)[hsp_id]
  names(hsp_n) <- hsp_id
  
  hsp_id <- as.numeric(names(rev(sort(hsp_n))))[-1]
  for(i in 1:length(hsp_id[1:6]))
  {
    tmp <- hsp_cor[[i]] <- windows_split(iso_host_metadata_hsp,
                                  hsp_id[i],Seq_split[[i]]) %>%
      group_by(cut_off) %>% 
      do(cor_res=cal_lm_p(log10(1/.$DoublingTime),.$host_vc_R2))
    tmp <- cbind.data.frame(tmp,list.rbind(tmp$cor_res))[,-2]
    hsp_cor[[i]] <- tmp
  }
  names(hsp_cor) <- colnames(iso_host_metadata_hsp)[hsp_id[1:6]]
  hsp_cor_long <- ldply(hsp_cor)
  hsp_cor_long$Adj_P <-
    p.adjust(hsp_cor_long$P,method="BH",
             n=length(hsp_cor_long$P))
  hsp_cor_long$sig <- ifelse(hsp_cor_long$Adj_P<0.05,"sig","ns")
  hsp_cor_long$.id <- 
    str_remove(hsp_cor_long$.id,"Heat shock protein[|]")
  hsp_cor_long$.id <- 
    factor(hsp_cor_long$.id,levels=c(
      "HSP20","HSP40","HSP60","HSP70","HSP90","HSP100","HSPIR"
    ))
  
  hsp_genes <- sort(unique(as.character.factor(hsp_cor_long$.id)))
  hsp_genes.labs <- c("hsp20","hsp40","hsp60","hsp70","hsp90","hsp100") 
  names(hsp_genes.labs) <- c("HSP20","HSP40","HSP60","HSP70","HSP90","HSP100")
  
  hsp_cor_long %>%
    #filter(N >= 10)%>%
    ggplot(aes(x=as.factor(cut_off),y=Coef))+
    geom_errorbar(
      aes(ymin = Coef-StdErr, ymax = Coef+StdErr),
      position = position_dodge(.5), width = 0.1 )+
    geom_point(aes(fill=sig),color="black",pch=21)+
    facet_wrap(~.id,ncol=2,
               labeller = labeller(.id=hsp_genes.labs)
               )+
    scale_fill_manual(values = c("white","black"))+mytheme+
    guides(fill=F,size=F)+
    labs(y="Slope",x="Number of HSP genes")+
    theme(strip.text = element_text(face="italic"))+
    geom_text(aes(x=as.factor(cut_off),
                    y=Coef+StdErr+2,label=N),angle=90,size=3) +
    scale_y_continuous(expand = c(0,0.5,0,0.5))-> p_hsp_cor 
  ggsave("../../04.Figure/01.Raw/p_hsp_cor.tiff",
         plot=p_hsp_cor,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
  hsp_id <- grep("Heat shock protein",
                 colnames(iso_host_metadata_hsp))
  iso_host_metadata_hsp$num_hsp <- 
    apply(iso_host_metadata_hsp[,hsp_id],1,sum,na.rm=T)
  hsp_id <- grep("Heat shock protein",
                 colnames(iso_host_metadata_hsp))
  hsp_n <- colnames(iso_host_metadata_hsp)[hsp_id]
  names(hsp_n) <- hsp_id
  hsp_n <- hsp_n[-7]
  hsp_id <- hsp_id[-7]
  hsp_df <- NULL
  for(i in 1:length(hsp_id)){
    tmp <- iso_host_metadata_hsp[
      iso_host_metadata_hsp[,hsp_id[i]]>=1,]
    tmp$hsp_type <- hsp_n[i]
    if(i==1){
      hsp_df <- tmp
    }
    hsp_df <- rbind.data.frame(hsp_df,tmp)
  }
  hsp_df %>%  
    group_by(hsp_type) %>%
    do(lm_res=cal_lm_p(.$host_vc_R2,log10(.$DoublingTime))) -> a
  a <- cbind.data.frame(a$hsp_type,ldply(a$lm_res))
  a$Adj_P <- p.adjust(a$P,method="BH",n=length(a$P))
  a %>% filter(Adj_P<0.05) -> b
  hsp_df$sig <- ifelse(hsp_df$hsp_type %in% b$`a$hsp_type`,"0","1")
  
  hsp_df$hsp_type <- as.character.factor(hsp_df$hsp_type)
  hsp_type.lev <- sort(unique(hsp_df$hsp_type))
  hsp_type.lev <- c(hsp_type.lev[-1],hsp_type.lev[1])
  hsp_df$hsp_type <- factor(hsp_df$hsp_type,
                            levels = rev(hsp_type.lev))
  hsp_df %>% 
    filter(log10(DoublingTime) < 1) %>%
    ggplot(aes(x=log10(1/DoublingTime),
               y=host_vc_R2,colour=hsp_type))+
    geom_point(alpha=0.7)+
    geom_smooth(aes(linetype=sig),method="lm",se=F)+
    mytheme+
    scale_colour_manual(
      values = c("#6c7197","#739211","#080300",
                 "#d92405","#3563eb","#eac124"),
                        labels=rev(c("hsp20","hsp40","hsp60",
                                 "hsp70","hsp90","hsp100")))+
    labs(y=expression("Specialization"*~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)"))+
    guides(linetype=F,
           colour=guide_legend(title="")
           )+
    theme(#legend.position = c(0.8,0.4),
          legend.background = element_blank(),
          legend.text = element_text(face="italic")) ->
    Figure_S1
  
  ggsave("../../04.Figure/01.Raw/Figure_S1.tiff",
         plot=Figure_S1,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  ggsave("../../04.Figure/01.Raw/Figure_S1.pdf",
         plot=Figure_S1,device = "pdf",
         width = 6.6,height = 4.8,units = "in")
  
  
}#Heat shock protein
#Cro/CI system
{
  iso_host_metadata_cro_ci <- merge.data.frame(
    iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0)
  iso_host_metadata_cro_ci %>% filter(!is.na(VirCI)) %>% 
    do(cor_res=cal_cor(log(.$DoublingTime),.$host_vc_R2)) %>% 
    as.data.frame()
  iso_host_metadata_cro_ci %>% filter(!is.na(VirCro)) %>% 
    do(cor_res=cal_cor(log(.$DoublingTime),.$host_vc_R2)) %>% 
    as.data.frame()
}

#4.6 CRISPRCas####
{
  iso_host_metadata_cas <- merge.data.frame(
    iso_host_metadata, iso_d_host,
    by.x="Assembly.Accession",by.y = 0,all.x = T)
  cas_id <- grep("CAS-[A-Za-z]+",colnames(iso_host_metadata_cas))
  cas_spacer_id <- grep("CAS-[A-Za-z]+_Spacer",colnames(iso_host_metadata_cas))
  cas_id <- cas_id[!cas_id %in% cas_spacer_id]
  tmp <- iso_host_metadata_cas[,cas_id]
  tmp[is.na(tmp)] <- 0
  n_cas <- rowSums(tmp,na.rm = T)
  tmp <- iso_host_metadata_cas[,cas_spacer_id]
  tmp[is.na(tmp)] <- 0
  n_spacer <- rowSums(tmp,na.rm = T)
  
  iso_host_metadata_cas$n_cas <- n_cas
  iso_host_metadata_cas$n_spacer <- n_spacer
  
  iso_host_metadata_cas$n_cas_ap <- 
    ifelse(iso_host_metadata_cas$n_cas == 0,
           "Absence","Presence")
  iso_host_metadata_cas %>%
    ggplot(aes(x=log(DoublingTime),
               y=host_vc_R2,colour=n_cas_ap))+
    geom_point()+
  geom_smooth(method="lm",se=F)+
    guides(colour=guide_legend(title="CRISPR-Cas system"))+
    mytheme+
    scale_colour_manual(values=c("grey","red"))+
    labs(y="Specialization d'")+
    theme(legend.position = c(0.85,0.25),
          legend.background = element_blank()
    ) -> p_Cas_DT_d
  ggsave("../../04.Figure/01.Raw/p_Cas_DT_d.tiff",
         plot= p_Cas_DT_d,device = "tiff",compression="lzw",
         width = 6.2,height = 4.8,units = "in")
  
    
    iso_host_metadata_cas %>% 
    group_by(n_cas_ap) %>% 
    do(lm_res=cal_lm_p(.$host_vc_R2,log(.$DoublingTime))) -> tmp
  cbind.data.frame(tmp$n_cas_ap,list.rbind(tmp$lm_res)[,c(1,3,2,5)])
    
    
  iso_host_metadata_cas.melt <- melt.data.frame(
    iso_host_metadata_cas,
    measure.vars = colnames(iso_host_metadata_cas)[cas_id],
    variable_name = "CasType"
  )
  iso_host_metadata_cas.melt %>% 
    filter(!is.na(value) & 
             !is.na(.$host_vc_R2)) %>%
    group_by(CasType) %>%
    do(cor_res=cal_cor(log10(1/.$DoublingTime),
                       .$host_vc_R2)) -> cas_cor
  cas_cor$Adj_P <- 
    p.adjust(cas_cor$P,method="BH",
             n=length(cas_cor$P))
  cas_cor <- cbind.data.frame(
    cas_cor$CasType,
    list.rbind(cas_cor$cor_res))
  iso_host_metadata_cas.melt %>% 
    filter(!is.na(value) & 
             !is.na(.$host_vc_R2)) -> cas_plot
  iso_host_metadata_cas.melt %>% 
    filter(
             !is.na(.$host_vc_R2) &
             n_cas==0
             ) -> cas_abs_plot
  cas_abs_plot <- cas_abs_plot[!duplicated(cas_abs_plot$Assembly.Accession),]
  cas_abs_plot$CasType <- "Absence"
  cas_plot <- rbind.data.frame(cas_plot,cas_abs_plot) 
  cas_plot$sig <- "ns"
  cas_plot$CasType <- str_remove(cas_plot$CasType,"pe")
  cas_plot$sig[cas_plot$CasType %in% c("Absence","CAS-IIIB")] <- "sig"
  cas_plot$sig <- factor(cas_plot$sig,levels = c("sig","ns"))
  
  CasType <- str_replace(cas_plot$CasType,"IU","I-U")
#  CasType <-  str_remove(CasType,"pe")
  CasType <- str_replace(CasType,"TypeU","U")
  CasType <- str_replace(CasType,"Unknown","U")
  CasType <- str_remove(CasType,"CAS-Type")
  CasType <- str_remove(CasType,"CAS-")
  for(i in LETTERS[1:6]){
    CasType <- str_replace(CasType,i,paste("-",i,sep=""))}
  CasType <- str_replace(CasType,"-Absence","Absence")
  CasType <- factor(CasType,levels = c(
    "I-A","I-B","I-C","I-D","I-E","I-F","I-U",
    "II-A","II-U",
    "III-A","III-B","III-U",
    "U","Absence"))
  
  cas_plot$CasType2 <- CasType
  
  cas_plot %>%
    filter(log10(DoublingTime) < 1) %>%
  ggplot(aes(x=log10(1/DoublingTime),y=host_vc_R2,colour=CasType2))+
    geom_point(alpha=0.15)+
    geom_smooth(aes(linetype=sig),method="lm",se=F)+
    scale_colour_manual(values = c(
      brewer.pal(12, "Paired"),"black"
      ,"grey"))+
    mytheme+
    labs(y=expression("Specialization"*~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)"))+
    guides(colour=guide_legend(title="CRISPR-Cas"),
           linetype=F)+
    theme(legend.position = c(0.15,0.5),
          legend.background = element_blank()
          ) -> p_CRISPRCas
  ggsave("../../04.Figure/02.Publish/p_CRISPRCas.tiff",
         plot=p_CRISPRCas,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  ggsave("../../04.Figure/02.Publish/p_CRISPRCas.pdf",
         plot=p_CRISPRCas,device = "pdf",
         width = 6.6,height = 4.8,units = "in")
  
  iso_host_metadata_cas %>%
    filter(n_spacer !=0) -> iso_host_metadata_cas.spacer
  iso_host_metadata_cas.spacer %>% 
    filter(OGT >=40)->iso_host_metadata_cas.spacer_h40
  iso_host_metadata_cas.spacer %>% 
    filter(OGT <40)->iso_host_metadata_cas.spacer_l40
  category <- c("≥40","Others")
  category <- factor(category,levels = category)
  Colour <- c("red","blue")
  df <- data.frame(category=category,Colour=Colour)
  df$x <- 3; df$y <- 100
  ggplot()+
    geom_point(data=df,aes(x=x,y=y,fill=category),pch=21,colour=NA)+
    geom_point(data=iso_host_metadata_cas.spacer_l40,
               aes(x=log(DoublingTime),y=log(n_spacer)),colour="grey")+
    geom_point(data=iso_host_metadata_cas.spacer_h40,
               aes(x=log(DoublingTime),y=log(n_spacer)),colour="red")+
#    geom_smooth(data=iso_host_metadata_cas.spacer,
#                aes(x=log(DoublingTime),y=n_spacer),
#                method = "lm",se=F,colour="black")+
    geom_smooth(data=iso_host_metadata_cas.spacer_l40,
                aes(x=log(DoublingTime),y=log(n_spacer)),
                method = "lm",se=F,colour="grey")+

    geom_smooth(data=iso_host_metadata_cas.spacer_h40,
                aes(x=log(DoublingTime),y=log(n_spacer)),
                method = "lm",se=F,colour="red")+
    mytheme+
    labs(y="Number of spacers",
         x="log(DoublingTime)")+
    guides(colour=guide_legend(title="OGT"))+
    scale_fill_manual(name="OGT",values = c("red","grey"))+
    theme(legend.position = c(0.8,0.6),
          legend.text = element_text(size=12)) +
    scale_y_continuous(limits = c(0,8))-> p_spacer_dt
  ggsave("../../04.Figure/01.Raw/p_spacer_dt.tiff",
         plot=p_spacer_dt,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
  iso_host_metadata_cas$OGT_cut <- 
    ifelse(iso_host_metadata_cas$OGT >=40,"≥40","Others")
  iso_host_metadata_cas %>%
    filter(n_spacer !=0)  %>%
    ggplot(aes(x=log(n_spacer),y=host_vc_R2,colour=OGT_cut))+
    geom_point()+geom_smooth(se=F,linetype="dashed")+
    labs(x="log(Spacers)",y="d'")+mytheme +
    scale_colour_manual(name="OGT",values = c("red","grey"))-> p_d_spacer
  ggsave("../../04.Figure/01.Raw/p_d_spacer.tiff",
         plot=p_d_spacer,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
}

iso_host_metadata_csp2 %>%
#  filter( !(!is.na(`Cold_shock_protein|7`) | 
#           !is.na(`Cold_shock_protein|ScoF`))) %>%
  select(DoublingTime,num_csp,`host_vc_R2`) -> sol
colnames(sol)[2] <- "num_csp"
sol %>% filter(!is.na(num_csp) & !is.na(host_vc_R2)) -> sol
sol$DoublingTime <- log10(sol$DoublingTime)
sol$num_csp <- sol$num_csp
sol %>% filter(DoublingTime < 4) -> sol
ordi<-ordisurf(sol[,1:2] ,
               sol$host_vc_R2,plot = T, bs="ds")
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ggplot(ordi.mite.na)+
   stat_contour(data = ordi.mite.na, 
  aes(x = x, y = y, z = z, fill = ..level..),
  geom="polygon", bins=1000)+
  labs(y="Number of CSP genes",
       x=expression(log[10]*italic(DT)))+
  scale_fill_gradientn(name="d'", colors=c("green", "blue", "yellow"))+
  #geom_point(data=sol,aes(x=DoublingTime,y=num_csp,colour=host_vc_R2))+
  mytheme -> p_heatmap_DT_CSP

ggsave("../../04.Figure/01.Raw/p_heatmap_DT_CSP.tiff",
       plot=p_heatmap_DT_CSP,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in")

#4.7 RM system####
{
    iso_host_metadata_RM <- merge.data.frame(
      iso_host_metadata, iso_d_host,
      by.x="Assembly.Accession",by.y = 0,all.x = T)
    iso_host_metadata_RM$RM_I[
      is.na(iso_host_metadata_RM$RM_I)] <- 0
    iso_host_metadata_RM$RM_II[
      is.na(iso_host_metadata_RM$RM_II)] <- 0
    iso_host_metadata_RM$RM_III[
      is.na(iso_host_metadata_RM$RM_III)] <- 0
    RMI <- iso_host_metadata_RM %>% 
      filter(RM_I != 0 & !is.na(host_vc_R2))
    RMII <- iso_host_metadata_RM %>% 
      filter(RM_II != 0 & !is.na(host_vc_R2))
    RMIII <- iso_host_metadata_RM %>% 
      filter(RM_III != 0 & !is.na(host_vc_R2))
    RMAbsence <- iso_host_metadata_RM %>% 
      filter(!(RM_I != 0 |  RM_II != 0 |
            RM_III != 0) & !is.na(host_vc_R2))
    RMI$RM="I"
    RMII$RM="II"
    RMIII$RM="III"
    RMAbsence$RM="Absence"
    RMI$Sig="ns"
    RMII$Sig="sig"
    RMIII$Sig="ns"
    RMAbsence$Sig="sig"
    RM <- rbind.data.frame(RMI,RMII,RMIII,RMAbsence)
    RM$Sig <- factor(RM$Sig,levels = c("sig","ns"))
    RM %>% group_by(RM) %>%
      filter(log10(DoublingTime) < 1) %>%
      do(lm_res=cal_lm_p(.$host_vc_R2,log10(1/.$DoublingTime))) -> RM_d_lm
    RM_d_lm_long <- ldply(RM_d_lm$lm_res)
    RM_d_lm_long$Adj_P <- 
      p.adjust(RM_d_lm_long$P,method="BH",
               n=length(RM_d_lm_long$P))
    RM %>% 
      filter(log10(DoublingTime) < 1) %>%
      ggplot(aes(x=log10(1/DoublingTime),
                      y=host_vc_R2,colour=RM)) + 
      geom_point(alpha=0.2)+
      geom_smooth(aes(linetype=Sig),se=F,method="lm",size=1)+
      mytheme +
      guides(linetype=F,
             colour=guide_legend(title="RM system"))+
      scale_colour_manual(values=c("grey","green","red","black"))+
      theme(legend.position = c(0.15,0.25),
            legend.background = element_blank())+
      labs(y=expression("Specialization"*~italic("d'")),
           x=expression(Log[10]~"host growth rate (doublings/day)")) -> p_RM_DT_d
    ggsave("../../04.Figure/01.Raw/p_RM_DT_d.tiff",
           plot=p_RM_DT_d,device = "tiff",compression="lzw",
           width = 6.6,height = 4.8,units = "in")
    ggsave("../../04.Figure/01.Raw/p_RM_DT_d.pdf",
           plot=p_RM_DT_d,device = "pdf",
           width = 6.6,height = 4.8,units = "in")
    
}

#4.8 Lytic switch####
{
  iso_host_metadata_LW <- merge.data.frame(
    iso_host_metadata, iso_d_host,
    by.x="Assembly.Accession",by.y = 0,all.x = T)
  LW_id <- which(colnames(iso_host_metadata_LW) %in% 
                   c("CI","Cro","fpsA","fpsR"))
  iso_host_metadata_LW[,LW_id][
    is.na(iso_host_metadata_LW[,LW_id])] <- 0
  iso_host_metadata_LW$CI_Cro <- 
    ifelse((iso_host_metadata_LW$CI > 0 & 
             iso_host_metadata_LW$Cro > 0),1,0)
  iso_host_metadata_LW$fpsA_R <- 
    ifelse(iso_host_metadata_LW$fpsA > 0 & 
             iso_host_metadata_LW$fpsR > 0,1,0)  
  iso_host_metadata_LW_Ab <- 
    iso_host_metadata_LW %>% filter(CI_Cro == 0 & fpsA_R == 0)
  iso_host_metadata_LW_CI <- 
    iso_host_metadata_LW %>% filter(CI != 0 )
  iso_host_metadata_LW_Cro <- 
    iso_host_metadata_LW %>% filter(Cro != 0 )
  iso_host_metadata_LW_CICro <- 
    iso_host_metadata_LW %>% filter(CI_Cro != 0 )
  iso_host_metadata_LW_fpsA <- 
    iso_host_metadata_LW %>% filter(fpsA != 0 )
  iso_host_metadata_LW_fpsR <- 
    iso_host_metadata_LW %>% filter(fpsR != 0 )
  iso_host_metadata_LW_fps <- 
    iso_host_metadata_LW %>% filter(fpsA_R != 0 )
  iso_host_metadata_LW_Ab$Type <- "Absence"
  iso_host_metadata_LW_CI$Type <- "CI"
  iso_host_metadata_LW_Cro$Type <- "Cro"
  iso_host_metadata_LW_fpsA$Type <- "FpsA"
  iso_host_metadata_LW_fpsR$Type <- "FpsR"
  iso_host_metadata_LW_CICro$Type <- "CI/Cro"
  iso_host_metadata_LW_all <- 
    rbind.data.frame(iso_host_metadata_LW_Ab,
                     iso_host_metadata_LW_CI,
                     iso_host_metadata_LW_Cro,
                     iso_host_metadata_LW_fpsA,
                     iso_host_metadata_LW_fpsR,
                     iso_host_metadata_LW_CICro
                     )
  iso_host_metadata_LW_all %>% group_by(Type) %>%
    do(lm_res=cal_lm_p(.$host_vc_R2,log(.$DoublingTime)))->a
  a <- cbind.data.frame(a$Type,list.rbind(a$lm_res))
  a$sig <- ifelse(a$P>=0.05,"ns","sig")
  a <- a[c(1,3,5,6),]
  a$Adj_P <- p.adjust(a$P,method="BH",n=length(a$P))
  iso_host_metadata_LW_all <- merge.data.frame(
    iso_host_metadata_LW_all,a,by.x = "Type",by.y = "a$Type"
  )
  iso_host_metadata_LW_all$sig <- 
    factor(iso_host_metadata_LW_all$sig,levels=c("sig","ns"))
  iso_host_metadata_LW_all %>% 
    filter(!Type %in% c("CI","Cro") &
             log10(DoublingTime) < 1
             ) %>%
    ggplot(aes(x=log10(1/DoublingTime),y=host_vc_R2,
               colour=Type))+
    geom_point(alpha=0.2)+
    geom_smooth(aes(colour=Type,linetype=sig),
                method="lm",se=F,
                size=1)+mytheme+
    guides(linetype=F,
           colour=guide_legend(title="Repressors"))+
    theme(legend.position = c(0.15,0.25),
          legend.background = element_blank()) + 
    scale_colour_manual(values = c(
      "#130b00","#f06e3d","#078e91","#a27ca1"
    )) +
    labs(y=expression("Specialization"*~italic("d'")),
         x=expression(Log[10]~"host growth rate (doublings/day)")) -> p_LW_DT_d
  ggsave("../../04.Figure/01.Raw/p_LW_DT_d.tiff",
         plot= p_LW_DT_d,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  ggsave("../../04.Figure/01.Raw/p_LW_DT_d.pdf",
         plot= p_LW_DT_d,device = "pdf",
         width = 6.6,height = 4.8,units = "in")
  
  
  csp <- colnames(iso_host_metadata_csp)[csp_id[1:16] ]
  csp_fpsR_com <- list()
  for(i in 1:length(csp))
  {
    tmp_csp <- iso_host_metadata_csp[,csp[i]]
    tmp_csp[is.na(tmp_csp)] <- 0
    tmp_csp <- ifelse(tmp_csp == 0,0,1)
    tmp_fpsR <- iso_host_metadata_csp$fpsR
    tmp_fpsR[is.na(tmp_fpsR)] <- 0
    tmp_fpsR <- ifelse(tmp_fpsR == 0,0,1)
    tmp_d <- iso_host_metadata_csp$host_vc_R2
    tmp_dt <- log(iso_host_metadata_csp$DoublingTime)
    tmp <- cbind.data.frame(tmp_csp,tmp_fpsR,tmp_d,tmp_dt)
    tmp <- tmp %>% filter(!is.na(tmp_d))
    csp_fpsR_com[[i]] <- aov(tmp_d~tmp_csp*tmp_fpsR+
                               tmp_csp:tmp_fpsR:tmp_dt,data=tmp)
  }
  aov_csp_fpsR <- lapply(csp_fpsR_com,function(x){
    summary_x <- summary(x)
    F_value <- summary_x[[1]]$`F value`[1:4]
    p_value <- summary_x[[1]]$`Pr(>F)`[1:4]
    item <- rownames(summary_x[[1]])[1:4]
    res <- data.frame(item,F_value,p_value)
    return(res)
  })
  names(aov_csp_fpsR) <- csp
  aov_csp_fpsR_df <- ldply(aov_csp_fpsR)
  aov_csp_fpsR_df$sig <- ""
  aov_csp_fpsR_df$sig <-
    ifelse(aov_csp_fpsR_df$p_value < 0.1,".",
           aov_csp_fpsR_df$sig)
  aov_csp_fpsR_df$sig <-
    ifelse(aov_csp_fpsR_df$p_value < 0.05,"*",
           aov_csp_fpsR_df$sig)
  aov_csp_fpsR_df$sig <-
    ifelse(aov_csp_fpsR_df$p_value < 0.01,"**",
           aov_csp_fpsR_df$sig)
  aov_csp_fpsR_df$sig <-
    ifelse(aov_csp_fpsR_df$p_value < 0.001,"***",
           aov_csp_fpsR_df$sig)
  aov_csp_fpsR_df$labels <- paste(
    format(aov_csp_fpsR_df$F_value,digits = 1),
    aov_csp_fpsR_df$sig,
    sep=""
  )
  aov_csp_fpsR_df$labels <- 
    str_remove(aov_csp_fpsR_df$labels,"[ ]+")
  aov_csp_fpsR_df$item <- 
    str_remove(aov_csp_fpsR_df$item,"tmp_")
  aov_csp_fpsR_df$.id <-
    str_remove(aov_csp_fpsR_df$.id,"Cold_shock_protein[|]")
  aov_csp_fpsR_df$item <- 
    str_remove_all(aov_csp_fpsR_df$item,"[ ]+")
  aov_csp_fpsR_df %>% filter(.id != "Unknown" & 
                               item != "Residuals") ->
    aov_csp_fpsR_df2
  aov_csp_fpsR_df2 %>% dcast(.id~item) -> csp_fpsR_aov
  csp_fpsR_aov$.id[1:3] <- paste("CSP",csp_fpsR_aov$.id[1:3],sep="")
  colnames(csp_fpsR_aov) <- c("Cold shock protein","CSP",
                              "CSP:FpsR","CSP:FpsR:DT","FpsR"
                              )
  csp_fpsR_aov <- csp_fpsR_aov[,c("Cold shock protein","CSP",
                                  "FpsR","CSP:FpsR",
                                  "CSP:FpsR:DT")]
  write.csv(csp_fpsR_aov,
              "../../05.Table/02.Publish/csp_fpsR_aov.csv",
            row.names = F)
}

#4.9 Receptors####
{
  iso_host_metadata_recep <- merge.data.frame(
    iso_host_metadata, iso_d_host,
    by.x="Assembly.Accession",by.y = 0,all.x = T)
  Prot_Recep_id <- 
    colnames(iso_host_metadata_recep)[
      grep("Prot_Recep",colnames(iso_host_metadata_recep))]
  iso_host_metadata_recep$tot_recep <- 
    rowSums(iso_host_metadata_recep[,Prot_Recep_id],na.rm = T)
  iso_host_metadata_recep$Absence <- 
    ifelse(iso_host_metadata_recep$tot_recep==0,-1,0)
  
  
  iso_host_metadata_recep2 <- 
    melt.data.frame(iso_host_metadata_recep,
                    measure.vars = c(Prot_Recep_id,"Absence"))
  
  iso_host_metadata_recep2$value[is.na(iso_host_metadata_recep2$value)] <- 0
  
  iso_host_metadata_recep2$variable <- 
    str_remove(iso_host_metadata_recep2$variable,"Prot_Recep[|]")
  
  iso_host_metadata_recep2 %>%
    filter(value != 0 ) %>% group_by(variable) %>%
    do(lm_res=cal_lm_p(.$host_vc_R2,log10(1/.$DoublingTime))) ->
  iso_host_metadata_recep2.lm
  iso_host_metadata_recep2.lm <- 
    cbind.data.frame(iso_host_metadata_recep2.lm,
                     list.rbind(iso_host_metadata_recep2.lm$lm_res))
  iso_host_metadata_recep2.lm <- iso_host_metadata_recep2.lm[,-2]
  iso_host_metadata_recep2.lm %>% arrange(P) -> iso_host_metadata_recep2.lm
   iso_host_metadata_recep2.lm %>% 
    arrange(variable) -> iso_host_metadata_recep2.lm
  Rcep_prot_anno <- 
    xlsx::read.xlsx("../../02.Data/01.Processing/HostReceptor/Receptor_gene.xlsx",
                    sheetIndex = 1)
  iso_host_metadata_recep2.lm <- 
    merge.data.frame(iso_host_metadata_recep2.lm,
                   Rcep_prot_anno[,-2],by.x="variable",
                   by.y="Gene.Symbol.of.viral.receptor",
                   all.x=T)
  colnames(iso_host_metadata_recep2.lm) <- 
    c("Gene","Coef","P","SD","Adj-R2","N","Receptor")
  iso_host_metadata_recep2.lm <- 
    iso_host_metadata_recep2.lm[,c("Gene","Receptor",
                                 "Coef","N","SD","Adj-R2","P")]
  write.table(iso_host_metadata_recep2.lm,
            "../../05.Table/01.Raw/Receptor_lm.tsv",row.names = F,
            quote=F,sep="\t")
  iso_host_metadata_recep2.lm.sig <- 
    iso_host_metadata_recep2.lm %>% filter(P < 0.05)
  iso_host_metadata_recep2 %>%
    filter(value != 0 & 
             variable %in% iso_host_metadata_recep2.lm.sig$Gene
           ) %>%
    ggplot(aes(x=log10(DoublingTime),y=host_vc_R2,colour=variable))+
    geom_point(alpha=0.5)+
    geom_smooth(method="lm",se=F)+
    mytheme+
    scale_colour_manual(
      values=c("grey","green","orange","blue","black"))+
    labs(y=expression(Specialisation~italic("d'")),
         x=expression(log[10]*italic(DT)))+
    guides(colour=guide_legend(title=""))+
    theme(legend.position = c(0.85,0.5)) -> fig.Recep_d_DT
  
  ggsave("../../04.Figure/01.Raw/fig.Recep_d_DT.tiff",
         plot= fig.Recep_d_DT,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
    
}

#5. Viral protein cluster####
{
  host_vrial_pc <- host_pc_M$iso_host_pc_R5
  iso_host_metadata_pc <- merge.data.frame(
    iso_host_metadata, iso_d_host,
    by.x="Assembly.Accession",by.y = 0)
  iso_host_metadata_pc <- merge.data.frame(
    iso_host_metadata_pc, host_vrial_pc,
    by.x="Assembly.Accession",by.y = 0)
  iso_host_metadata_pc %>% 
    filter(!is.na(.$host_vc_R2)) -> iso_host_metadata_pc
  lm_res <- lapply(colnames(host_vrial_pc),function(x){
    res <- iso_host_metadata_pc[iso_host_metadata_pc[,x] !=0,] %>%
      do(lm_res=cal_cor(log(.$DoublingTime),.$host_vc_R2))
    return(res)
  })
  lm_res2 <- lapply(lm_res,function(x){
    tmp <- x$lm_res[[1]];return(tmp)})
  names(lm_res2) <- colnames(host_vrial_pc)
  tmp <- lapply(lm_res2,length)
  lm_res3 <- lm_res2[unlist(lapply(tmp,function(x){x==4}))]
  ldply(lm_res3) %>% 
    filter(P < 0.05 & N > 10) %>% 
    arrange(Cor) -> lm_res3
  lm_res3 <- merge.data.frame(lm_res3,iso_pc_funcs,by.x=".id",
                              by.y = "pc")
  lm_res3 %>%　arrange(Cor) -> lm_res3
  iso_host_metadata_pc %>% filter(CNOHJPCP_04743!=0) %>%
    ggplot(aes(x=log(DoublingTime),y=host_vc_R2))+
    geom_point(aes(colour=OGT))+geom_smooth(method = "lm")
  iso_host_metadata_pc %>% 
    filter(CNOHJPCP_04743!=0) -> iso_host_metadata_pc.tmp
  Fig_p_OGT_DT_d2+
    geom_point(data=iso_host_metadata_pc.tmp,aes(x=log(DoublingTime),
                                                 y=host_vc_R2),colour="blue")+
    geom_smooth(data=iso_host_metadata_pc.tmp,aes(x=log(DoublingTime),
                                                  y=host_vc_R2),
                method = "lm",colour="blue",se=F)
  
}
#6. Net module####
{
  iso_host_metadata_mod <- merge.data.frame(
    iso_host_metadata, iso_d_host,
    by.x="Assembly.Accession",by.y = 0)
  iso_host_metadata_mod <- merge.data.frame(
    iso_host_metadata_pc, mod_members,
    by.x="Assembly.Accession",by.y = "GenomeID")
  iso_host_metadata_mod %>% group_by(ModuleID) %>% 
    do(cor_res=cal_cor(log(.$DoublingTime),
                       .$host_vc_R2)) ->iso_host_metadata_mod.cor 
  iso_mod.cor <- iso_host_metadata_mod.cor$cor_res
  names(iso_mod.cor) <- iso_host_metadata_mod.cor$ModuleID
  iso_mod.cor <- iso_mod.cor[!is.na(iso_mod.cor)]
  tmp <- lapply(iso_mod.cor,length)
  iso_mod.cor <- iso_mod.cor[tmp==4]
  ldply(iso_mod.cor) %>%
    filter(P < 0.05) %>%
    arrange(Cor) -> iso_mod.cor.df
  
  iso_host_metadata_mod %>%
    filter(ModuleID %in%
             #  rev(iso_mod.cor.df$.id)[1:4]
             c("56")
    ) -> iso_host_metadata_mod.tmp
  
  Fig_p_OGT_DT_d2+
    geom_point(data=iso_host_metadata_mod.tmp,aes(x=log(DoublingTime),
                                                  y=host_vc_R2),colour="blue")+
    geom_smooth(data=iso_host_metadata_mod.tmp,aes(x=log(DoublingTime),
                                                   y=host_vc_R2),
                method = "lm",colour="blue",se=F)
  
}

#7. frequency of Num VC####
{
  library(vcd)
  iso_host_metadata_freq_vc <- iso_host_metadata
  csp_id <- grep("Cold_shock_protein",
                 colnames(iso_host_metadata_freq_vc))
  iso_host_metadata_freq_vc$num_csp <- 
    rowSums(iso_host_metadata_freq_vc[,csp_id[1:16]],na.rm = T)
  iso_host_metadata_freq_vc$num_csp[
    is.na(iso_host_metadata_freq_vc$num_csp)
  ] <- 0
  
  psp_id <- grep("Phage shock protein",
                 colnames(iso_host_metadata_freq_vc))
  iso_host_metadata_freq_vc$num_psp <- 
    rowSums(iso_host_metadata_freq_vc[,psp_id],na.rm = T)
  iso_host_metadata_freq_vc$num_psp[
    is.na(iso_host_metadata_freq_vc$num_psp)
    ] <- 0
  
  iso_host_metadata_all <- iso_host_metadata_freq_vc 
  iso_host_metadata_all$CateType <- "All"
  iso_host_metadata_all$Presence <- "-"
    
  iso_host_metadata_cas %>% 
    filter(n_cas ==0) -> iso_host_metadata_nocas
  iso_host_metadata_nocas$CateType <- "CRISPR-Cas system"
  iso_host_metadata_nocas$Presence <- 0
  iso_host_metadata_cas %>% 
    filter(n_cas >0) -> iso_host_metadata_withcas
  iso_host_metadata_withcas$CateType <- "CRISPR-Cas system"
  iso_host_metadata_withcas$Presence <- 1
  
  iso_host_metadata_RM %>%
    filter(num_RM ==0) -> iso_host_metadata_norm
  iso_host_metadata_norm$CateType <- "RM system"
  iso_host_metadata_norm$Presence <- 0
  iso_host_metadata_RM %>%
    filter(num_RM >0) -> iso_host_metadata_withrm
  iso_host_metadata_withrm$CateType <- "RM system"
  iso_host_metadata_withrm$Presence <- 1
  
  iso_host_metadata_freq_vc %>% 
    filter(num_csp == 0) -> iso_host_metadata_nocsp
  iso_host_metadata_nocsp$CateType <- "Cold shock protein"
  iso_host_metadata_nocsp$Presence <- 0
  iso_host_metadata_freq_vc %>% 
    filter(num_csp > 0) -> iso_host_metadata_withcsp
  iso_host_metadata_withcsp$CateType <- "Cold shock protein"
  iso_host_metadata_withcsp$Presence <- 1
  
  iso_host_metadata_freq_vc %>% 
    filter(is.na(`Acid shock protein|Asp`)) -> iso_host_metadata_noacid
  iso_host_metadata_noacid$CateType <- "Acid shock protein"
  iso_host_metadata_noacid$Presence <- 0
  iso_host_metadata_freq_vc %>% 
    filter(`Acid shock protein|Asp` > 0) -> iso_host_metadata_withacid
  iso_host_metadata_withacid$CateType <- "Acid shock protein"
  iso_host_metadata_withacid$Presence <- 1
  
  iso_host_metadata_freq_vc %>% 
    filter(is.na(`Alkaline shock protein|Asp23`)) -> iso_host_metadata_noalkaline
  iso_host_metadata_noalkaline$CateType <- "Alkaline shock protein"
  iso_host_metadata_noalkaline$Presence <- 0
  iso_host_metadata_freq_vc %>% 
    filter(!is.na(`Alkaline shock protein|Asp23`)) -> iso_host_metadata_withalkaline
  iso_host_metadata_withalkaline$CateType <- "Alkaline shock protein"
  iso_host_metadata_withalkaline$Presence <- 1
  
  iso_host_metadata_freq_vc %>% 
    filter(num_psp ==0) -> iso_host_metadata_nopsp
  iso_host_metadata_nopsp$CateType <- "Phage shock protein"
  iso_host_metadata_nopsp$Presence <- 0
  iso_host_metadata_freq_vc %>% 
    filter(num_psp > 0) -> iso_host_metadata_withpsp
  iso_host_metadata_withpsp$CateType <- "Phage shock protein"
  iso_host_metadata_withpsp$Presence <- 1
  
  iso_host_metadata_freq_vc %>% 
    filter(OGT >=40) -> iso_host_metadata_L40
  iso_host_metadata_L40$CateType <- "OGT"
  iso_host_metadata_L40$Presence <- "≥ 40 ℃"
  iso_host_metadata_freq_vc %>% 
    filter(OGT < 40) -> iso_host_metadata_l40
  iso_host_metadata_l40$CateType <- "OGT"
  iso_host_metadata_l40$Presence <- "< 40 ℃"

  iso_host_metadata_freq_vc %>% 
    filter(DoublingTime < 2.5) -> iso_host_metadata_fg
  iso_host_metadata_fg$CateType <- "Growth rate"
  iso_host_metadata_fg$Presence <- "Fast"
  iso_host_metadata_freq_vc %>% 
    filter(DoublingTime >= 2.5) -> iso_host_metadata_sg
  iso_host_metadata_sg$CateType <- "Growth rate"
  iso_host_metadata_sg$Presence <- "Slow"
  
  
    selected_var <- c("num_VC","CateType","Presence")
  poisson <- rbind.data.frame(
    iso_host_metadata_all[,selected_var],
    iso_host_metadata_nocas[,selected_var],
    iso_host_metadata_withcas[,selected_var],
    iso_host_metadata_norm[,selected_var],
    iso_host_metadata_withrm[,selected_var],
    iso_host_metadata_nocsp[,selected_var],
    iso_host_metadata_withcsp[,selected_var],
    iso_host_metadata_noacid[,selected_var],
    iso_host_metadata_withacid[,selected_var],
    iso_host_metadata_noalkaline[,selected_var],
    iso_host_metadata_withalkaline[,selected_var],
    iso_host_metadata_nopsp[,selected_var],
    iso_host_metadata_withpsp[,selected_var],
    iso_host_metadata_L40[,selected_var],
    iso_host_metadata_l40[,selected_var],
    iso_host_metadata_fg[,selected_var],
    iso_host_metadata_sg[,selected_var]
  )
  poisson$num_VC[is.na(poisson$num_VC)] <- 0
  
  poisson %>% group_by(CateType,Presence) %>%
    do(poisson_test=poisson_test(.$num_VC)) -> poisson_tmp
  poisson_df <- cbind.data.frame(poisson_tmp[,1:2],
                   list.rbind(poisson_tmp$poisson_test))
  poisson_df$CateType <- 
    factor(poisson_df$CateType,
           levels = c("All","Cold shock protein","Acid shock protein",
                      "Alkaline shock protein","Phage shock protein",
                      "CRISPR-Cas system","RM system","Growth rate","OGT"
                      ))
  poisson_df %>% arrange(CateType) -> poisson_df
  poisson_df$`P value` <- "< 0.001"
  poisson_df <- poisson_df[,c(1:5,7)]
  colnames(poisson_df) <- c("Items","Prensence","N","Lambda","SD","P value")
  poisson_df$Lambda <- round(poisson_df$Lambda,3)
  poisson_df$SD <- round(poisson_df$SD,3)
  write.table(poisson_df,"../../05.Table/01.Raw/Poisson_test.txt",
              quote=F,row.names = F,sep="\t")
  
  poisson %>% group_by(CateType) %>% 
    filter(CateType != "All") %>%
    do(ks_test=ks_test_ldf(.$num_VC,.$Presence)) -> ks_test
  
  ttestsd(as.numeric(poisson_df[1,c(4,5,3)]),
          as.numeric(poisson_df[13,c(4,5,3)]))
  
  #iso_host_metadata_tmp <- iso_host_metadata 
  tmp_nvc <- iso_host_metadata$num_VC
  tmp_nvc[is.na(tmp_nvc)] <- 0
  poisson.lambda <- MASS::fitdistr(tmp_nvc,"Poisson")
  gf = goodfit(tmp_nvc,type= "poisson",method= "ML")
  gf.summary = capture.output(summary(gf))[[5]]
  pvalue = unlist(strsplit(gf.summary, split = " "))
  pvalue = as.numeric(pvalue[length(pvalue)])
  tmp_nvc <- count(tmp_nvc)
  tmp_nvc$Pro <- tmp_nvc$freq/sum(tmp_nvc$freq)
  tmp_nvc$Pro_pre <- dpois(tmp_nvc$x,poisson.lambda$estimate)
  tmp_nvc %>% 
    #filter(!is.na(host_vc_R2)) %>%
    ggplot(aes(x=x,y=Pro))+
    geom_col(color="black",fill="white")+
    mytheme+
    scale_x_continuous(expand = c(0.01,0,0,0),breaks = seq(0,12))+
    scale_y_continuous(expand = c(0,0,0,0.1))+
    geom_point(aes(x=x,y=Pro_pre),color="blue")+
    geom_text(x=7,y=0.25,label="Poisson distribution\nLambda=1.54, p = 0.017",size=6)+
    labs(x="Number of virus clusters",y="Proportion") -> p_vc_density
  ggsave("../../04.Figure/01.Raw/p_vc_no_crisprcas_density.tiff",
         plot=p_vc_density,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
}
#8. frequency of OGT for host with Csp####
iso_host_metadata_csp %>% 
  filter(`Cold_shock_protein|CspD` >=3) %>%
  ggplot(aes(x=OGT))+
  geom_histogram(stat = "bin",
                 position = "identity",
                 bins=8,fill="white",
                 colour="black")+
  scale_y_continuous(expand=c(0,0,0.1,0))+
  mytheme+
  labs(y="Number of genomes\n (with CspD ≥ 3)") -> hist_OGT_cspD

iso_host_metadata_csp %>% 
  filter(`Cold_shock_protein|CspC` >=3) %>%
  ggplot(aes(x=OGT))+
  geom_histogram(stat = "bin",
                 position = "identity",
                 bins=8,fill="white",
                 colour="black")+
  scale_y_continuous(expand=c(0,0,0.1,0))+
  mytheme+
  labs(y="Number of genomes\n (with CspC ≥ 3)") -> hist_OGT_cspC

cowplot::plot_grid(hist_OGT_cspC,hist_OGT_cspD,
                   labels = c('A', 'B')) -> p_hist_OGT_csp
ggsave("../../04.Figure/01.Raw/p_hist_OGT_csp.tiff",
       plot=p_hist_OGT_csp,device = "tiff",compression="lzw",
       width = 6.2*2,height = 4.8,units = "in")

#9. pc functions####
pc_host_metadata <- 
  iso_host_metadata_csp[,colnames(iso_host_metadata_csp) %in% 
  c("Assembly.Accession","OGT","DoublingTime","host_vc_R2")]
pc_host_csp <- 
  iso_host_metadata_csp.pos[,colnames(iso_host_metadata_csp.pos) %in%
  c("Assembly.Accession","num_csp")]
iso_vpc_final %>% group_by(host,pc) %>% dplyr::count() -> host_pc


kegg_ref <- read.table("../../02.Data/01.Processing/Maps/kegg_ref",
                       stringsAsFactors = F,sep="\t")
colnames(kegg_ref) <- c("KO","funcs")
kegg_ref_meta <- list.rbind(lapply(strsplit(kegg_ref$funcs,split="_"),
                  function(x){
                    tmp <- data.frame(M1=x[1],M2=x[2])
                    return(tmp)}))
kegg_ref <- cbind.data.frame(kegg_ref,kegg_ref_meta)
iso_pc_funcs_COG_KEGG <- 
  iso_pc_funcs[,c("pc","COG.Functional.cat.","KEGG_ko")]
iso_pc_funcs_COG_KEGG$KEGG_ko <- 
  str_extract(iso_pc_funcs_COG_KEGG$KEGG_ko,"K[0-9]+")

iso_pc_funcs_COG_KEGG <- 
  merge.data.frame(iso_pc_funcs_COG_KEGG,
                   kegg_ref,
                   by.x = "KEGG_ko",by.y = "KO",all.x = T)
iso_pc_funcs_COG_KEGG <- apply(iso_pc_funcs_COG_KEGG,2,function(x){
  tmp <- x
  if(class(x) == "factor") tmp <- as.character.factor(x)
  return(tmp)
})
iso_pc_funcs_COG_KEGG <- as.data.frame(iso_pc_funcs_COG_KEGG)
iso_pc_funcs_COG_KEGG$COG.Functional.cat. <-
  str_extract(iso_pc_funcs_COG_KEGG$COG.Functional.cat.,"[A-Z]")
iso_pc_funcs_COG_KEGG %>% 
  mutate_if(is.factor, as.character.factor) -> iso_pc_funcs_COG_KEGG
iso_pc_funcs_COG_KEGG[is.na(iso_pc_funcs_COG_KEGG)] <- "Unknown"


host_pc_meta <- merge.data.frame(host_pc,pc_host_metadata,
                 by.x="host",by.y = "Assembly.Accession",all.x=T)
host_pc_meta <- merge.data.frame(host_pc_meta,pc_host_csp,
                 by.x="host",by.y="Assembly.Accession",all.x=T)
host_pc_meta <- merge.data.frame(host_pc_meta,iso_pc_funcs_COG_KEGG,
                 by="pc",all.x=T)
host_pc_meta$num_csp <- ifelse(is.na(host_pc_meta$num_csp),
                               0,host_pc_meta$num_csp)
host_pc_meta$trophic <-
  ifelse(host_pc_meta$OGT >= 40, "Thermophilic",
         "Others")
host_pc_meta$growth <- 
  ifelse(host_pc_meta$DoublingTime>=2.5,"Slow","Fast")
host_pc_meta$num_csp <- ifelse(is.na(host_pc_meta$num_csp),
                 0,host_pc_meta$num_csp)
host_pc_meta$trophic2 <- 
  ifelse(host_pc_meta$num_csp!=0,"Csp","No csp")
#10.1 cog for OGT ####
host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
  group_by(COG.Functional.cat.,trophic) %>%
  do(data.frame(n=sum(.$n))) %>%
  group_by(trophic) -> host_pc_ogt_cog

host_pc_ogt_cog %>% 
#  filter(COG.Functional.cat. != "Unknown") %>%
  group_by(trophic) %>% 
  do(data.frame(S=sum(.$n))) -> cog_sum
cog <- unique(host_pc_ogt_cog$COG.Functional.cat.)
host_pc_ogt_cog %>% 
  filter(trophic == "Others" & 
           COG.Functional.cat. %in% cog[-21]) -> a

host_pc_ogt_cog %>% 
  filter(trophic != "Others" & 
           COG.Functional.cat. %in% cog[-21]) -> b
M <- as.table(rbind(a$n[-c(1,2)], b$n))
colnames(M) <- b$COG.Functional.cat.

cog2 <- cog[-c(1,2,21)]
chisq_res <- list()
for(i in 1:length(cog2))
{
  M2 <- as.table(rbind(M[,cog2[i]],cog_sum$S-M[,cog2[i]]))
  chisq_res[[i]] <- chisq.test(M2)$p.value
}
names(chisq_res) <- cog2
chisq_res <- ldply(chisq_res)
colnames(chisq_res) <- c("COG","p.value")
chisq_res <- rbind.data.frame(
  data.frame(
    COG=c("A","B"),
    p.value=1
  ),
  chisq_res
)
host_pc_ogt_cog %>% 
#  filter(COG.Functional.cat. != "Unknown") %>%
  group_by(trophic) %>% 
  do(data.frame(S=sum(.$n),
                COG.Functional.cat.=.$COG.Functional.cat.,
                trophic=.$trophic,
                n=.$n
  )) -> host_pc_ogt_cog2
host_pc_ogt_cog2$freq <- host_pc_ogt_cog2$n/host_pc_ogt_cog2$S
host_pc_ogt_cog2 %>%
  group_by(COG.Functional.cat.) %>%
  do(data.frame(Max=max(.$freq))) -> host_pc_ogt_cog2.maxfreq
chisq_res$y <- host_pc_ogt_cog2.maxfreq$Max[-21]
chisq_res$sig <- ""
chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
chisq_res$COG <- factor(chisq_res$COG,
                        levels = levels(host_pc_ogt_cog2.pl$COG.Functional.cat.))
chisq_res %>% filter(!COG %in% c("S")) -> chisq_res2

host_pc_ogt_cog2 %>%
  filter(!COG.Functional.cat. %in%
           c("S","Unknown")) -> host_pc_ogt_cog2.pl
host_pc_ogt_cog2.pl$trophic <- 
  factor(host_pc_ogt_cog2.pl$trophic,
         levels = c("Thermophilic","Others"))

  ggplot(host_pc_ogt_cog2.pl,
         aes(x=as.character.factor(COG.Functional.cat.),y=freq*100))+
  geom_col(aes(fill=trophic),position = "dodge",colour="black")+
  #geom_text(x="A",y=0.02,label="11")
  geom_text(data=chisq_res2,aes(x=COG,y=(y+0.002)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="COG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)) -> Fig.COG.ogt

  merge.data.frame(host_pc_ogt_cog2.pl,chisq_res2,
        by.x="COG.Functional.cat.",by.y="COG") -> Fig.COG.ogt.pl
  
#10.2 cog for growth ####
  host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
    group_by(COG.Functional.cat.,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(growth) -> host_pc_ogt_cog
  
  host_pc_ogt_cog %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> cog_sum
  cog <- unique(host_pc_ogt_cog$COG.Functional.cat.)
  host_pc_ogt_cog %>% 
    filter(growth == "Fast" & 
             COG.Functional.cat. %in% cog[-21]) -> a
  
  host_pc_ogt_cog %>% 
    filter(growth != "Fast" & 
             COG.Functional.cat. %in% cog[-21]) -> b
  M <- as.table(rbind(a$n[-c(1,2)], b$n))
  colnames(M) <- b$COG.Functional.cat.
  
  cog2 <- cog[-c(21)]
  chisq_res <- list()
  for(i in 1:length(cog2))
  {
    M2 <- as.table(rbind(M[,cog2[i]],cog_sum$S-M[,cog2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- cog2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("COG","p.value")
  
  host_pc_ogt_cog %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  COG.Functional.cat.=.$COG.Functional.cat.,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_ogt_cog2
  host_pc_ogt_cog2$freq <- host_pc_ogt_cog2$n/host_pc_ogt_cog2$S
  host_pc_ogt_cog2 %>%
    group_by(COG.Functional.cat.) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_cog2.maxfreq
  chisq_res$y <- host_pc_ogt_cog2.maxfreq$Max[-c(21)]
  chisq_res$sig <- ""
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  chisq_res$COG <- factor(chisq_res$COG,
                          levels = levels(host_pc_ogt_cog2.pl$COG.Functional.cat.))
  chisq_res %>% filter(!COG %in% c("S")) -> chisq_res2
  host_pc_ogt_cog2 %>%
    filter(!COG.Functional.cat. %in%
             c("S","Unknown")) -> host_pc_ogt_cog2.pl
  host_pc_ogt_cog2.pl$growth <- 
    factor(host_pc_ogt_cog2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_ogt_cog2.pl,
         aes(x=as.character.factor(COG.Functional.cat.),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(data=chisq_res2,aes(x=COG,y=(y+0.002)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="COG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("yellow","white"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)) -> Fig.COG.growth
  merge.data.frame(host_pc_ogt_cog2.pl,chisq_res2,
       by.x="COG.Functional.cat.",by.y="COG") -> Fig.COG.growth.pl
#10.3 cog for psychrophilic####
  host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
    group_by(COG.Functional.cat.,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_ogt_cog
  
  host_pc_ogt_cog %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> cog_sum
  cog <- unique(host_pc_ogt_cog$COG.Functional.cat.)
  host_pc_ogt_cog %>% 
    filter(trophic2 == "Csp") -> a
  
  host_pc_ogt_cog %>% 
    filter(trophic2 != "Csp" ) -> b
  M <- as.table(rbind(a$n[a$COG.Functional.cat. %in% b$COG.Functional.cat. ],
                      b$n[b$COG.Functional.cat. %in% a$COG.Functional.cat.]))
    colnames(M) <- b$COG.Functional.cat.[b$COG.Functional.cat. %in% a$COG.Functional.cat.]
  
  cog2 <- cog[-c(21)]
  chisq_res <- list()
  for(i in 1:length(cog2))
  {
    M2 <- as.table(rbind(M[,cog2[i]],cog_sum$S-M[,cog2[i]]))
    chisq_res[[i]] <- chisq.test(M2,B=10000)$p.value
  }
  names(chisq_res) <- cog2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("COG","p.value")
  
  host_pc_ogt_cog %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  COG.Functional.cat.=.$COG.Functional.cat.,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_ogt_cog2
  host_pc_ogt_cog2$freq <- host_pc_ogt_cog2$n/host_pc_ogt_cog2$S
  host_pc_ogt_cog2 %>%
    group_by(COG.Functional.cat.) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_cog2.maxfreq
  host_pc_ogt_cog2.maxfreq <- as.data.frame(host_pc_ogt_cog2.maxfreq)
  rownames( host_pc_ogt_cog2.maxfreq) <-  host_pc_ogt_cog2.maxfreq$COG.Functional.cat.
  chisq_res$y <- host_pc_ogt_cog2.maxfreq[chisq_res$COG,"Max"]
  chisq_res$sig <- ""
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  chisq_res$COG <- factor(chisq_res$COG,
                          levels = levels(host_pc_ogt_cog2.pl$COG.Functional.cat.))
  chisq_res %>% filter(!COG %in% c("S")) -> chisq_res2
  host_pc_ogt_cog2 %>%
    filter(!COG.Functional.cat. %in%
             c("S","Unknown")) -> host_pc_ogt_cog2.pl
  host_pc_ogt_cog2.pl$trophic2 <- 
    factor(host_pc_ogt_cog2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_ogt_cog2.pl,
         aes(x=as.character.factor(COG.Functional.cat.),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(data=chisq_res2,aes(x=COG,y=(y+0.002)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="COG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Presence","Absence"))+
    theme(legend.position = c(0.85,0.8)) -> Fig.COG.csp
  
  merge.data.frame(host_pc_ogt_cog2.pl,chisq_res2,
                   by.x="COG.Functional.cat.",by.y="COG") -> Fig.COG.csp.pl
#10.4 COG all####
Fig.COG.ogt.pl$type="OGT"
Fig.COG.growth.pl$type="Growth rate"
Fig.COG.csp.pl$type="CSP"
colnames(Fig.COG.ogt.pl)[3] <- "Var"
colnames(Fig.COG.growth.pl)[3] <- "Var"
colnames(Fig.COG.csp.pl)[3] <- "Var"
Fig.COG.pl <- rbind.data.frame(Fig.COG.ogt.pl,Fig.COG.growth.pl)
Fig.COG.pl <- rbind.data.frame(Fig.COG.pl,Fig.COG.csp.pl)
Fig.COG.pl$type <- 
  factor(Fig.COG.pl$type,levels = c("OGT","Growth rate","CSP"))
Fig.COG.pl %>% ggplot(aes(x=COG.Functional.cat.,y=freq*100))+
  geom_col(aes(fill=Var),position = "dodge",colour="black")+
  #geom_text(x="A",y=0.02,label="11")
  geom_text(aes(x=COG.Functional.cat.,y=(y+0.002)*100,label=sig))+
  mytheme +
  scale_y_continuous(expand = c(0,0,0,1)) + 
  facet_wrap(~type,nrow=3)+
  labs(x="COG categories",y="Proportion %")+
  scale_fill_manual(values = c("red","grey","yellow",
                               "white","green","black"))+
  guides(fill=F) -> Fig.COG.all

ggsave("../../04.Figure/01.Raw/Fig.COG.all.tiff",
       plot= Fig.COG.all,device = "tiff",compression="lzw",
       width = 6.2,height = 3.5*2,units = "in")
#11.1 kegg for OGT ####
  host_pc_meta$M2 <- as.character.factor(host_pc_meta$M2)
  host_pc_meta$M2 <-
    ifelse(is.na(host_pc_meta$M2),"Unknown",  host_pc_meta$M2)
  
  host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
    group_by(M2,trophic) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic) -> host_pc_ogt_kegg
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n))) -> kegg_sum
  kegg <- unique(host_pc_ogt_kegg$M2)
  host_pc_ogt_kegg %>% 
    filter(trophic == "Others" & 
             M2 %in% kegg) -> a
  
  host_pc_ogt_kegg %>% 
    filter(trophic != "Others" & 
             M2 %in% kegg) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$M2 %in% b$M2], as.numeric(b$n)))
  colnames(M) <- b$M2
  
  kegg2 <- kegg[a$M2 %in% b$M2]
  chisq_res <- list()
  for(i in 1:length(kegg2))
  {
    M2 <- as.table(rbind(M[,kegg2[i]],kegg_sum$S-M[,kegg2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- kegg2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("M2","p.value")
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n),
                  M2=.$M2,
                  trophic=.$trophic,
                  n=.$n
    )) -> host_pc_ogt_kegg2
  host_pc_ogt_kegg2$freq <- host_pc_ogt_kegg2$n/host_pc_ogt_kegg2$S
  rownames(host_pc_ogt_kegg2.maxfreq) <- host_pc_ogt_kegg2.maxfreq$M2
  host_pc_ogt_kegg2 %>%
    group_by(M2) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_kegg2.maxfreq
  host_pc_ogt_kegg2.maxfreq <-
    as.data.frame(host_pc_ogt_kegg2.maxfreq)
  rownames(host_pc_ogt_kegg2.maxfreq) <-  host_pc_ogt_kegg2.maxfreq$M2
  chisq_res$y <- host_pc_ogt_kegg2.maxfreq[chisq_res$M2,"Max"]
  chisq_res$sig <- ""
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  host_pc_ogt_kegg2 %>%
    filter(!M2 %in%
             c("Overview","Unknown") &
             M2 %in% b$M2
             ) -> host_pc_ogt_kegg2.pl
  
  chisq_res$M2 <- factor(chisq_res$M2,
                         levels = levels(host_pc_ogt_kegg2.pl$M2))
  
  host_pc_ogt_kegg2.pl$trophic <- 
    factor(host_pc_ogt_kegg2.pl$trophic,
           levels = c("Thermophilic","Others"))
  
  ggplot(host_pc_ogt_kegg2.pl,
         aes(x=as.character.factor(M2),y=freq*100))+
    geom_col(aes(fill=trophic),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=M2,y=(y+0.001)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)
          ) -> Fig.KEGG.ogt
  merge.data.frame(host_pc_ogt_kegg2.pl,chisq_res2,
                   by.x="M2",by.y="M2") -> Fig.KEGG.ogt.pl  
#11.2 kegg for growth ####
  host_pc_meta$M2 <- as.character.factor(host_pc_meta$M2)
  host_pc_meta$M2 <-
    ifelse(is.na(host_pc_meta$M2),"Unknown",  host_pc_meta$M2)
  
  host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
    group_by(M2,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(growth) -> host_pc_ogt_kegg
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> kegg_sum
  kegg <- unique(host_pc_ogt_kegg$M2)
  host_pc_ogt_kegg %>% 
    filter(growth == "Fast" & 
             M2 %in% kegg) -> a
  
  host_pc_ogt_kegg %>% 
    filter(growth != "Fast" & 
             M2 %in% kegg) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$M2 %in% b$M2], 
                      as.numeric(b$n)[b$M2 %in% a$M2] ))
  colnames(M) <- b$M2[b$M2 %in% a$M2]
  
  kegg2 <- unique(b$M2[b$M2 %in% a$M2])
  chisq_res <- list()
  for(i in 1:length(kegg2))
  {
    M2 <- as.table(rbind(M[,kegg2[i]],kegg_sum$S-M[,kegg2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- kegg2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("M2","p.value")
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  M2=.$M2,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_ogt_kegg2
  host_pc_ogt_kegg2$freq <- host_pc_ogt_kegg2$n/host_pc_ogt_kegg2$S
  host_pc_ogt_kegg2 %>%
    group_by(M2) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_kegg2.maxfreq
  host_pc_ogt_kegg2.maxfreq <- 
    as.data.frame(host_pc_ogt_kegg2.maxfreq)
  rownames(host_pc_ogt_kegg2.maxfreq) <-host_pc_ogt_kegg2.maxfreq$M2
  chisq_res$y <- host_pc_ogt_kegg2.maxfreq[chisq_res$M2,"Max"]
  chisq_res$sig <- ""
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  host_pc_ogt_kegg2 %>%
    filter(!M2 %in%
             c("Overview","Unknown") &
             M2 %in% b$M2
    ) -> host_pc_ogt_kegg2.pl
  
  chisq_res$M2 <- factor(chisq_res$M2,
                         levels = levels(host_pc_ogt_kegg2.pl$M2))
  
  host_pc_ogt_kegg2.pl$growth <- 
    factor(host_pc_ogt_kegg2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_ogt_kegg2.pl,
         aes(x=as.character.factor(M2),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=M2,y=(y+0.001)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("yellow","white"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)
    ) -> Fig.KEGG.growth

  merge.data.frame(host_pc_ogt_kegg2.pl,chisq_res2,
                   by.x="M2",by.y="M2") -> Fig.KEGG.growth.pl 
#11.3 kegg for psychrophilic ####
  host_pc_meta$M2 <- as.character.factor(host_pc_meta$M2)
  host_pc_meta$M2 <-
    ifelse(is.na(host_pc_meta$M2),"Unknown",  host_pc_meta$M2)
  
  host_pc_meta %>% filter(!is.na(host_vc_R2)) %>%
    group_by(M2,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_ogt_kegg
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> kegg_sum
  kegg <- unique(host_pc_ogt_kegg$M2)
  host_pc_ogt_kegg %>% 
    filter(trophic2 == "Csp") -> a
  
  host_pc_ogt_kegg %>% 
    filter(trophic2 != "Csp") -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$M2 %in% b$M2], 
                      as.numeric(b$n)[b$M2 %in% a$M2] ))
  colnames(M) <- b$M2[b$M2 %in% a$M2]
  
  kegg2 <- unique(b$M2[b$M2 %in% a$M2])
  chisq_res <- list()
  for(i in 1:length(kegg2))
  {
    M2 <- as.table(rbind(M[,kegg2[i]],kegg_sum$S-M[,kegg2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- kegg2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("M2","p.value")
  
  host_pc_ogt_kegg %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  M2=.$M2,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_ogt_kegg2
  host_pc_ogt_kegg2$freq <- host_pc_ogt_kegg2$n/host_pc_ogt_kegg2$S
  host_pc_ogt_kegg2 %>%
    group_by(M2) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_kegg2.maxfreq
  host_pc_ogt_kegg2.maxfreq <- 
    as.data.frame(host_pc_ogt_kegg2.maxfreq)
  rownames(host_pc_ogt_kegg2.maxfreq) <-host_pc_ogt_kegg2.maxfreq$M2
  chisq_res$y <- host_pc_ogt_kegg2.maxfreq[chisq_res$M2,"Max"]
  chisq_res$sig <- ""
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  host_pc_ogt_kegg2 %>%
    filter(!M2 %in%
             c("Overview","Unknown") &
             M2 %in% b$M2
    ) -> host_pc_ogt_kegg2.pl
  
  chisq_res$M2 <- factor(chisq_res$M2,
                         levels = levels(host_pc_ogt_kegg2.pl$M2))
  
  host_pc_ogt_kegg2.pl$trophic2 <- 
    factor(host_pc_ogt_kegg2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_ogt_kegg2.pl,
         aes(x=as.character.factor(M2),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=M2,y=(y+0.0008)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Presence","Absence"))+
    theme(legend.position = c(0.85,0.8)
    ) -> Fig.KEGG.csp

  merge.data.frame(host_pc_ogt_kegg2.pl,chisq_res2,
                   by.x="M2",by.y="M2") -> Fig.KEGG.csp.pl
  
#11.4 KEGG all####
  Fig.KEGG.ogt.pl$type="OGT"
  Fig.KEGG.growth.pl$type="Growth rate"
  Fig.KEGG.csp.pl$type="CSP"
  colnames(Fig.KEGG.ogt.pl)[3] <- "Var"
  colnames(Fig.KEGG.growth.pl)[3] <- "Var"
  colnames(Fig.KEGG.csp.pl)[3] <- "Var"
  Fig.KEGG.pl <- rbind.data.frame(Fig.KEGG.ogt.pl,Fig.KEGG.growth.pl)
  Fig.KEGG.pl <- rbind.data.frame(Fig.KEGG.pl,Fig.KEGG.csp.pl)
  Fig.KEGG.pl$type <- 
    factor(Fig.KEGG.pl$type,levels = c("OGT","Growth rate","CSP"))
  Fig.KEGG.pl %>% 
    filter(freq >=0.0005 &
             !M2 %in% c("Biosynthesis of other secondary metabolites",
                        "Cancers","Drug resistance",
                        "Metabolism of terpenoids and polyketides",
                        "Infectious diseases"
                        )
             ) %>% 
    ggplot(aes(x=M2,y=freq*100))+
    geom_col(aes(fill=Var),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(aes(x=M2,y=(y+0.008)*100,label=sig),angle=90)+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,1)) + 
    coord_flip()+
    facet_wrap(~type,nrow=1)+
    labs(x="KEGG functional categories",y="Proportion %")+
    scale_fill_manual(values = c("red","grey","yellow",
                                 "white","green","black"))+
    guides(fill=F) +
    theme(axis.text = element_text(size=10)#,
          #axis.text.x=element_text(angle=65,
          #                         hjust=1,vjust=1)
          )-> Fig.KEGG.all
  
  ggsave("../../04.Figure/01.Raw/Fig.KEGG.all.tiff",
         plot= Fig.KEGG.all,device = "tiff",compression="lzw",
         width = 6.6,height = 3.5*1.3,units = "in")

#12.1 NCycDB####
  iso_pc_funcs_NCyc <- 
    iso_pc_funcs[,c("pc","Ncyc_ID","Ncyc_gene")]
  host_pc_metaN <- merge.data.frame(host_pc,pc_host_metadata,
                                   by.x="host",by.y = "Assembly.Accession",all.x=T)
  host_pc_metaN <- merge.data.frame(host_pc_metaN,pc_host_csp,
                                   by.x="host",by.y="Assembly.Accession",all.x=T)
  host_pc_metaN <- merge.data.frame(host_pc_metaN,iso_pc_funcs_NCyc,
                                   by="pc",all.x=T)
  host_pc_metaN$num_csp <- ifelse(is.na(host_pc_metaN$num_csp),
                                 0,host_pc_metaN$num_csp)
  host_pc_metaN$trophic <-
    ifelse(host_pc_metaN$OGT >= 40, "Thermophilic",
           "Others")
  host_pc_metaN$growth <- 
    ifelse(host_pc_metaN$DoublingTime>=2.5,"Slow","Fast")
  host_pc_metaN$num_csp <- ifelse(is.na(host_pc_metaN$num_csp),
                                 0,host_pc_metaN$num_csp)
  host_pc_metaN$trophic2 <- 
    ifelse(host_pc_metaN$num_csp!=0,"Csp","No csp")
#12.2 NCyc for OGT ####
  host_pc_metaN %>% filter(!is.na(Ncyc_ID)) %>% dim()
  
  host_pc_metaN$Ncyc_gene <- 
    as.character.factor(host_pc_metaN$Ncyc_gene)
  host_pc_metaN$Ncyc_gene <-
    ifelse(is.na(host_pc_metaN$Ncyc_gene),"Unknown", 
           host_pc_metaN$Ncyc_gene)
  
  host_pc_metaN %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Ncyc_gene,trophic) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic) -> host_pc_ogt_ncyc
  
  host_pc_ogt_ncyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n))) -> ncyc_sum
  ncyc <- unique(host_pc_ogt_ncyc$Ncyc_gene)
  host_pc_ogt_ncyc %>% 
    filter(trophic == "Others" & 
             Ncyc_gene %in% ncyc) -> a
  
  host_pc_ogt_ncyc %>% 
    filter(trophic != "Others" & 
             Ncyc_gene %in% ncyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Ncyc_gene %in% b$Ncyc_gene],
                      as.numeric(b$n)))
  colnames(M) <- b$Ncyc_gene
  
  ncyc2 <- ncyc[a$Ncyc_gene %in% b$Ncyc_gene]
  chisq_res <- list()
  for(i in 1:length(ncyc2))
  {
    M2 <- as.table(rbind(M[,ncyc2[i]],ncyc_sum$S-M[,ncyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- ncyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Ncyc_gene","p.value")
  
  host_pc_ogt_ncyc %>% 
 #     filter(Ncyc_gene != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n),
                  Ncyc_gene=.$Ncyc_gene,
                  trophic=.$trophic,
                  n=.$n
    )) -> host_pc_ogt_ncyc2
  host_pc_ogt_ncyc2$freq <- host_pc_ogt_ncyc2$n/host_pc_ogt_ncyc2$S
  
  host_pc_ogt_ncyc2 %>%
    group_by(Ncyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_ncyc2.maxfreq
  host_pc_ogt_ncyc2.maxfreq <-
    as.data.frame(host_pc_ogt_ncyc2.maxfreq)
  rownames(host_pc_ogt_ncyc2.maxfreq) <- 
    host_pc_ogt_ncyc2.maxfreq$Ncyc_gene
  chisq_res$y <- host_pc_ogt_ncyc2.maxfreq[chisq_res$Ncyc_gene,"Max"]
  chisq_res$sig <- ""
#  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Ncyc_gene != "Unknown")
  host_pc_ogt_ncyc2 %>%
    filter( Ncyc_gene != "Unknown"
           # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_ncyc2.pl
  
  chisq_res$Ncyc_gene <- factor(chisq_res$Ncyc_gene,
                         levels = levels(host_pc_ogt_ncyc2.pl$Ncyc_gene))
  
  host_pc_ogt_ncyc2.pl$trophic <- 
    factor(host_pc_ogt_ncyc2.pl$trophic,
           levels = c("Thermophilic","Others"))
  
  ggplot(host_pc_ogt_ncyc2.pl,
         aes(x=as.character.factor(Ncyc_gene),y=freq*100))+
    geom_col(aes(fill=trophic),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Ncyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)#,
#          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.ncyc.ogt
  
  merge.data.frame(host_pc_ogt_ncyc2.pl,chisq_res2,
                   by.x="Ncyc_gene",
                   by.y="Ncyc_gene",all.x=T) -> Fig.ncyc.ogt.pl  

#12.3 NCyc for growth ####
  #host_pc_metaN$Ncyc_gene <- 
  #  as.character.factor(host_pc_metaN$Ncyc_gene)
  host_pc_metaN$Ncyc_gene <-
    ifelse(is.na(host_pc_metaN$Ncyc_gene),"Unknown", 
           host_pc_metaN$Ncyc_gene)
  
  host_pc_metaN %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Ncyc_gene,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(growth) -> host_pc_ogt_ncyc
  
  host_pc_ogt_ncyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> ncyc_sum
  ncyc <- unique(host_pc_ogt_ncyc$Ncyc_gene)
  host_pc_ogt_ncyc %>% 
    filter(growth == "Fast" & 
             Ncyc_gene %in% ncyc) -> a
  
  host_pc_ogt_ncyc %>% 
    filter(growth != "Fast" & 
             Ncyc_gene %in% ncyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Ncyc_gene %in% b$Ncyc_gene],
                      as.numeric(b$n)[b$Ncyc_gene %in% a$Ncyc_gene] ))
  colnames(M) <- b$Ncyc_gene[b$Ncyc_gene %in% a$Ncyc_gene]
  
  ncyc2 <- b$Ncyc_gene[b$Ncyc_gene %in% a$Ncyc_gene]
  chisq_res <- list()
  for(i in 1:length(ncyc2))
  {
    M2 <- as.table(rbind(M[,ncyc2[i]],ncyc_sum$S-M[,ncyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- ncyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Ncyc_gene","p.value")
  
  host_pc_ogt_ncyc %>% 
  #  filter(Ncyc_gene != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  Ncyc_gene=.$Ncyc_gene,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_ogt_ncyc2
  host_pc_ogt_ncyc2$freq <- host_pc_ogt_ncyc2$n/host_pc_ogt_ncyc2$S
  
  host_pc_ogt_ncyc2 %>%
    group_by(Ncyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_ncyc2.maxfreq
  host_pc_ogt_ncyc2.maxfreq <-
    as.data.frame(host_pc_ogt_ncyc2.maxfreq)
  rownames(host_pc_ogt_ncyc2.maxfreq) <- 
    host_pc_ogt_ncyc2.maxfreq$Ncyc_gene
  chisq_res$y <- host_pc_ogt_ncyc2.maxfreq[chisq_res$Ncyc_gene,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Ncyc_gene != "Unknown")
  host_pc_ogt_ncyc2 %>%
    filter( Ncyc_gene != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_ncyc2.pl
  
  chisq_res$Ncyc_gene <- factor(chisq_res$Ncyc_gene,
                                levels = levels(host_pc_ogt_ncyc2.pl$Ncyc_gene))
  
  host_pc_ogt_ncyc2.pl$gr <- 
    factor(host_pc_ogt_ncyc2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_ogt_ncyc2.pl,
         aes(x=as.character.factor(Ncyc_gene),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Ncyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("yellow","white"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.ncyc.growth
  
  merge.data.frame(host_pc_ogt_ncyc2.pl,chisq_res2,
                   by.x="Ncyc_gene",
                   by.y="Ncyc_gene",all.x=T) -> Fig.ncyc.growth.pl  

#12.4 NCyc for psychrophilic ####
 # host_pc_metaN$Ncyc_gene <- 
 #   as.character.factor(host_pc_metaN$Ncyc_gene)
  host_pc_metaN$Ncyc_gene <-
    ifelse(is.na(host_pc_metaN$Ncyc_gene),"Unknown", 
           host_pc_metaN$Ncyc_gene)
  
  host_pc_metaN %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Ncyc_gene,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_ogt_ncyc
  
  host_pc_ogt_ncyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> ncyc_sum
  ncyc <- unique(host_pc_ogt_ncyc$Ncyc_gene)
  host_pc_ogt_ncyc %>% 
    filter(trophic2 == "Csp" & 
             Ncyc_gene %in% ncyc) -> a
  
  host_pc_ogt_ncyc %>% 
    filter(trophic2 != "Csp" & 
             Ncyc_gene %in% ncyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Ncyc_gene %in% b$Ncyc_gene],
                      as.numeric(b$n)[b$Ncyc_gene %in% a$Ncyc_gene] ))
  colnames(M) <- b$Ncyc_gene[b$Ncyc_gene %in% a$Ncyc_gene]
  
  ncyc2 <- b$Ncyc_gene[b$Ncyc_gene %in% a$Ncyc_gene]
  chisq_res <- list()
  for(i in 1:length(ncyc2))
  {
    M2 <- as.table(rbind(M[,ncyc2[i]],ncyc_sum$S-M[,ncyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- ncyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Ncyc_gene","p.value")
  
  host_pc_ogt_ncyc %>% 
  #  filter(Ncyc_gene != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  Ncyc_gene=.$Ncyc_gene,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_ogt_ncyc2
  host_pc_ogt_ncyc2$freq <- host_pc_ogt_ncyc2$n/host_pc_ogt_ncyc2$S
  
  host_pc_ogt_ncyc2 %>%
    group_by(Ncyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_ncyc2.maxfreq
  host_pc_ogt_ncyc2.maxfreq <-
    as.data.frame(host_pc_ogt_ncyc2.maxfreq)
  rownames(host_pc_ogt_ncyc2.maxfreq) <- 
    host_pc_ogt_ncyc2.maxfreq$Ncyc_gene
  chisq_res$y <- host_pc_ogt_ncyc2.maxfreq[chisq_res$Ncyc_gene,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Ncyc_gene != "Unknown")
  host_pc_ogt_ncyc2 %>%
    filter( Ncyc_gene != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_ncyc2.pl
  
  chisq_res$Ncyc_gene <- factor(chisq_res$Ncyc_gene,
                                levels = levels(host_pc_ogt_ncyc2.pl$Ncyc_gene))
  
  host_pc_ogt_ncyc2.pl$gr <- 
    factor(host_pc_ogt_ncyc2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_ogt_ncyc2.pl,
         aes(x=as.character.factor(Ncyc_gene),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Ncyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Presence","Absence"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.ncyc.csp
  
  merge.data.frame(host_pc_ogt_ncyc2.pl,chisq_res2,
                   by.x="Ncyc_gene",
                   by.y="Ncyc_gene",all.x=T) -> Fig.ncyc.csp.pl 

#12.5 NCyc all####
  Fig.ncyc.ogt.pl$type="OGT"
  Fig.ncyc.growth.pl$type="Growth rate"
  Fig.ncyc.csp.pl$type="CSP"
  colnames(Fig.ncyc.ogt.pl)[3] <- "Var"
  colnames(Fig.ncyc.growth.pl)[3] <- "Var"
  colnames(Fig.ncyc.csp.pl)[3] <- "Var"
  Fig.ncyc.pl <- rbind.data.frame(Fig.ncyc.ogt.pl,
                                  Fig.ncyc.growth.pl[,-6])
  Fig.ncyc.pl <- rbind.data.frame(Fig.ncyc.pl,
                                  Fig.ncyc.csp.pl[,-6])
  Fig.ncyc.pl$type <- 
    factor(Fig.ncyc.pl$type,levels = c("OGT","Growth rate","CSP"))
  
  Fig.ncyc.pl %>% 
    filter(type=="OGT") %>%
    group_by(Ncyc_gene) %>% 
    dplyr::summarise(gs=mean(freq) ) -> na
  na %>% filter(gs >=0.0001) -> nb
  Fig.ncyc.pl$Ncyc_gene <- factor(Fig.ncyc.pl$Ncyc_gene,
    levels=levels(Fig.ncyc.pl$Ncyc_gene)[
    levels(Fig.ncyc.pl$Ncyc_gene) %in%
    unique(Fig.ncyc.pl$Ncyc_gene)])
  ng1 <- str_extract(levels(Fig.ncyc.pl$Ncyc_gene),"[a-z]+")
  ng2 <- str_remove(levels(Fig.ncyc.pl$Ncyc_gene),"[a-z]+")
  ng3 <- NULL
  ng3.2 <- NULL
  for(i in 1:length(ng1))
  {
    
    ng3 <- c(ng3,
             substitute(expression(italic(g1)*g2),
                        list(g1=ng1[i],g2=ng2[i]))
    )
  }
  
  Fig.ncyc.pl %>% 
    #filter(Ncyc_gene %in% nb$Ncyc_gene) %>%
    ggplot(aes(x=Ncyc_gene,y=freq*100))+
    geom_col(aes(fill=Var),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(aes(x=Ncyc_gene,y=(y+0.0001)*100,label=sig),
              angle=90)+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.03)) + 
    coord_flip()+
    facet_wrap(~type,nrow=1)+
    labs(x="Nitrogen cycle genes",y="Proportion %")+
    scale_fill_manual(values = c("red","grey","yellow",
                                 "white","green","black"))+
    guides(fill=F) +
    #scale_x_discrete(labels=ng3)+
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(face="italic"))-> Fig.Ncyc.all
  
  ggsave("../../04.Figure/01.Raw/Fig.Ncyc.all.tiff",
         plot= Fig.Ncyc.all,device = "tiff",compression="lzw",
         width = 6.2,height = 3.5*2.5,units = "in")
  ggsave("../../04.Figure/01.Raw/Fig.Ncyc.all.pdf",
         plot= Fig.Ncyc.all,
         width = 6.2,height = 3.5*2.5)
#13.1 BacMet####
  iso_pc_funcs_BacMet <- 
    iso_pc_funcs[,c("pc","BacMet_ID","Compound")]
  host_pc_metaBacMet <- merge.data.frame(host_pc,pc_host_metadata,
                                    by.x="host",by.y = "Assembly.Accession",all.x=T)
  host_pc_metaBacMet <- merge.data.frame(host_pc_metaBacMet,pc_host_csp,
                                    by.x="host",by.y="Assembly.Accession",all.x=T)
  host_pc_metaBacMet <- merge.data.frame(host_pc_metaBacMet,iso_pc_funcs_BacMet,
                                    by="pc",all.x=T)
  host_pc_metaBacMet$num_csp <- ifelse(is.na(host_pc_metaBacMet$num_csp),
                                  0,host_pc_metaBacMet$num_csp)
  host_pc_metaBacMet$trophic <-
    ifelse(host_pc_metaBacMet$OGT >= 40, "Thermophilic",
           "Others")
  host_pc_metaBacMet$growth <- 
    ifelse(host_pc_metaBacMet$DoublingTime>=2.5,"Slow","Fast")
  host_pc_metaBacMet$num_csp <- ifelse(is.na(host_pc_metaBacMet$num_csp),
                                  0,host_pc_metaBacMet$num_csp)
  host_pc_metaBacMet$trophic2 <- 
    ifelse(host_pc_metaBacMet$num_csp!=0,"Csp","No csp")
  host_pc_metaBacMet$Metal <- 
    str_extract(host_pc_metaBacMet$Compound,
                     "\\([A-Za-z]+\\)")
  Metal <- str_remove(
    str_remove(host_pc_metaBacMet$Metal,"[(]"),
    "[)]")
  NonMetal <- c("HCl","SDS","CTM","SDC","BAC","TPA","TPP",
                "CCCP","DAPI","CTAB","CPC","TBT")
  Metal[Metal %in% NonMetal] <- NA
  host_pc_metaBacMet$Metal <- Metal
  antibio <- str_remove(str_extract(
    host_pc_metaBacMet$Compound,
  "class: [\\w\\s]+"),"class: ")
  host_pc_metaBacMet$antibio <- antibio
#13.2 Metal for OGT####
  host_pc_metaBacMet$Metal <- 
    as.character.factor(host_pc_metaBacMet$Metal)
  host_pc_metaBacMet$Metal <-
    ifelse(is.na(host_pc_metaBacMet$Metal),"Unknown", 
           host_pc_metaBacMet$Metal)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Metal,trophic) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic) -> host_pc_ogt_metal
  
  host_pc_ogt_metal %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n))) -> metal_sum
  metal <- unique(host_pc_ogt_metal$Metal)
  host_pc_ogt_metal %>% 
    filter(trophic == "Others" & 
             Metal %in% metal) -> a
  
  host_pc_ogt_metal %>% 
    filter(trophic != "Others" & 
             Metal %in% metal) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Metal %in%
                                        b$Metal],
                      as.numeric(b$n)[b$Metal %in%
                                        a$Metal] ))
  colnames(M) <- b$Metal[b$Metal %in% a$Metal]
  
  metal2 <- metal[a$Metal %in% b$Metal]
  chisq_res <- list()
  for(i in 1:length(metal2))
  {
    M2 <- as.table(rbind(M[,metal2[i]],
                         metal_sum$S-M[,metal2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- metal2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Metal","p.value")
  
  host_pc_ogt_metal %>% 
   # filter(Metal != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n),
                  Metal=.$Metal,
                  trophic=.$trophic,
                  n=.$n
    )) -> host_pc_ogt_metal2
  host_pc_ogt_metal2$freq <-
    host_pc_ogt_metal2$n/host_pc_ogt_metal2$S
  
  host_pc_ogt_metal2 %>%
    group_by(Metal) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_ogt_metal2.maxfreq
  host_pc_ogt_metal2.maxfreq <-
    as.data.frame(host_pc_ogt_metal2.maxfreq)
  rownames(host_pc_ogt_metal2.maxfreq) <- 
    host_pc_ogt_metal2.maxfreq$Metal
  chisq_res$y <- host_pc_ogt_metal2.maxfreq[chisq_res$Metal,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Metal != "Unknown")
  host_pc_ogt_metal2 %>%
    filter( Metal != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_metal2.pl
  
  chisq_res$Metal <- factor(chisq_res$Metal,
        levels = levels(host_pc_ogt_metal2.pl$Metal))
  
  host_pc_ogt_metal2.pl$trophic <- 
    factor(host_pc_ogt_metal2.pl$trophic,
           levels = c("Thermophilic","Others"))
  
  ggplot(host_pc_ogt_metal2.pl,
         aes(x=as.character.factor(Metal),y=freq*100))+
    geom_col(aes(fill=trophic),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Metal,
                                  y=(y+0.001)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.5)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.metal.ogt
  
  merge.data.frame(host_pc_ogt_metal2.pl,chisq_res2,
                   by.x="Metal",
                   by.y="Metal",all.x=T) -> Fig.metal.ogt.pl  

#13.3 Metal for growth####
  host_pc_metaBacMet$Metal <-
    ifelse(is.na(host_pc_metaBacMet$Metal),"Unknown", 
           host_pc_metaBacMet$Metal)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Metal,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(growth) -> host_pc_growth_metal
  
  host_pc_growth_metal %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> metal_sum
  metal <- unique(host_pc_growth_metal$Metal)
  host_pc_growth_metal %>% 
    filter(growth == "Fast" & 
             Metal %in% metal) -> a
  
  host_pc_growth_metal %>% 
    filter(growth != "Fast" & 
             Metal %in% metal) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Metal %in%
                                        b$Metal],
                      as.numeric(b$n)))
  colnames(M) <- b$Metal[b$Metal %in% a$Metal]
  
  metal2 <- metal[a$Metal %in% b$Metal]
  chisq_res <- list()
  for(i in 1:length(metal2))
  {
    M2 <- as.table(rbind(M[,metal2[i]],
                         metal_sum$S-M[,metal2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- metal2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Metal","p.value")
  
  host_pc_growth_metal %>% 
    #filter(Metal != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  Metal=.$Metal,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_growth_metal2
  host_pc_growth_metal2$freq <-
    host_pc_growth_metal2$n/host_pc_growth_metal2$S
  
  host_pc_growth_metal2 %>%
    group_by(Metal) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_growth_metal2.maxfreq
  host_pc_growth_metal2.maxfreq <-
    as.data.frame(host_pc_growth_metal2.maxfreq)
  rownames(host_pc_growth_metal2.maxfreq) <- 
    host_pc_growth_metal2.maxfreq$Metal
  chisq_res$y <- host_pc_growth_metal2.maxfreq[
    chisq_res$Metal,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Metal != "Unknown")
  host_pc_growth_metal2 %>%
    filter( Metal != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_growth_metal2.pl
  
  chisq_res$Metal <- factor(chisq_res$Metal,
      levels = levels(host_pc_growth_metal2.pl$Metal))
  
  host_pc_growth_metal2.pl$growth <- 
    factor(host_pc_growth_metal2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_growth_metal2.pl,
         aes(x=as.character.factor(Metal),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",
             colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Metal,
                                  y=(y+0.0005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.metal.growth
  
  merge.data.frame(host_pc_growth_metal2.pl,chisq_res2,
                   by.x="Metal",
                   by.y="Metal",all.x=T) -> Fig.metal.growth.pl  

#13.4 Metal for psychrophilic####
  host_pc_metaBacMet$Metal <-
    ifelse(is.na(host_pc_metaBacMet$Metal),"Unknown", 
           host_pc_metaBacMet$Metal)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Metal,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_csp_metal
  
  host_pc_csp_metal %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> metal_sum
  metal <- unique(host_pc_csp_metal$Metal)
  host_pc_csp_metal %>% 
    filter(trophic2 == "Csp" & 
             Metal %in% metal) -> a
  
  host_pc_csp_metal %>% 
    filter(trophic2 != "Csp" & 
             Metal %in% metal) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Metal %in%
                                        b$Metal],
                      as.numeric(b$n)))
  colnames(M) <- b$Metal[b$Metal %in% a$Metal]
  
  metal2 <- metal[a$Metal %in% b$Metal]
  chisq_res <- list()
  for(i in 1:length(metal2))
  {
    M2 <- as.table(rbind(M[,metal2[i]],
                         metal_sum$S-M[,metal2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- metal2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Metal","p.value")
  
  host_pc_csp_metal %>% 
    #filter(Metal != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  Metal=.$Metal,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_csp_metal2
  host_pc_csp_metal2$freq <-
    host_pc_csp_metal2$n/host_pc_csp_metal2$S
  
  host_pc_csp_metal2 %>%
    group_by(Metal) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_csp_metal2.maxfreq
  host_pc_csp_metal2.maxfreq <-
    as.data.frame(host_pc_csp_metal2.maxfreq)
  rownames(host_pc_csp_metal2.maxfreq) <- 
    host_pc_csp_metal2.maxfreq$Metal
  chisq_res$y <- host_pc_csp_metal2.maxfreq[
    chisq_res$Metal,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Metal != "Unknown")
  host_pc_csp_metal2 %>%
    filter( Metal != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_csp_metal2.pl
  
  chisq_res$Metal <- factor(chisq_res$Metal,
          levels = levels(host_pc_csp_metal2.pl$Metal))
  
  host_pc_csp_metal2.pl$trophic2 <- 
    factor(host_pc_csp_metal2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_csp_metal2.pl,
         aes(x=as.character.factor(Metal),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",
             colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Metal,
                                  y=(y+0.0005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("green","black"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.metal.trophic2
  
  merge.data.frame(host_pc_csp_metal2.pl,chisq_res2,
                   by.x="Metal",
                   by.y="Metal",all.x=T) ->
    Fig.metal.csp.pl  

#13.5 Metal all####
  Fig.metal.ogt.pl$type="OGT"
  Fig.metal.growth.pl$type="Growth rate"
  Fig.metal.csp.pl$type="CSP"
  colnames(Fig.metal.ogt.pl)[3] <- "Var"
  colnames(Fig.metal.growth.pl)[3] <- "Var"
  colnames(Fig.metal.csp.pl)[3] <- "Var"
  Fig.metal.pl <- rbind.data.frame(Fig.metal.ogt.pl,
                                  Fig.metal.growth.pl)
  Fig.metal.pl <- rbind.data.frame(Fig.metal.pl,
                                  Fig.metal.csp.pl)
  Fig.metal.pl$type <- 
    factor(Fig.metal.pl$type,
           levels = c("OGT","Growth rate","CSP"))
  Fig.metal.pl %>% ggplot(aes(x=Metal,y=freq*100))+
    geom_col(aes(fill=Var),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(aes(Metal,y=(y+0.0005)*100,label=sig),angle=90)+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1),breaks = seq(0,1,0.2)) + 
    coord_flip()+
    facet_wrap(~type,nrow=1)+
    labs(x="Genes of heavy metal resistances",
         y="Proportion %")+
    scale_fill_manual(values = c("red","grey","yellow",
                                 "white","green","black"))+
    guides(fill=F) +
    theme(axis.text.x = element_text(size=10#,angle=45,
                                     #hjust=1,vjust=1
                                     )
    )-> Fig.Metal.all
  
  ggsave("../../04.Figure/01.Raw/Fig.Metal.all.tiff",
         plot= Fig.Metal.all,device = "tiff",compression="lzw",
         width = 6.2,height = 3.5*2,units = "in")

#14.1 Antibio for OGT####
  host_pc_metaBacMet$antibio <- 
    as.character.factor(host_pc_metaBacMet$antibio)
  host_pc_metaBacMet$antibio <-
    ifelse(is.na(host_pc_metaBacMet$antibio),"Unknown", 
           host_pc_metaBacMet$antibio)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(antibio,trophic) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic) -> host_pc_ogt_antibio
  
  host_pc_ogt_antibio %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n))) -> antibio_sum
  antibio <- unique(host_pc_ogt_antibio$antibio)
  host_pc_ogt_antibio %>% 
    filter(trophic == "Others" & 
             antibio %in% antibio) -> a
  
  host_pc_ogt_antibio %>% 
    filter(trophic != "Others" & 
             antibio %in% antibio) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$antibio %in%
                                        b$antibio],
                      as.numeric(b$n)[b$antibio %in%
                                        a$antibio] ))
  colnames(M) <- b$antibio[b$antibio %in% a$antibio]
  
  antibio2 <- antibio[a$antibio %in% b$antibio]
  chisq_res <- list()
  for(i in 1:length(antibio2))
  {
    M2 <- as.table(rbind(M[,antibio2[i]],
                         antibio_sum$S-M[,antibio2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- antibio2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("antibio","p.value")
  
  host_pc_ogt_antibio %>% 
    # filter(Metal != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n),
                  antibio=.$antibio,
                  trophic=.$trophic,
                  n=.$n
    )) -> host_pc_ogt_antibio2
  host_pc_ogt_antibio2$freq <-
    host_pc_ogt_antibio2$n/host_pc_ogt_antibio2$S
  
  host_pc_ogt_antibio2 %>%
    group_by(antibio) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_ogt_antibio2.maxfreq
  host_pc_ogt_antibio2.maxfreq <-
    as.data.frame(host_pc_ogt_antibio2.maxfreq)
  rownames(host_pc_ogt_antibio2.maxfreq) <- 
    host_pc_ogt_antibio2.maxfreq$antibio
  chisq_res$y <- 
    host_pc_ogt_antibio2.maxfreq[chisq_res$antibio,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(antibio != "Unknown")
  host_pc_ogt_antibio2 %>%
    filter( antibio != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_antibio2.pl
  
  chisq_res$antibio <- factor(chisq_res$antibio,
                            levels = levels(host_pc_ogt_antibio2.pl$antibio))
  
  host_pc_ogt_antibio2.pl$trophic <- 
    factor(host_pc_ogt_antibio2.pl$trophic,
           levels = c("Thermophilic","Others"))
  
  ggplot(host_pc_ogt_antibio2.pl,
         aes(x=as.character.factor(antibio),y=freq*100))+
    geom_col(aes(fill=trophic),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=antibio,
                                  y=(y+0.001)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.5)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.antibio.ogt
  
  merge.data.frame(host_pc_ogt_antibio2.pl,chisq_res2,
                   by.x="antibio",
                   by.y="antibio",all.x=T) ->
    Fig.antibio.ogt.pl  
  
#14.2  Antibio for growth####
  host_pc_metaBacMet$antibio <-
    ifelse(is.na(host_pc_metaBacMet$antibio),"Unknown", 
           host_pc_metaBacMet$antibio)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(antibio,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(antibio) -> host_pc_growth_antibio
  
  host_pc_growth_antibio %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> antibio_sum
  antibio <- unique(host_pc_growth_antibio$antibio)
  host_pc_growth_antibio %>% 
    filter(growth == "Fast" & 
             antibio %in% antibio) -> a
  
  host_pc_growth_antibio %>% 
    filter(growth != "Fast" & 
             antibio %in% antibio) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$antibio %in%
                                        b$antibio],
                      as.numeric(b$n)))
  colnames(M) <- b$antibio[b$antibio %in% a$antibio]
  
  antibio2 <- antibio[a$antibio %in% b$antibio]
  chisq_res <- list()
  for(i in 1:length(antibio2))
  {
    M2 <- as.table(rbind(M[,antibio2[i]],
                         antibio_sum$S-M[,antibio2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- antibio2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("antibio","p.value")
  
  host_pc_growth_antibio %>% 
    #filter(Metal != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  antibio=.$antibio,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_growth_antibio2
  host_pc_growth_antibio2$freq <-
    host_pc_growth_antibio2$n/host_pc_growth_antibio2$S
  
  host_pc_growth_antibio2 %>%
    group_by(antibio) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_growth_antibio2.maxfreq
  host_pc_growth_antibio2.maxfreq <-
    as.data.frame(host_pc_growth_antibio2.maxfreq)
  rownames(host_pc_growth_antibio2.maxfreq) <- 
    host_pc_growth_antibio2.maxfreq$antibio
  chisq_res$y <- host_pc_growth_antibio2.maxfreq[
    chisq_res$antibio,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(antibio != "Unknown")
  host_pc_growth_antibio2 %>%
    filter( antibio != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_growth_antibio2.pl
  
  chisq_res$antibio <- factor(chisq_res$antibio,
          levels = levels(host_pc_growth_antibio2.pl$antibio))
  
  host_pc_growth_antibio2.pl$growth <- 
    factor(host_pc_growth_antibio2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_growth_antibio2.pl,
         aes(x=as.character.factor(antibio),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",
             colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=antibio,
                                  y=(y+0.0005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.antibio.growth
  
  merge.data.frame(host_pc_growth_antibio2.pl,chisq_res2,
                   by.x="antibio",
                   by.y="antibio",all.x=T) -> 
    Fig.antibio.growth.pl  
  
#14.3  Antibio for psychrophilic####
  host_pc_metaBacMet$antibio <-
    ifelse(is.na(host_pc_metaBacMet$antibio),"Unknown", 
           host_pc_metaBacMet$antibio)
  
  host_pc_metaBacMet %>% filter(!is.na(host_vc_R2)) %>%
    group_by(antibio,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_csp_antibio
  
  host_pc_csp_antibio %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> antibio_sum
  antibio <- unique(host_pc_csp_antibio$antibio)
  host_pc_csp_antibio %>% 
    filter(trophic2 == "Csp" & 
             antibio %in% antibio) -> a
  
  host_pc_csp_antibio %>% 
    filter(trophic2 != "Csp" & 
             antibio %in% antibio) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$antibio %in%
                                        b$antibio],
                      as.numeric(b$n)))
  colnames(M) <- b$antibio[b$antibio %in% a$antibio]
  
  antibio2 <- antibio[a$antibio %in% b$antibio]
  chisq_res <- list()
  for(i in 1:length(antibio2))
  {
    M2 <- as.table(rbind(M[,antibio2[i]],
                         antibio_sum$S-M[,antibio2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- antibio2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("antibio","p.value")
  
  host_pc_csp_antibio %>% 
    #filter(Metal != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  antibio=.$antibio,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_csp_antibio2
  host_pc_csp_antibio2$freq <-
    host_pc_csp_antibio2$n/host_pc_csp_antibio2$S
  
  host_pc_csp_antibio2 %>%
    group_by(antibio) %>%
    do(data.frame(Max=max(.$freq))) ->
    host_pc_csp_antibio2.maxfreq
  host_pc_csp_antibio2.maxfreq <-
    as.data.frame(host_pc_csp_antibio2.maxfreq)
  rownames(host_pc_csp_antibio2.maxfreq) <- 
    host_pc_csp_antibio2.maxfreq$antibio
  chisq_res$y <- host_pc_csp_antibio2.maxfreq[
    chisq_res$antibio,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(antibio != "Unknown")
  host_pc_csp_antibio2 %>%
    filter( antibio != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_csp_antibio2.pl
  
  chisq_res$antibio <- factor(chisq_res$antibio,
       levels = levels(host_pc_csp_antibio2.pl$antibio))
  
  host_pc_csp_antibio2.pl$trophic2 <- 
    factor(host_pc_csp_antibio2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_csp_antibio2.pl,
         aes(x=as.character.factor(antibio),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",
             colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=antibio,
                                  y=(y+0.0005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Presence","Absence"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.antibio.trophic2
  
  merge.data.frame(host_pc_csp_antibio2.pl,chisq_res2,
                   by.x="antibio",
                   by.y="antibio",all.x=T) ->
    Fig.antibio.csp.pl  
  
#14.4  Antibio all####
  Fig.antibio.ogt.pl$type="OGT"
  Fig.antibio.growth.pl$type="Growth rate"
  Fig.antibio.csp.pl$type="CSP"
  colnames(Fig.antibio.ogt.pl)[3] <- "Var"
  colnames(Fig.antibio.growth.pl)[3] <- "Var"
  colnames(Fig.antibio.csp.pl)[3] <- "Var"
  Fig.antibio.pl <- rbind.data.frame(Fig.antibio.ogt.pl,
                                   Fig.antibio.growth.pl)
  Fig.antibio.pl <- rbind.data.frame(Fig.antibio.pl,
                                   Fig.antibio.csp.pl)
  Fig.antibio.pl$type <- 
    factor(Fig.antibio.pl$type,
           levels = c("OGT","Growth rate","CSP"))
  Fig.antibio.pl %>% 
    ggplot(aes(x=antibio,y=freq*100))+
    geom_col(aes(fill=Var),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(aes(antibio,y=(y+0.002)*100,label=sig),angle=90)+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.1),
                       breaks = seq(0,1,0.2)) + 
    coord_flip()+
    facet_wrap(~type,nrow=1)+
    labs(x="Genes of antibiotic resistances",y="Proportion %")+
    scale_fill_manual(values = c("red","grey","yellow",
                                 "white","green","black"))+
    guides(fill=F) +
    theme(axis.text = element_text(size=10#,angle=45,
                                     #hjust=1,vjust=1
                                     )
    )-> Fig.antibio.all
  
  ggsave("../../04.Figure/01.Raw/Fig.antibio.all.tiff",
         plot= Fig.antibio.all,device = "tiff",compression="lzw",
         width = 6.6*1.3,height = 3.5*2.3,units = "in")
  
#15.1 SCycDB####
  iso_pc_funcs_SCyc <- 
    iso_pc_funcs[,c("pc","Scyc_ID","Scyc_gene")]
  host_pc_metaS <- merge.data.frame(host_pc,pc_host_metadata,
                                    by.x="host",by.y = "Assembly.Accession",all.x=T)
  host_pc_metaS <- merge.data.frame(host_pc_metaS,pc_host_csp,
                                    by.x="host",by.y="Assembly.Accession",all.x=T)
  host_pc_metaS <- merge.data.frame(host_pc_metaS,iso_pc_funcs_SCyc,
                                    by="pc",all.x=T)
  host_pc_metaS$num_csp <- ifelse(is.na(host_pc_metaS$num_csp),
                                  0,host_pc_metaS$num_csp)
  host_pc_metaS$trophic <-
    ifelse(host_pc_metaS$OGT >= 40, "Thermophilic",
           "Others")
  host_pc_metaS$growth <- 
    ifelse(host_pc_metaS$DoublingTime>=2.5,"Slow","Fast")
  host_pc_metaS$num_csp <- ifelse(is.na(host_pc_metaS$num_csp),
                                  0,host_pc_metaS$num_csp)
  host_pc_metaS$trophic2 <- 
    ifelse(host_pc_metaS$num_csp!=0,"Csp","No csp")
#15.2 SCyc for OGT ####
  host_pc_metaS %>% filter(!is.na(Scyc_ID)) %>% dim()
  
  host_pc_metaS$Scyc_gene <- 
    as.character.factor(host_pc_metaS$Scyc_gene)
  host_pc_metaS$Scyc_gene <-
    ifelse(is.na(host_pc_metaS$Scyc_gene),"Unknown", 
           host_pc_metaS$Scyc_gene)
  
  host_pc_metaS %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Scyc_gene,trophic) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic) -> host_pc_ogt_scyc
  
  host_pc_ogt_scyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n))) -> scyc_sum
  scyc <- unique(host_pc_ogt_scyc$Scyc_gene)
  host_pc_ogt_scyc %>% 
    filter(trophic == "Others" & 
             Scyc_gene %in% scyc) -> a
  
  host_pc_ogt_scyc %>% 
    filter(trophic != "Others" & 
             Scyc_gene %in% scyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Scyc_gene %in% b$Scyc_gene],
                      as.numeric(b$n)))
  colnames(M) <- b$Scyc_gene
  
  scyc2 <- scyc[a$Scyc_gene %in% b$Scyc_gene]
  chisq_res <- list()
  for(i in 1:length(scyc2))
  {
    M2 <- as.table(rbind(M[,scyc2[i]],scyc_sum$S-M[,scyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- scyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Scyc_gene","p.value")
  
  host_pc_ogt_scyc %>% 
    #     filter(Ncyc_gene != "Unknown") %>%
    group_by(trophic) %>% 
    do(data.frame(S=sum(.$n),
                  Scyc_gene=.$Scyc_gene,
                  trophic=.$trophic,
                  n=.$n
    )) -> host_pc_ogt_scyc2
  host_pc_ogt_scyc2$freq <- host_pc_ogt_scyc2$n/host_pc_ogt_scyc2$S
  
  host_pc_ogt_scyc2 %>%
    group_by(Scyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_scyc2.maxfreq
  host_pc_ogt_scyc2.maxfreq <-
    as.data.frame(host_pc_ogt_scyc2.maxfreq)
  rownames(host_pc_ogt_scyc2.maxfreq) <- 
    host_pc_ogt_scyc2.maxfreq$Scyc_gene
  chisq_res$y <- host_pc_ogt_scyc2.maxfreq[chisq_res$Scyc_gene,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Scyc_gene != "Unknown")
  host_pc_ogt_scyc2 %>%
    filter( Scyc_gene != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_scyc2.pl
  
  chisq_res$Scyc_gene <- factor(chisq_res$Scyc_gene,
                                levels = levels(host_pc_ogt_scyc2.pl$Scyc_gene))
  
  host_pc_ogt_scyc2.pl$trophic <- 
    factor(host_pc_ogt_scyc2.pl$trophic,
           levels = c("Thermophilic","Others"))
  
  ggplot(host_pc_ogt_scyc2.pl,
         aes(x=as.character.factor(Scyc_gene),y=freq*100))+
    geom_col(aes(fill=trophic),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Scyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="OGT"))+
    scale_fill_manual(values = c("red","grey"),
                      labels=c("≥40","Others"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.scyc.ogt
  
  merge.data.frame(host_pc_ogt_scyc2.pl,chisq_res2,
                   by.x="Scyc_gene",
                   by.y="Scyc_gene",all.x=T) -> Fig.scyc.ogt.pl  
  
  #15.3 SCyc for growth ####
  #host_pc_metaN$Ncyc_gene <- 
  #  as.character.factor(host_pc_metaN$Ncyc_gene)
  host_pc_metaS$Scyc_gene <-
    ifelse(is.na(host_pc_metaS$Scyc_gene),"Unknown", 
           host_pc_metaS$Scyc_gene)
  
  host_pc_metaS %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Scyc_gene,growth) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(growth) -> host_pc_ogt_scyc
  
  host_pc_ogt_scyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n))) -> scyc_sum
  scyc <- unique(host_pc_ogt_scyc$Scyc_gene)
  host_pc_ogt_scyc %>% 
    filter(growth == "Fast" & 
             Scyc_gene %in% scyc) -> a
  
  host_pc_ogt_scyc %>% 
    filter(growth != "Fast" & 
             Scyc_gene %in% scyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Scyc_gene %in% b$Scyc_gene],
                      as.numeric(b$n)[b$Scyc_gene %in% a$Scyc_gene] ))
  colnames(M) <- b$Scyc_gene[b$Scyc_gene %in% a$Scyc_gene]
  
  scyc2 <- b$Scyc_gene[b$Scyc_gene %in% a$Scyc_gene]
  chisq_res <- list()
  for(i in 1:length(scyc2))
  {
    M2 <- as.table(rbind(M[,scyc2[i]],scyc_sum$S-M[,scyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- scyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Scyc_gene","p.value")
  
  host_pc_ogt_scyc %>% 
    #  filter(Ncyc_gene != "Unknown") %>%
    group_by(growth) %>% 
    do(data.frame(S=sum(.$n),
                  Scyc_gene=.$Scyc_gene,
                  growth=.$growth,
                  n=.$n
    )) -> host_pc_ogt_scyc2
  host_pc_ogt_scyc2$freq <- host_pc_ogt_scyc2$n/host_pc_ogt_scyc2$S
  
  host_pc_ogt_scyc2 %>%
    group_by(Scyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_scyc2.maxfreq
  host_pc_ogt_scyc2.maxfreq <-
    as.data.frame(host_pc_ogt_scyc2.maxfreq)
  rownames(host_pc_ogt_scyc2.maxfreq) <- 
    host_pc_ogt_scyc2.maxfreq$Scyc_gene
  chisq_res$y <- host_pc_ogt_scyc2.maxfreq[chisq_res$Scyc_gene,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Scyc_gene != "Unknown")
  host_pc_ogt_scyc2 %>%
    filter( Scyc_gene != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_scyc2.pl
  
  chisq_res$Scyc_gene <- factor(chisq_res$Scyc_gene,
                                levels = levels(host_pc_ogt_scyc2.pl$Scyc_gene))
  
  host_pc_ogt_scyc2.pl$gr <- 
    factor(host_pc_ogt_scyc2.pl$growth,
           levels = c("Fast","Slow"))
  
  ggplot(host_pc_ogt_scyc2.pl,
         aes(x=as.character.factor(Scyc_gene),y=freq*100))+
    geom_col(aes(fill=growth),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Scyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="Growth rate"))+
    scale_fill_manual(values = c("yellow","white"),
                      labels=c("Fast","Slow"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.scyc.growth
  
  merge.data.frame(host_pc_ogt_scyc2.pl,chisq_res2,
                   by.x="Scyc_gene",
                   by.y="Scyc_gene",all.x=T) -> Fig.scyc.growth.pl  
  
  #15.4 SCyc for psychrophilic ####
  # host_pc_metaN$Ncyc_gene <- 
  #   as.character.factor(host_pc_metaN$Ncyc_gene)
  host_pc_metaS$Scyc_gene <-
    ifelse(is.na(host_pc_metaS$Scyc_gene),"Unknown", 
           host_pc_metaS$Scyc_gene)
  
  host_pc_metaS %>% filter(!is.na(host_vc_R2)) %>%
    group_by(Scyc_gene,trophic2) %>%
    do(data.frame(n=sum(.$n))) %>%
    group_by(trophic2) -> host_pc_ogt_scyc
  
  host_pc_ogt_scyc %>% 
    #  filter(COG.Functional.cat. != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n))) -> scyc_sum
  scyc <- unique(host_pc_ogt_scyc$Scyc_gene)
  host_pc_ogt_scyc %>% 
    filter(trophic2 == "Csp" & 
             Scyc_gene %in% scyc) -> a
  
  host_pc_ogt_scyc %>% 
    filter(trophic2 != "Csp" & 
             Scyc_gene %in% scyc) -> b
  
  M <- as.table(rbind(as.numeric(a$n)[a$Scyc_gene %in% b$Scyc_gene],
                      as.numeric(b$n)[b$Scyc_gene %in% a$Scyc_gene] ))
  colnames(M) <- b$Scyc_gene[b$Scyc_gene %in% a$Scyc_gene]
  
  scyc2 <- b$Scyc_gene[b$Scyc_gene %in% a$Scyc_gene]
  chisq_res <- list()
  for(i in 1:length(scyc2))
  {
    M2 <- as.table(rbind(M[,scyc2[i]],scyc_sum$S-M[,scyc2[i]]))
    chisq_res[[i]] <- chisq.test(M2)$p.value
  }
  names(chisq_res) <- scyc2
  chisq_res <- ldply(chisq_res)
  colnames(chisq_res) <- c("Scyc_gene","p.value")
  
  host_pc_ogt_scyc %>% 
    #  filter(Ncyc_gene != "Unknown") %>%
    group_by(trophic2) %>% 
    do(data.frame(S=sum(.$n),
                  Scyc_gene=.$Scyc_gene,
                  trophic2=.$trophic2,
                  n=.$n
    )) -> host_pc_ogt_scyc2
  host_pc_ogt_scyc2$freq <- host_pc_ogt_scyc2$n/host_pc_ogt_scyc2$S
  
  host_pc_ogt_scyc2 %>%
    group_by(Scyc_gene) %>%
    do(data.frame(Max=max(.$freq))) -> host_pc_ogt_scyc2.maxfreq
  host_pc_ogt_scyc2.maxfreq <-
    as.data.frame(host_pc_ogt_scyc2.maxfreq)
  rownames(host_pc_ogt_scyc2.maxfreq) <- 
    host_pc_ogt_scyc2.maxfreq$Scyc_gene
  chisq_res$y <- host_pc_ogt_scyc2.maxfreq[chisq_res$Scyc_gene,"Max"]
  chisq_res$sig <- ""
  #  chisq_res$sig <- ifelse(chisq_res$p.value < 0.1,"·",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.05,"*",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.01,"**",chisq_res$sig)
  chisq_res$sig <- ifelse(chisq_res$p.value < 0.001,"***",chisq_res$sig)
  
  #chisq_res %>% filter(!M2 %in% c("Overview","Unknown"))  -> chisq_res2
  chisq_res2 <- chisq_res %>% filter(Scyc_gene != "Unknown")
  host_pc_ogt_scyc2 %>%
    filter( Scyc_gene != "Unknown"
            # & Ncyc_gene %in% b$Ncyc_gene
    ) -> host_pc_ogt_scyc2.pl
  
  chisq_res$Scyc_gene <- factor(chisq_res$Scyc_gene,
                                levels = levels(host_pc_ogt_scyc2.pl$Scyc_gene))
  
  host_pc_ogt_scyc2.pl$gr <- 
    factor(host_pc_ogt_scyc2.pl$trophic2,
           levels = c("Csp","No csp"))
  
  ggplot(host_pc_ogt_scyc2.pl,
         aes(x=as.character.factor(Scyc_gene),y=freq*100))+
    geom_col(aes(fill=trophic2),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    coord_flip()+
    geom_text(data=chisq_res2,aes(x=Scyc_gene,
                                  y=(y+0.00005)*100,label=sig))+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    labs(x="KEGG categories",y="Porprotion %")+
    guides(fill=guide_legend(title="CSP"))+
    scale_fill_manual(values = c("green","black"),
                      labels=c("Presence","Absence"))+
    theme(legend.position = c(0.85,0.8)#,
          #          axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    ) -> Fig.scyc.csp
  
  merge.data.frame(host_pc_ogt_scyc2.pl,chisq_res2,
                   by.x="Scyc_gene",
                   by.y="Scyc_gene",all.x=T) -> Fig.scyc.csp.pl 
  
  #15.5 SCyc all####
  Fig.scyc.ogt.pl$type="OGT"
  Fig.scyc.growth.pl$type="Growth rate"
  Fig.scyc.csp.pl$type="CSP"
  colnames(Fig.scyc.ogt.pl)[3] <- "Var"
  colnames(Fig.scyc.growth.pl)[3] <- "Var"
  colnames(Fig.scyc.csp.pl)[3] <- "Var"
  Fig.scyc.pl <- rbind.data.frame(Fig.scyc.ogt.pl,
                                  Fig.scyc.growth.pl[,-6])
  Fig.scyc.pl <- rbind.data.frame(Fig.scyc.pl,
                                  Fig.scyc.csp.pl[,-6])
  Fig.scyc.pl$type <- 
    factor(Fig.scyc.pl$type,levels = c("OGT","Growth rate","CSP"))
  Fig.scyc.pl %>% 
    filter(type=="OGT") %>%
    group_by(Scyc_gene) %>% 
    dplyr::summarise(gs=mean(freq) ) -> sa
  sa %>% filter(gs >=0.0001) -> sb

  sg1 <- str_extract(levels(Fig.scyc.pl$Scyc_gene),"[a-z]+")
  sg2 <- str_remove(levels(Fig.scyc.pl$Scyc_gene),"[a-z]+")
  sg3 <- NULL
  sg3.2 <- NULL
  for(i in 1:length(sg1))
  {
    
    sg3 <- c(sg3,
             substitute(expression(italic(g1)*g2),
                        list(g1=sg1[i],g2=sg2[i]))
             )
    sg3.2 <- c(sg3.2,
               paste("expression(italic(",Fig.scyc.pl$g1[i],")*",
                       Fig.scyc.pl$g2[i],")",sep=""
               ))
  }
  Fig.scyc.pl %>% 
    #filter(p.value <0.05) %>%
    filter(Scyc_gene %in% sb$Scyc_gene) %>%
    ggplot(aes(x=Scyc_gene,y=freq*100))+
    geom_col(aes(fill=Var),position = "dodge",colour="black")+
    #geom_text(x="A",y=0.02,label="11")
    geom_text(aes(Scyc_gene,y=(y+0.001)*100,label=sig),
              angle=90,size=4)+
    mytheme +
    scale_y_continuous(expand = c(0,0,0,0.05)) + 
    coord_flip()+
    facet_wrap(~type,nrow=1)+
    labs(x="Genes of sulfur cycle",y="Proportion %")+
    scale_fill_manual(values = c("red","grey","yellow",
                                 "white","green","black"))+
    guides(fill=F) +
    #scale_x_discrete(labels=sg3)+
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(face="italic")
    )-> Fig.Scyc.all
  
  ggsave("../../04.Figure/01.Raw/Fig.Scyc.all.tiff",
         plot= Fig.Scyc.all,device = "tiff",compression="lzw",
         width = 6.6,height = 3.5*2*2.5,units = "in")
   ggsave("../../04.Figure/01.Raw/Fig.Scyc.all.pdf",
         plot= Fig.Scyc.all,device = "pdf",
         width = 6.6,height = 3.5*2*2.5)

#Multiplicity of infection vs.  d'####
  ggplot(iso_host_metadata_csp,
         aes(x=`Num Temperate`+`Num Virulent`,y=host_vc_R2))+
    geom_point(alpha=0.1)+
    geom_smooth(method="lm",se=F,linetype="dashed")+mytheme+
    labs(x="Multiplicity of infections",
         y=expression(Specialisation~italic("d'")))+
    scale_x_continuous(breaks = seq(1,15,3)) -> Fig.MOI_d
  ggsave("../../04.Figure/01.Raw/Fig.MOI_d.tiff",
         plot= Fig.MOI_d,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  ggplot(iso_host_metadata_csp,
         aes(x=num_VC,y=host_vc_R2))+
    geom_point(alpha=0.1)+
    geom_smooth(method="lm",se=F,linetype="dashed")+mytheme+
    labs(x="Number of viral clusters",
         y=expression(Specialisation~italic("d'")))+
    scale_x_continuous(breaks = seq(1,15,3)) -> Fig.VC_d
  ggsave("../../04.Figure/01.Raw/Fig.VC_d.tiff",
         plot= Fig.VC_d,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  ggplot(iso_host_metadata_csp,
         aes(x=`Num Temperate`+`Num Virulent`,y=num_VC))+
    geom_point(alpha=0.1)+
    geom_smooth(method="lm",se=F)+mytheme+
    labs(x="Multiplicity of infections",
         y="Number of viral clusters")+
    scale_x_continuous(breaks = seq(1,15,3)) +
    scale_y_continuous(breaks = seq(1,15,3)) + 
    geom_text(aes(x=5,y=11),
              label=expression("y = 1.02x - 0.02,"~R^2~"= 0.93, p < 0.001"),
              size=5)-> Fig.MOI_VC
  ggsave("../../04.Figure/01.Raw/Fig.MOI_VC.tiff",
         plot= Fig.MOI_VC,device = "tiff",compression="lzw",
         width = 6.6,height = 4.8,units = "in")
  
#Net dis vs. Trait dis ####
load("../../02.Data/02.Rdata/1.host_vc_M_muti.Rdata")
save(host_vc_M,file="host_vc_M.Rdata",version = 2)
host_vc_dis <- vegdist(host_vc_M$iso_host_vc_R2,method="gower",binary = T)
host_vc_dis_m <- as.matrix(host_vc_dis)

#iso_host_metadata2 <- iso_host_metadata[iso_host_metadata$DoublingTime < 320,]
iso_host_metadata2 <- iso_host_metadata
rownames(iso_host_metadata2) <- iso_host_metadata2$Assembly.Accession
iso_host_metadata2 <- 
  iso_host_metadata2[rownames(host_vc_M$iso_host_vc_R2),]
iso_host_metadata2 <- 
  iso_host_metadata2[-which(is.na(iso_host_metadata2$Size_Mb)),]
iso_host_vc_R2.2 <-
  host_vc_M$iso_host_vc_R2[rownames(iso_host_metadata2),]

genome_traits <- c(
  "Size_Mb","GC","num_5S","num_16S","num_23S",
  "num_sigma_factor","OGT","DoublingTime","num_prophage")
Genetic_csp <- colnames(iso_host_metadata2)[
  grep("Cold_shock_protein",
                    colnames(iso_host_metadata2))]
Genetic_hsp <- colnames(iso_host_metadata2)[
  grep("Heat shock protein",
                    colnames(iso_host_metadata2))]
Genetic_receptor <- colnames(iso_host_metadata2)[
  grep("Prot_Recep",
                    colnames(iso_host_metadata2))]
Genetic_Cas <- colnames(iso_host_metadata2)[
  grep("CAS-",
              colnames(iso_host_metadata2))[c(1:15)]]
Genetic_RM <- colnames(iso_host_metadata2)[
  grep("RM",
            colnames(iso_host_metadata2))]

Genetic_switches <- c("CI","Cro","fpsA","fpsR")

Genetic_antiCRISPR <- colnames(iso_host_metadata2)[
  grep("AntiCRISPR",
                   colnames(iso_host_metadata2))]

Genetic_traits <- list(
  Genome_traits=genome_traits,
  Genetic_csp=Genetic_csp,
  Genetic_hsp=Genetic_hsp,
  Genetic_receptor=Genetic_receptor,
  Genetic_immunity=c(Genetic_Cas,Genetic_RM),
  Genetic_antiCRISPR=Genetic_antiCRISPR,
  Genetic_switches=Genetic_switches
)
ZScores <- c(T,rep(F,6))

mantel_res <- list()
for(i in 1:length(Genetic_traits)){
  mantel_res[[i]] <-
    mantel_genomic_trait_net(Metadata=iso_host_metadata2,
                             NetMatrix=iso_host_vc_R2.2,
                             used.traits=Genetic_traits[[i]],
                             group.trait=NULL,
                             method.1="gower",
                             method.2="euclidean",
                             Zscores=ZScores[i])
}
names(mantel_res) <- names(Genetic_traits)
save(mantel_res,file="../../02.Data/02.Rdata/3.mantel_res.Rdata")
ldply(mantel_res)


group.traits <- c("phylum","gram_stain","cell_shape",
                  "motility","sporulation","salinity","oxygen_requirement",
                  "empo_2")

group_mantel_res1 <- list()
k <- 1
for(i in 1:length(group.traits)){
  group_mantel_res2 <- list()
  for(j in 1:length(Genetic_traits))
  {
    group_mantel_res2[[j]] <- 
      mantel_genomic_trait_net(Metadata=iso_host_metadata2,
                               NetMatrix=iso_host_vc_R2.2,
                               used.traits=Genetic_traits[[j]],
                               group.trait=group.traits[i],
                               method.1="gower",
                               method.2="euclidean",
                               Zscores=ZScores[j])
    
  }
  names(group_mantel_res2) <- names(Genetic_traits)
  group_mantel_res1[[i]] <- ldply(group_mantel_res2,.id="Genetic_traits")
}
names(group_mantel_res1) <- group.traits
group_mantel_res1.df <- ldply(group_mantel_res1,.id="Groups")
group_mantel_res1.df %>% filter(!is.na(p)) -> group_mantel_res1.df2

group_mantel_res1.df2$sig <- 
  ifelse(group_mantel_res1.df2$p > 0.05,"","*")
group_mantel_res1.df2$sig <- 
  ifelse(group_mantel_res1.df2$p > 0.01,
         group_mantel_res1.df2$sig,"**")
group_mantel_res1.df2$sig <- 
  ifelse(group_mantel_res1.df2$p > 0.001,
         group_mantel_res1.df2$sig,"***")

group_mantel_res1.df2 %>%
  ggplot(aes(y=Groups,x=Genetic_traits))+
  geom_tile(aes(fill=r))+
  geom_text(aes(label=sig))+
  mytheme+
  labs(x="")+
  scale_fill_gradient2(low="#99CCFF",high="#FF6666")+
  theme(axis.text.x = element_text(angle=60,hjust=1))


#dist of d' & genomic traits
rownames(iso_host_stat_all) <- iso_host_stat_all$Assembly.Accession
d_R2 <- iso_host_stat_all[rownames(iso_host_metadata2),
                          "host_vc_R2"]
d_R2 <- data.frame(d=d_R2)
rownames(d_R2) <- rownames(iso_host_metadata2)

#scatter plot


group_dis1 <- list()
k <- 1
for(i in 1:length(group.traits)){
  group_dis2 <- list()
  for(j in 1:length(Genetic_traits))
  {
    group_dis2[[j]] <- 
      trait_net_dis(Metadata=iso_host_metadata2,
                              # NetMatrix=iso_host_vc_R2.2,
                               NetMatrix=d_R2,
                               used.traits=Genetic_traits[[j]],
                               group.trait=group.traits[i],
                               #method.1="altGower",
                               method.1="euclidean",
                               method.2="euclidean",
                               Zscores=ZScores[j],
                    m1.binary=F)
    
  }
  names(group_dis2) <- names(Genetic_traits)
  group_dis1[[i]] <- ldply(group_dis2,.id="Genetic_traits")
}
names(group_dis1) <- group.traits
group_dis1.df <- ldply(group_dis1,.id="Groups0")

group_dis1.df %>% 
  group_by(Genetic_traits,Groups) %>%
  do(data.frame(
    p=cal_lm_p(.$net.dis,.$trait.dis)[2]
    )) -> dis.lm.p
dis.lm.p <- as.data.frame(dis.lm.p)
rownames(dis.lm.p) <- paste(dis.lm.p$Genetic_traits,
                            dis.lm.p$Groups,sep="_")
dis.lm.p$sig <- ifelse(dis.lm.p$p > 0.05,"ns","sig")

group_dis1.df$ID <- paste(group_dis1.df$Genetic_traits,
                          group_dis1.df$Groups,sep="_")
group_dis1.df$sig <- dis.lm.p[group_dis1.df$ID,"sig"]
group_dis1.df$sig <- 
  factor(group_dis1.df$sig,levels = c("sig","ns"))

group_dis1.df %>% 
  ggplot(aes(x=log(trait.dis),y=net.dis))+
  geom_point(size=0.05,alpha=0.3)+
  geom_smooth(aes(linetype=sig),method="lm")+
  mytheme+
  facet_grid(Groups~Genetic_traits,
             scales = "free_x") -> p.trait_net_dis

ggsave("trait_net_dis.tif2",
     plot=p.trait_net_dis,device = "tiff",
     compression="lzw",
     width = 6.6*1.5,height = 4.8*8,units = "in")



#Negative DT-d'####
rownames(iso_host_metadata_csp) <- 
  iso_host_metadata_csp$Assembly.Accession
csp_m <- iso_host_metadata_csp[,csp_id[1:16]]

csp_m[is.na(csp_m)] <- 0
rownames(csp_m) <- iso_host_metadata_csp$Assembly.Accession
csp_m <- csp_m[-which(rowSums(csp_m)==0),]
#csp_pca <- rda(csp_m)
csp_pca <- decorana(csp_m)
#csp_pca.plot <- as.data.frame(csp_pca$CA$u)
csp_pca.plot <- as.data.frame(csp_pca$rproj)

csp_pca.plot$order <- iso_host_metadata_csp[rownames(csp_m),"order"]
csp_pca.plot$scof <- iso_host_metadata_csp[rownames(csp_m),
                                           "Cold_shock_protein|ScoF"]
csp_pca.plot$csp7 <- iso_host_metadata_csp[rownames(csp_m),
                                           "Cold_shock_protein|7"]
csp_pca.plot[is.na(csp_pca.plot)] <- 0
csp_pca.plot$d <- iso_host_metadata_csp[rownames(csp_m),
                                        "host_vc_R2"]
csp_pca.plot$DT <- iso_host_metadata_csp[rownames(csp_m),
                                         "DoublingTime"]



csp_pca.plot$scof.2 <-
  ifelse(csp_pca.plot$scof >= 1, "p","a")
csp_pca.plot$csp7.2 <-
  ifelse(csp_pca.plot$csp7 >= 1, "p","a")
csp_pca.plot$negative_csp <- "Absence"
csp_pca.plot$negative_csp <- 
  ifelse(csp_pca.plot$scof.2=="p","Scof","Absence")
csp_pca.plot$negative_csp <- 
  ifelse(csp_pca.plot$csp7.2=="p","Csp7",
         csp_pca.plot$negative_csp)
csp_pca.plot$negative_csp <- 
  ifelse(csp_pca.plot$csp7.2=="p" &
           csp_pca.plot$scof.2=="p"
           ,"Both",csp_pca.plot$negative_csp)


csp_pca.plot[csp_pca.plot$negative_csp != "Absence","order"]

csp_pca.plot$order2 <- 
  ifelse(csp_pca.plot$order %in% 
           c(#"Cytophagales",
             "Lactobacillales","Micrococcales"),
         csp_pca.plot$order,"Other")
csp_pca.plot %>% 
  ggplot(aes(x=DCA1,y=DCA2))+
  geom_point(aes(fill=negative_csp),pch=21)+
  mytheme

csp_pca.plot %>%
  #filter(d!=0 & negative_csp=="Both" & log10(DT) <=2) %>%
  filter(d!=0 & log10(DT) <=2 & 
           order == "Micrococcales" &
           csp7.2=="p") %>%
  #filter(d!=0 & log10(DT) <=2 & csp7.2=="p") %>%
  #group_by(order == "Streptomycetales") %>%
  #group_by(negative_csp) %>%
  do(res=cal_lm_p(.$d,log10(.$DT))) -> tmp

csp_pca.plot %>%
  filter(d!=0 & log10(DT) <=2 & 
          csp7.2=="p"
           ) -> csp_pca.plot.csp7
csp_pca.plot.csp7.all <- csp_pca.plot.csp7
csp_pca.plot.csp7.all$order <- "All"
csp_pca.plot.csp7.all$sig <- ".sig"
csp_pca.plot.csp7 %>% 
  filter(order %in% c("Streptomycetales",
                      "Streptosporangiales",
             "Rhizobiales","Micrococcales")) -> csp_pca.plot.csp7.abund
csp_pca.plot.csp7.abund$sig <-"sig"

#rbind.data.frame(csp_pca.plot.csp7.all,csp_pca.plot.csp7.abund) %>%
#  filter(order != "All") %>%
csp_pca.plot.csp7.abund %>%
  ggplot(aes(x=log10(1/DT),y=d))+
  geom_point(colour="black",alpha=0.5)+
  geom_smooth(se=F,linetype="dashed",
              method="lm",colour="#FF0033")+
  facet_wrap(~order)+mytheme +
  labs(y=expression("Specialization"*~italic("d'")),
       x=expression(Log[10]~"host growth rate (doublings/day)"))+
  guides(linetype=F) -> p_csp_7
ggsave("../../04.Figure/01.Raw/p_csp_7.tiff",
       plot=p_csp_7,device = "tiff",
       compression="lzw",
       width = 6.6,height = 4.8,units = "in")
ggsave("../../04.Figure/01.Raw/p_csp_7.pdf",
       plot=p_csp_7,device = "pdf",
       width = 6.6,height = 4.8,units = "in")



csp_pca.plot %>%
  #filter(d!=0 & negative_csp=="Both" & log10(DT) <=2) %>%
  filter(d!=0 & log10(DT) <=2 & 
           order == "Micrococcales" 
  #       & negative_csp == "Scof"
           ) -> all_Micrococcales

csp_pca.plot %>%
  #filter(d!=0 & negative_csp=="Both" & log10(DT) <=2) %>%
  filter(d!=0 & log10(DT) <=2 & 
           order == "Micrococcales" 
                & negative_csp == "Scof"
  ) -> Micrococcales_scof

csp_pca.plot %>%
  #filter(d!=0 & negative_csp=="Both" & log10(DT) <=2) %>%
  filter(d!=0 & log10(DT) <=2 & 
           order == "Micrococcales" 
         & negative_csp == "Csp7"
  ) -> Micrococcales_Csp7

all_Micrococcales$type="Micrococcales"
Micrococcales_scof$type="Micrococcales with ScoF"
Micrococcales_Csp7$type="Micrococcales with Csp7"
all_Micrococcales$sig <- ".sig"
Micrococcales_scof$sig <- ".sig"
Micrococcales_Csp7$sig <- "ns"

rbind.data.frame(all_Micrococcales,
                 Micrococcales_scof,Micrococcales_Csp7) %>%
  ggplot(aes(x=log10(DT),y=d))+
  geom_point(colour="black",alpha=0.5)+
  geom_smooth(aes(linetype=sig),method="lm",
              se=F,colour="#FF0033")+
  mytheme+
  facet_wrap(~type)+
  labs(y=expression(Specialization*~italic("d'")),
       x=expression(log[10]*italic(DT)))+
  guides(linetype=F) -> p_Micrococcales
ggsave("../../04.Figure/01.Raw/p_Micrococcales.tiff",
       plot=p_Micrococcales,device = "tiff",
       compression="lzw",
       width = 6.6,height = 4.8*0.6,units = "in")

csp_pca.plot %>%
  #filter(d!=0 & negative_csp=="Both" & log10(DT) <=2) %>%
  filter(d!=0 & log10(DT) <=2 & 
           order == "Micrococcales" & 
           negative_csp == "Scof"
  ) %>%
  #do(res=cal_lm_p(.$d,log10(.$DT))) -> tmp
  do(res=cal_cor(.$d,log10(.$DT))) -> tmp
  ggplot(aes(x=log10(DT),y=d))+
  geom_point(colour="black",alpha=0.5)+
    geom_smooth(method="lm",
                colour="#FF0033",se=F)+mytheme


csp_pca.plot$csp7 <-
  iso_host_metadata_csp$`Cold_shock_protein|7`
csp_pca.plot$csp7 <- 
  ifelse(is.na(csp_pca.plot$csp7),0,csp_pca.plot$csp7)
csp_pca.plot$scof <-
  iso_host_metadata_csp$`Cold_shock_protein|ScoF`
csp_pca.plot$scof <- 
  ifelse(is.na(csp_pca.plot$scof),0,csp_pca.plot$scof)

csp_pca.plot$csp7.2 <- 
  ifelse(csp_pca.plot$csp7 !=0,1,0)
csp_pca.plot$scof.2 <- 
  ifelse(csp_pca.plot$scof !=0,1,0)
csp_pca.plot$csp7_scof.2 <- 
  ifelse(csp_pca.plot$scof !=0 &
           csp_pca.plot$csp7 !=0,1,0)

csp_pca.plot %>% group_by(order) %>%
  do(data.frame(S=sum(.$scof.2))) %>%
  arrange(-1*S) %>% 
  filter(S > 10) %>%
  ggplot(aes(x=order,y=S))+
  geom_col()+
  coord_flip()+
  labs(y="Copy number of scof")+
  mytheme

csp_pca.plot %>% group_by(order) %>%
  do(data.frame(S=sum(.$csp7_scof.2))) %>%
  arrange(-1*S)

csp_pca.plot %>% group_by(order) %>%
  do(data.frame(S=sum(.$csp7.2))) %>%
  arrange(-1*S) %>% 
  filter(S > 5) %>%
  ggplot(aes(x=order,y=S))+
  geom_col()+
  coord_flip()+
  labs(y="Copy number of csp7")+
  mytheme

#scof#
csp_pca.plot %>% filter(scof !=0 & !is.na(d)) -> tmp
used.orders <- names(sort(-1*table(tmp$order)))[1:11]
csp_pca.plot %>% filter(scof !=0 & !is.na(d) &
                          order %in% used.orders) -> csp_pca.plot2
csp_pca.plot2 %>% 
  group_by(order) %>% do(lm_res=cal_lm_p(.$d,log10(.$DT))) -> lm.scof
lm.scof <- cbind.data.frame(order=lm.scof$order,ldply(lm.scof$lm_res))
lm.scof$sig <- "ns"
lm.scof$sig <- ifelse(lm.scof$P < 0.05,".sig",lm.scof$sig)

csp_pca.plot2 <- merge.data.frame(csp_pca.plot2,lm.scof[,c(1,7)])

csp_pca.plot2 %>%
  filter(order %in% used.orders &
           log10(DT) < 0.5) %>%
  ggplot(aes(x=log10(1/DT),y=d))+
  geom_point(colour="black",alpha=0.5)+
  geom_smooth(aes(linetype=sig),method="lm",
              se=F,colour="#FF0033")+
  mytheme+
  facet_wrap(~order)+
  labs(y=expression("Specialization"*~italic("d'")),
       x=expression(Log[10]~"host growth rate (doublings/day)"))+
  guides(linetype=F) -> p_scof_order

ggsave("../../04.Figure/01.Raw/p_scof_order.tiff",
       plot=p_scof_order,device = "tiff",
       compression="lzw",
       width = 6.6*1.5,height = 4.8*1.5,units = "in")
ggsave("../../04.Figure/01.Raw/p_scof_order.pdf",
       plot=p_scof_order,device = "pdf",
       width = 6.6*1.5,height = 4.8*1.5,units = "in")


iso_host_metadata_hsp %>% 
  ggplot(aes(x=num_VC,y=host_vc_R2))+ 
  geom_point(alpha=0.3)+ 
  geom_smooth(method="lm",linetype="dashed")+ mytheme +
  labs(x="Number of viral clusters",
       y=expression("Host specialization"*~italic("d'")))+
  scale_x_continuous(breaks = seq(1,12)) -> Fig.VC_d 

ggsave("../../04.Figure/01.Raw/Fig.VC_d.tiff",
       plot=Fig.VC_d,device = "tiff",
       compression="lzw",
       width = 6.6,height = 4.8,units = "in")

nVC_R2 <- rowSums(host_vc_M$iso_host_vc_R2)
d_nVC <- data.frame(Host=names(nVC_R2),nVC_R2=nVC_R2)

Genome_traits <- c("Size_Mb","GC","Genes",
                   "num_5S","num_16S","num_23S",
                   "OGT","DoublingTime","host_vc_R2","num_VC")

sel_df <- iso_host_metadata_csp[,Genome_traits]
rownames(sel_df) <- iso_host_metadata_csp$Assembly.Accession
sel_df <- merge.data.frame(sel_df,d_nVC,
                               by.x=0,by.y="Host")


sel_df <-  apply(sel_df[,c(2:10,12)], 2,as.numeric)
sel_df <- as.data.frame(na.omit(sel_df))
sel_df %>% filter(log10(DoublingTime) < 0.8) -> sel_df
sel_df$logDT <- log10(1/sel_df$DoublingTime)

sel_df$`1/NVC` <- 1/sel_df$nVC_R2

sel_df2 <- sel_df[,c(9,1:6,12,11,7)]
colnames(sel_df2) <- c("d'","GS (Mb)","GC %","CDS","5S rRNA","16S rRNA",
                       "23S rRNA","1/NVC","Log10 Gr","OGT")

my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(alpha=0.3) + 
    geom_smooth(method=method, ...)
  p
}

ggpairs(sel_df2,
        lower = list(continuous = my_fn))+
  mytheme -> Fig.varCor
ggsave("../../04.Figure/01.Raw/Fig.varCor.tiff",
       plot=Fig.varCor,device = "tiff",
       compression="lzw",
       width = 6.6*2.4,height = 4.8*2.4,units = "in")
ggsave("../../04.Figure/01.Raw/Fig.varCor.pdf",
       plot=Fig.varCor,device = "pdf",
       #compression="lzw",
       width = 6.6*2.4,height = 4.8*2.4,units = "in")
#Virus size####
used.VC <- colnames(host_vc_M$iso_host_vc_R2)
iso_vpc_final %>%
  filter(VC %in% used.VC) %>%
  select(virus,VC) -> VC_prophage
VC_prophage <- VC_prophage[!duplicated(VC_prophage),]
Virus_size <- 
  read.csv("../../02.Data/01.Processing/VirusMeta/iso_virus_seq_features.csv",
           row.names = 1,header=T)
Virus_size$VC <- str_remove(rownames(Virus_size),".fasta")
Virus_size2 <- merge.data.frame(Virus_size,VC_prophage,
                               by.x="VC",by.y="virus")

#Virus stat####
iso_virus_meta2 %>% filter(Size > 0.003) -> iso_virus_meta2_l3
iso_virus_meta2_l3[!duplicated(iso_virus_meta2_l3$VC),] -> tmp

#an
iso_host_metadata_csp %>% 
  filter(num_csp !=0) -> csp_tmp0
g <- ifelse(csp_tmp0$phylum == "Actinobacteria",
            "Actinobacteria","Other")
summary(aov(csp_tmp0$`Cold_shock_protein|ScoF`~g))
csp_tmp0 %>% group_by(phylum) %>%
  dplyr::summarise(S=sum(`Cold_shock_protein|ScoF`,na.rm = T)) %>%
  arrange(-1*S)
csp_tmp0 %>% filter(`Cold_shock_protein|ScoF` > 0 &
                      !is.na(host_vc_R2)
                      ) %>%
  dplyr::count(phylum)

csp_tmp0 %>% group_by(phylum) %>%
  dplyr::summarise(S=sum(`Cold_shock_protein|7`,na.rm = T)) %>%
  arrange(-1*S)
csp_tmp0 %>% filter(`Cold_shock_protein|7` > 0 &
                      !is.na(host_vc_R2)
) %>%
  dplyr::count(phylum)

#16.Virus metadata####
save("iso_vc_metadata_cont",
     "iso_vc_metadata_disc",
     file="../../02.Data/02.Rdata/iso_vc_metadata.Rdata")

iso_virus_seq <- read.csv("../../02.Data/01.Processing/VirusMeta/iso_virus_seq_features.csv",
                          stringsAsFactors = F)
iso_virus_seq$X <- str_remove(iso_virus_seq$X,".fasta")
iso_vc <- iso_vpc_final[!duplicated(iso_vpc_final[,c("virus","VC")]),
                        c("virus","VC")]
host_vc_M_R2 <- host_vc_M$iso_host_vc_R2
host_vc_M_R2[host_vc_M_R2 != 0] <- 1
vc_hostrange <- as.data.frame(colSums(host_vc_M_R2))
colnames(vc_hostrange) <- "num_host"

iso_vc_medata <- merge.data.frame(iso_vc_metadata_cont,
                 data.frame(dIndex$host_vc_d_vc$iso_host_vc_R2),
                 by.x="VC",by.y=0)
colnames(iso_vc_medata)[5] <- "d"
iso_vc_medata2 <- merge.data.frame(iso_virus_seq,
                                   iso_vc,
                                  by.x="X",by.y="virus")
iso_vc_medata2 <- merge.data.frame(iso_vc_medata2,
                                   data.frame(dIndex$host_vc_d_vc$iso_host_vc_R2),
                                   by.x="VC",by.y=0)
colnames(iso_vc_medata2)[6] <- "d"
iso_vc_medata2 <- merge.data.frame(iso_vc_medata2,
                                   vc_hostrange,
                                   by.x="VC",by.y=0
                                   )

iso_vc_medata <- merge.data.frame(iso_vc_medata,
                                  iso_vc_metadata_disc[,-6],
                                  by="VC")
iso_vc_medata <- merge.data.frame(iso_vc_medata,
                                   vc_hostrange,
                                   by.x="VC",by.y=0
)

host_vc_link <- iso_vpc_final[!duplicated(iso_vpc_final[,c("host","VC")]),
                              c("host","VC")]
iso_vc_medata <- merge.data.frame(iso_vc_medata,host_vc_link,by="VC")
iso_vc_medata <- merge.data.frame(iso_vc_medata,
                                  iso_host_metadata,
                                  by.x="host",
                                  by.y="Assembly.Accession")
iso_vc_medata$OGT_cut40 <- 
  ifelse(iso_vc_medata$OGT >= 40,"y","n")

iso_vc_medata %>% 
  filter(empo_3 != "Plant_Secretion") %>% 
  group_by(empo_3) %>% 
  do(lm_res = cal_lm_p(.$d,log(.$Mean.Size))) -> tmp
tmp$P <- ldply(tmp$lm_res)$P
tmp$sig <- ifelse(tmp$P< 0.05,"y","n")

iso_vc_medata <- merge.data.frame(
  iso_vc_medata,
  tmp,
  by="empo_3")
iso_vc_medata$sig <- factor(iso_vc_medata$sig,levels = c("y","n"))
iso_vc_medata %>% ggplot(aes(x=log(Mean.Size),y=d,linetype=sig))+
  geom_point()+geom_smooth(method="lm")+
  facet_wrap(~empo_3)+guides(linetype=F)+
  mytheme -> p_d_size.vc
ggsave("../../04.Figure/01.Raw/p_d_size.vc.tiff",
       plot=p_d_size.vc,device = "tiff",compression="lzw",
       width = 6.2*2,height = 4.8*2,units = "in")

save("iso_host_metadata",file="../../02.Data/02.Rdata/iso_host_metadata.Rdata")

#Overall fitting####
iso_host_metadata_cas %>% 
  filter( 
    DoublingTime < 1 &
    host_vc_R2 > 0 &
      (
        `Heat shock protein|HSP20` > 0 |
        `Heat shock protein|HSP40` > 0 |
        `Heat shock protein|HSP70` > 0 |
        `Heat shock protein|HSP100` > 0
      ) &
    (
     `Cold_shock_protein|CspC` > 1 |
     `Cold_shock_protein|CspD` > 1 |
     `Cold_shock_protein|CspLA` > 1
       ) &
      (
        `Prot_Recep|fliC` > 0 |
        `Prot_Recep|fljB` > 0 |
        `Prot_Recep|mshA` > 0
      ) &
      (
        (CI > 0 & Cro > 0) | fpsR > 0
      )
          ) -> IC1

IC1 %>%
  ggplot(aes(x=log10(1/IC1$DoublingTime),y=host_vc_R2))+
  geom_point()+
  geom_smooth(method="lm")+
  mytheme +
  labs(y=expression("Specialization"*~italic("d'")),
       x=expression(Log[10]~"host growth rate (doublings/day)")) -> Figure_S17
ggsave("../../04.Figure/01.Raw/Figure_S17.tiff",
       plot= Figure_S17,device = "tiff",compression="lzw",
       width = 6.2,height = 4.8,units = "in")
