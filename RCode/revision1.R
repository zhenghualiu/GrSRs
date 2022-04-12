setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(lme4)
library(dplyr)
library(ggeffects)
library(ggplot2)
RMSE <- function(x){mean(x^2)^0.5}
MAE <- function(x){mean(abs(x))}
load("../../04.Figure/01.Raw/figure1_df.Rdata")
f1_df %>% filter(log10(DoublingTime) < 1) -> f1_df
f1_df$log10Gr <- log10(1/f1_df$DoublingTime)
lm.fit <- lm(d~log10(1/DoublingTime),data = f1_df)
lmer.fit1 <- lmer(d~log10(1/DoublingTime)+
                    (1+log10(1/DoublingTime)|Discrete_Var_value),
                  data = f1_df)
lmer.fit2 <- lmer(d~(1+log10(1/DoublingTime)|Discrete_Var_value),data = f1_df)


anova(lmer.fit1,lmer.fit2,test="LRT")

lmer.fit3 <- lmerTest::lmer(d~log10Gr+(1+log10Gr|Discrete_Var_value),data = f1_df)
pred.mm <- ggpredict(lmer.fit3,terms=c("log10Gr")) 
ggplot(pred.mm) + 
  geom_point(data = f1_df,aes(x = log10(1/DoublingTime), 
                              y = d, colour = Discrete_Var_value),position = "jitter") + 
  geom_abline(slope=-6.079e-02,intercept=6.219e-01,size=0.8,colour="red") +          # slope
#  scale_x_continuous(limits = c(0,3))+
  geom_smooth(data = f1_df,aes(x = log10(1/DoublingTime), 
                               y = d),method = "lm",size=0.8)+
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "red", alpha = 0.3) + # error band
  theme_bw()+
  guides(colour=guide_legend(title=""))+
  mytheme+
  labs(x=expression(Log[10]~"host growth rate (doublings/day)"),
       y="Virus specificity of host d'") -> mix_effect_model

ggsave("../../04.Figure/02.Publish/p_mix_effect_model.tiff",
       device = "tiff",plot = mix_effect_model,
       width = 6.6,height = 4.8,units = "in",dpi=300,
       compression = "lzw")  

fit.s <- summary(lmer.fit3)
y.pre <- predict(lmer.fit3)
r2 <- 1- sum((f1_df$d-y.pre)^2)/sum((f1_df$d-mean(f1_df$d))^2)
Adj_R <-1- sum((f1_df$d-y.pre)^2)/sum((f1_df$d-mean(f1_df$d))^2)*(nrow(f1_df)-1)/(nrow(f1_df)-2-1)
fix_effct <- data.frame(Discrete_Var_value="All",
                        estimate=fit.s$coefficients[2,1],
                        std=fit.s$coefficients[2,2],
                        Adj_R=Adj_R,
                        p=fit.s$coefficients[2,5])
cal_lm_est <- function(y,x){
  lm.fit <- lm(y~x)
  fit.s <- summary(lm.fit)
  df <- data.frame(
                   estimate=fit.s$coefficients[2,1],
                   std=fit.s$coefficients[2,2],
                   Adj_R=fit.s$adj.r.squared,
                   p=fit.s$coefficients[2,4])
  return(df)
}

f1_df %>%
  group_by(Discrete_Var_value) %>%
  do(data.frame(fit.m=cal_lm_est(.$d,log10(1/.$DoublingTime)))) %>%
  as.data.frame() -> fix_effct.lm
colnames(fix_effct.lm) <- colnames(fix_effct)
fix_effct.com <-
  rbind.data.frame(fix_effct.lm,fix_effct)
fix_effct.com$Discrete_Var_value <-
  factor(fix_effct.com$Discrete_Var_value,
         levels = c(rev(levels(fix_effct.lm$Discrete_Var_value)),"All"))
fix_effct.com$sig <-
  ifelse(fix_effct.com$p >= 0.05,"ns","sig")
fix_effct.com$sig <- factor(fix_effct.com$sig,levels = c("ns","sig"))
fix_effct.com$Adj_R[fix_effct.com$Adj_R < 0] <- 0
fix_effct.com %>%
  ggplot(aes(x=Discrete_Var_value,y=estimate))+
  geom_errorbar( 
    aes(ymin = estimate+std, ymax = estimate-std),
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
  scale_fill_manual(values=c("white","black")) -> p_slop_df3.habitat.mix

ggsave("../../04.Figure/01.Raw/p_slop_df3.habitat.mix.tiff",
       plot=p_slop_df3.habitat.mix,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in")

load("../../02.Data/02.Rdata/iso_host_stat_cont.Rdata")
iso_host_stat_cont %>%
  filter(log10(iso_host_stat_cont$DoublingTime) < 5) -> tmp

summary(lm(host_vc_R2~log10(1/DoublingTime)+OGT+log10(1/DoublingTime)*OGT,data=tmp))
summary(lm(host_vc_R2~
             log10(1/DoublingTime)+
             OGT+
             log10(1/DoublingTime)*OGT*Size+
             Size,data=tmp)) -> lm.fit
write.csv(lm.fit$coefficients,
          "../../05.Table/01.Raw/MixEffectModel_DT_OGT_Size.csv",
          quote = F)
iso_host_stat_cont %>%
  filter(OGT >= 40) -> OGT40

summary(lm(host_vc_R2~
             log10(1/DoublingTime)+
             Size+
             OGT+
             log10(1/DoublingTime)*OGT*Size,data=OGT40)) -> lm.fit
write.csv(lm.fit$coefficients,
          "../../05.Table/01.Raw/MixEffectModel_DT_OGT_Size_40.csv",
          quote = F)

summary(lm(host_vc_R2~
             poly(log10(1/DoublingTime),2)*OGT*Size,
           data=OGT40)) -> poly40.lm.fit
write.csv(poly40.lm.fit$coefficients,
          "../../05.Table/01.Raw/MixEffectModel_DT_OGT_Size_40_poly.lm.fit.csv",
          quote = F)
summary(lm(host_vc_R2~
             poly(log10(1/DoublingTime),2)*OGT*Size,
           data=iso_host_stat_cont)) -> poly.lm.fit
write.csv(poly.lm.fit$coefficients,
          "../../05.Table/01.Raw/MixEffectModel_DT_OGT_Size_all_poly.lm.fit.csv",
          quote = F)
#phylogenetic regression####
library(ape)
library(phylolm)
library(phytools)
tre <- read.tree("../../02.Data/01.Processing/GToTree.tre")
#run 3. statistical analyses L1-254 first
tre_n0 <- tre
el <- sum(tre_n0$edge.length==0)
for(i in tre_n0$tip.label)
{
    tre.tmp <- drop.tip(tre_n0,i)
    el.tmp <- sum(tre.tmp$edge.length==0)
    if(el.tmp < el) {tre_n0 <- tre.tmp; el <- el.tmp}
}
  
iso_host_metadata_cont.metadata  <- 
  merge.data.frame(
                   iso_host_metadata[,-c(114,10,118)],
                   iso_host_metadata_cont,
                   by.x="Assembly.Accession",
                   by.y="GenomeID")
rownames(iso_host_metadata_cont.metadata) <- 
  iso_host_metadata_cont.metadata$Assembly.Accession
iso_host_metadata_cont.metadata %>% 
  filter(Assembly.Accession %in% tre_n0$tip.label) -> iso_host_metadata_cont.metadata.phylm
iso_host_metadata_cont.metadata.phylm <-
  merge.data.frame(iso_host_metadata_cont.metadata.phylm,
                   iso_d_host,
                   by.x="Assembly.Accession",
                   by.y=0)
iso_host_metadata_cont.metadata.phylm %>%
  filter(!is.na(host_vc_R2)) -> tmp
rownames(tmp) <- tmp$Assembly.Accession
{
#BM model
phylolm.fit.BM = phylolm(
  host_vc_R2~log10(1/DoublingTime),
  model = "BM",
  #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
               phy=tre_n0,
               data=tmp,
 # measurement_error=TRUE,
               boot=100)
RMSE(phylolm.fit.BM$residuals)
MAE(phylolm.fit.BM$residuals)
mean(phylolm.fit.BM$residuals)
phylolm.fit.s.BM <- summary(phylolm.fit.BM)
write.csv(phylolm.fit.s.BM$coefficients,
          "../../02.Data/01.Processing/revision1/phylolm.fit.s.BM.csv",quote = F)
#EB model
phylolm.fit.EB = phylolm(
  host_vc_R2~log10(1/DoublingTime),
  model = "EB",
  #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
  phy=tre_n0,
  data=tmp,
  boot=100#,lower.bound = -1, upper.bound = 5
  )
phylolm.fit.s.EB <- summary(phylolm.fit.EB)
write.csv(phylolm.fit.s.EB$coefficients,
          "../../02.Data/01.Processing/revision1/phylolm.fit.s.EB.csv",quote = F)
#lambda model
phylolm.fit.lambda = phylolm(
  host_vc_R2~log10(1/DoublingTime),
  model = "lambda",
  #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
  phy=tre_n0,
  data=tmp,
  boot=100#,lower.bound = -1, upper.bound = 5
)
phylolm.fit.s.lambda <- summary(phylolm.fit.lambda)
write.csv(phylolm.fit.s.lambda$coefficients,
          "../../02.Data/01.Processing/revision1/phylolm.fit.s.lambda.csv",quote = F)
#kappa model
phylolm.fit.kappa = phylolm(
  host_vc_R2~log10(1/DoublingTime),
  model = "kappa",
  #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
  phy=tre_n0,
  data=tmp,
  boot=100#,lower.bound = -1, upper.bound = 5
)
phylolm.fit.s.kappa <- summary(phylolm.fit.kappa)
write.csv(phylolm.fit.s.kappa$coefficients,
          "../../02.Data/01.Processing/revision1/phylolm.fit.s.kappa.csv",quote = F)
#delta model
phylolm.fit.delta = phylolm(
  host_vc_R2~log10(1/DoublingTime),
  model = "delta",
  #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
  phy=tre_n0,
  data=tmp,
  boot=100#,lower.bound = -1, upper.bound = 5
)
phylolm.fit.s.delta <- summary(phylolm.fit.delta)
write.csv(phylolm.fit.s.delta$coefficients,
          "../../02.Data/01.Processing/revision1/phylolm.fit.s.delta.csv",quote = F)
}#phylolm

lm.fit = lm(
  host_vc_R2~log10(1/DoublingTime),
  # host_vc_R2~log10(1/DoublingTime)+OGT+Size,
                      data=tmp)
lm.fit.s <- summary(lm.fit)
write.csv(lm.fit.s$coefficients,
          "../../02.Data/01.Processing/revision1/lm.fit.s.csv",quote = F)

RMSE(lm.fit$residuals)
MAE(lm.fit$residuals)
mean(lm.fit$residuals)
#phylosig test
#error
log10Gr <- log10(1/tmp$DoublingTime)
names(log10Gr) <- tmp$Assembly.Accession

e_physig_K <- phylosig(tre_n0, lm.fit$residuals, method="K", test = TRUE,nsim=1000)
e_physig_lambda <- phylosig(tre_n0, lm.fit$residuals, method="lambda", test = TRUE,nsim=1000)

log10Gr_physig_K <- phylosig(tre_n0, log10Gr, method="K", test = TRUE,nsim=1000)
log10Gr_physig_lambda <- phylosig(tre_n0, log10Gr, method="lambda", test = TRUE,nsim=1000)


tmp$phylum[tmp$phylum=="delta/epsilon subdivisions"] <- "delta_epsilon_subdivisions"
phylum.count <- table(tmp$phylum)
phylum.n <- names(phylum.count)[phylum.count >= 5]

output_path <- "../../02.Data/01.Processing/revision1/"
BM.fit <- list()
k <- 1
for(i in phylum.n)
{
  tmp %>% filter(phylum == i) -> dta.tmp
  phylolm.fit.tmp = phylolm(
    #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
    host_vc_R2~log10(1/DoublingTime),
                        phy=tre_n0,
                        model = "BM",
                        data=dta.tmp,
                        boot=100)
  BM.fit[[k]] <- phylolm.fit.tmp
  k <- k + 1
  phylolm.fit.tmp.s <- summary(phylolm.fit.tmp)
  file.tmp <- paste(output_path,i,"phylolm.fit.s.BM.csv",sep="")
  out_con <- file(file.tmp,"w")
  write(t(phylolm.fit.tmp.s$coefficients),
        out_con,
        ncolumns=6,
        append = T,sep = ",")
  write(paste("Adjusted R-squared: ",phylolm.fit.tmp$adj.r.squared,sep=""),out_con,append = T)
  write(paste("RMSE: ",RMSE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("MAE: ",MAE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("Centre of resudials: ",mean(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  close(out_con)
}
names(BM.fit) <- phylum.n
delta.fit <- list()
k <- 1
for(i in phylum.n)
{
  tmp %>% filter(phylum == i) -> dta.tmp
  phylolm.fit.tmp = phylolm(
    #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
    host_vc_R2~log10(1/DoublingTime),
    phy=tre_n0,
    model = "delta",
    data=dta.tmp,
    boot=100)
  delta.fit[[k]] <- phylolm.fit.tmp
  k <- k + 1
  phylolm.fit.tmp.s <- summary(phylolm.fit.tmp)
  file.tmp <- paste(output_path,i,"phylolm.fit.s.delta.csv",sep="")
  out_con <- file(file.tmp,"w")
  write(t(phylolm.fit.tmp.s$coefficients),
        out_con,
        ncolumns=6,
        append = T,sep = ",")
  write(paste("Adjusted R-squared: ",phylolm.fit.tmp$adj.r.squared,sep=""),out_con,append = T)
  write(paste("RMSE: ",RMSE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("MAE: ",MAE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("Centre of resudials: ",mean(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  close(out_con)
}
names(delta.fit) <- phylum.n
for(i in phylum.n)
{
  tmp %>% filter(phylum == i) -> dta.tmp
  phylolm.fit.tmp = phylolm(
    #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
    host_vc_R2~log10(1/DoublingTime),
    phy=tre_n0,
    model = "EB",
    data=dta.tmp,
    boot=100)
  phylolm.fit.tmp.s <- summary(phylolm.fit.tmp)
  file.tmp <- paste(output_path,i,"phylolm.fit.s.EB.csv",sep="")
  out_con <- file(file.tmp,"w")
  write(t(phylolm.fit.tmp.s$coefficients),
        out_con,
        ncolumns=6,
        append = T,sep = ",")
  write(paste("Adjusted R-squared: ",phylolm.fit.tmp$adj.r.squared,sep=""),out_con,append = T)
  write(paste("RMSE: ",RMSE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("MAE: ",MAE(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  write(paste("Centre of resudials: ",mean(phylolm.fit.tmp$residuals),sep=""),out_con,append = T)
  close(out_con)
}



phylum.count <- table(tmp$phylum)
phylum.n <- names(phylum.count)[phylum.count >= 8]

output_path <- "../../02.Data/01.Processing/revision1/"
for(i in phylum.n)
{
  tmp %>% filter(phylum == i) -> dta.tmp
  lm.fit.tmp = lm(
    #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
    host_vc_R2~log10(1/DoublingTime),
    data=dta.tmp)
  lm.fit.tmp.s <- summary(lm.fit.tmp)
  file.tmp <- paste(output_path,i,"lm.fit.s.csv",sep="")
  out_con <- file(file.tmp,"w")
  write(t(lm.fit.tmp.s$coefficients),
        out_con,
        ncolumns=4,
        append = T,sep = ",")
  write(paste("Adjusted R-squared: ",lm.fit.tmp.s$adj.r.squared,sep=""),out_con,append = T)
  write(paste("RMSE: ",RMSE(lm.fit.tmp.s$residuals),sep=""),out_con,append = T)
  write(paste("MAE: ",MAE(lm.fit.tmp.s$residuals),sep=""),out_con,append = T)
  write(paste("Centre of resudials: ",mean(lm.fit.tmp.s$residuals),sep=""),out_con,append = T)
  close(out_con)
}

#Get P value from table S3 for each phylum
lm.p.val <- scan()
p.adjust(lm.p.val, method = "BH", n = length(lm.p.val))

Discrete_Var_value <- c(
  "Actinobacteria","Bacteroidetes","Cyanobacteria",
  ""
)

output_path <- "../../02.Data/01.Processing/revision1/"
lm.size.fit <- list()
for(i in 1:length(phylum.n))
{
  tmp %>% filter(phylum == phylum.n[i]) -> dta.tmp
  lm.fit.tmp = lm(
    host_vc_R2~log10(1/DoublingTime)+Size,
    #host_vc_R2~log10(1/DoublingTime),
    data=dta.tmp)
  lm.fit.tmp.s <- summary(lm.fit.tmp)
  df <- as.data.frame(lm.fit.tmp.s$coefficients)
  df$adj.r.squared <- lm.fit.tmp.s$adj.r.squared
  df$term <- rownames(df)
  lm.size.fit[[i]] <- df
}
names(lm.size.fit) <- phylum.n

lm.size.fit.df <- plyr::ldply(lm.size.fit)
lm.size.fit.df %>% 
  filter(term != "(Intercept)") -> lm.size.fit.df
write.csv(lm.size.fit.df,
          "../../02.Data/01.Processing/revision1/lm.size.fit.df.csv")


phylolm.size.fit <- list()
for(i in phylum.n)
{
  tmp %>% filter(phylum == i) -> dta.tmp
  phylolm.fit.tmp = phylolm(
    #host_vc_R2~log10(1/DoublingTime)+OGT+Size,
    host_vc_R2~log10(1/DoublingTime)+Size,
    model="BM",
   # measurement_error=TRUE,
    phy=tre_n0,
    data=dta.tmp,
    boot=100)
  phylolm.fit.tmp.s <- summary(phylolm.fit.tmp)
  df <- as.data.frame(phylolm.fit.tmp.s$coefficients)
  df$adj.r.squared <- phylolm.fit.tmp.s$adj.r.squared
  df$term <- rownames(df)
  phylolm.size.fit[[i]] <- df
}


for(i in 1:length(phylum.n))
{
  tmp %>% filter(phylum == phylum.n[i]) -> dta.tmp
  lm.fit.tmp = lm(
    host_vc_R2~log10(1/DoublingTime)+Size,
    #host_vc_R2~log10(1/DoublingTime),
    data=dta.tmp)
  lm.fit.tmp.s <- summary(lm.fit.tmp)
  df <- as.data.frame(lm.fit.tmp.s$coefficients)
  df$adj.r.squared <- lm.fit.tmp.s$adj.r.squared
  df$term <- rownames(df)
  lm.size.fit[[i]] <- df
}
names(lm.size.fit) <- phylum.n


#virus information#
library(ggplot2)
library(brew)
library(scales)
iso_gbg <- read.csv("../../02.Data/01.Processing/VirusCluster/iso_genome_by_genome_overview.csv")
iso_gbg %>% 
  filter(Order != "Unassigned") %>% 
  filter(!duplicated(VC.Subcluster))-> iso_vc_assigned
iso_vc_assigned <- unique(iso_vc_assigned$VC.Subcluster)
iso_vc.o <- read.csv("../../02.Data/01.Processing/VirusCluster/iso_viral_cluster_overview.csv")
iso_vc <- iso_vc.o[,c("VC","Members")]
iso_v <- str_split(iso_vc$Members,",")
names(iso_v) <- iso_vc$VC
iso_v <- lapply(iso_v,function(x){
  df <- data.frame(Vir=x)
  return(df)
})
iso_v <- plyr::ldply(iso_v,.id = "VC")
iso_v$Vir <- as.character.factor(iso_v$Vir)
iso_v$Vir <- str_remove(iso_v$Vir,".faa")
iso_v_r5 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R5) )
iso_v_r4 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R4) )
iso_v_r3 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R3) )
iso_v_r2 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R2) )
iso_v_r1 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R1) ) 
iso_v_r0 <- iso_v %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R0) ) 
iso_v_remove <- colnames(host_vc_M$iso_host_vc_R0)[
  !colnames(host_vc_M$iso_host_vc_R0) %in% colnames(host_vc_M$iso_host_vc_R2)]
vir_completeness <- 
  read.csv("../../02.Data/01.Processing/completeness.tsv",sep="\t")
quality_summary <-
  read.csv("../../02.Data/01.Processing/quality_summary.tsv",sep="\t")
iso_vpc_final %>% filter(VC %in% colnames(host_vc_M$iso_host_vc_R2)) -> iso_vir_vc
iso_gbg %>% 
  filter(VC.Subcluster %in% colnames(host_vc_M$iso_host_vc_R2)) %>% 
  filter(Order != "Unassigned") -> iso_gbg_assigned_R2
iso_gbg %>% 
  filter(VC.Subcluster %in% colnames(host_vc_M$iso_host_vc_R0)) %>% 
  filter(Order != "Unassigned") -> iso_gbg_assigned_R0
iso_gbg %>% 
#  filter(VC.Subcluster %in% colnames(host_vc_M$iso_host_vc_R0)) %>% 
  filter(Order != "Unassigned") -> iso_gbg_assigned

iso_vir_vc <- unique(iso_vir_vc$virus)
quality_summary %>% filter(contig_id %in% iso_vir_vc)
phage_fna_info <- 
  read.table("../../02.Data/01.Processing/phage_fna_info.txt",sep="\t",header=F)
phage_fna_info$V1 <-
  str_remove(phage_fna_info$V1,".fasta")
phage_fna_info$V2 <-
  str_remove(phage_fna_info$V2,">")


phage_fna_info %>%
  filter(V1 %in% iso_v_r1$Vir) -> phage_fna_info_r1
phage_fna_info %>%
  filter(V1 %in% iso_v_r2$Vir) -> phage_fna_info_r2
phage_fna_info %>%
  filter(V1 %in% iso_v_r3$Vir) -> phage_fna_info_r3
phage_fna_info %>%
  filter(V1 %in% iso_v_r3$Vir) -> phage_fna_info_r4
phage_fna_info %>%
  filter(V1 %in% iso_v_r3$Vir) -> phage_fna_info_r5


phage_fna_info %>%
  filter(V1 %in% iso_v_remove) -> phage_fna_info_remove

quality_summary %>%
  filter(contig_length >= 3000) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> quality_count

quality_count$pro <-
  quality_count$n/sum(quality_count$n)
quality_count$text <- paste(quality_count$checkv_quality," ",
                            percent(quality_count$pro,0.01),
                            " (",
                            quality_count$n,")",sep="")
quality_count$pro <- quality_count$n/sum(quality_count$n)
quality_count$checkv_quality <- 
  factor(quality_count$checkv_quality,
         levels = c("Complete","High-quality","Medium-quality","Low-quality","Not-determined") )
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

quality_count %>%
  arrange(checkv_quality) -> quality_count

quality_count %>%
  ggplot(aes(x="",y=pro,fill=checkv_quality))+
  geom_bar(width=1,stat="identity")+
  coord_polar("y",start=0)+
  blank_theme+
  scale_fill_manual(labels=quality_count$text,
                    values = brewer.pal(5,"Set1"))+
  guides(fill=guide_legend(title=""))+
  theme(axis.text.x = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,vjust=0,size=20))+
  labs(title="All proviruses") -> quality_count.p1

quality_summary %>%
  filter(contig_id %in% phage_fna_info_r1$V2) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> r1_quality

quality_summary %>%
  filter(contig_id %in% phage_fna_info_r2$V2) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> r2_quality

quality_summary %>%
  filter(contig_id %in% phage_fna_info_r3$V2) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> r3_quality

quality_summary %>%
  filter(contig_id %in% phage_fna_info_r4$V2) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> r4_quality

quality_summary %>%
  filter(contig_id %in% phage_fna_info_r5$V2) %>%
  group_by(checkv_quality) %>%
  dplyr::count() -> r5_quality

remove_rate <- 1-unlist(c(r1_quality[5,2],r2_quality[5,2],
           r3_quality[5,2],r4_quality[5,2],
           r5_quality[5,2]))/5803
N_unvai <- unlist(c(r1_quality[5,2],r2_quality[5,2],
             r3_quality[5,2],r4_quality[5,2],
             r5_quality[5,2]))
N_viralclusters <-
  c(
    ncol(host_vc_M$iso_host_vc_R1),
    ncol(host_vc_M$iso_host_vc_R2),
    ncol(host_vc_M$iso_host_vc_R3),
    ncol(host_vc_M$iso_host_vc_R4),
    ncol(host_vc_M$iso_host_vc_R5)
  )
N_vir <- 
  c(
    nrow(phage_fna_info_r1),
    nrow(phage_fna_info_r2),
    nrow(phage_fna_info_r3),
    nrow(phage_fna_info_r4),
    nrow(phage_fna_info_r5)
  )

data.frame(N=seq(1,5),rate=1-remove_rate,VC=N_unvai/N_vir) %>%
  melt.data.frame(id.vars = "N") %>%
  ggplot(aes(x=N,y=value*100,colour=variable))+
  geom_point(size=3)+
  mytheme+
  labs(x="Network degree of removed nodes",
       y="Proportion of network proviruses\n without available completeness  (%)")+
  guides(colour=guide_legend(title=""))+
  scale_colour_manual(values=c("red","blue"),
                      labels=c("In all proviruses",
                              "In host-provirus network"))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=15)) -> remove_proviruses

ggsave("../../04.Figure/02.Publish/remove_proviruses.tiff",
       plot=remove_proviruses,device = "tiff",compression="lzw",
       width = 6.6,height = 4.8,units = "in",scale = 1
       )  

#aov
library(agricolae)
iso_host_metadata_aov <- merge.data.frame(
  iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0)
iso_host_metadata_aov$phylum[
  iso_host_metadata_aov$phylum == "delta/epsilon subdivisions"] <- "Delta/epsilon subdivisions"
  
iso_host_metadata_aov %>% 
  filter(empo_3 != "Unknown"
         & !is.na(host_vc_R2)
         ) -> df_env
summary(aov(host_vc_R2~empo_3,data=df_env))
iso_host_metadata_aov %>% 
  filter(phylum != ""
         & !is.na(host_vc_R2)
  ) -> df_phylum
count(df_phylum$phylum) %>% filter(freq >= 5) -> df_phylum2
df_phylum %>% filter(phylum %in% df_phylum2$x) -> df_phylum3
model <- aov(host_vc_R2~phylum,data=df_phylum3)
out <- LSD.test(model,"phylum", p.adj="bonferroni")
out$means$ID <- rownames(out$means)
out$means$ID <- factor(
  out$means$ID,levels = rownames(out$groups))
out$groups$ID <- rownames(out$groups)
ggplot(data=out$means,
       aes(x=ID,y=host_vc_R2))+
  geom_point()+
  coord_flip()+
  geom_pointrange(aes(ymax=host_vc_R2+std,
                      ymin=host_vc_R2-std))+
  geom_text(data=out$groups,aes(x=ID,y=0.95,label=groups))+
  mytheme+
  scale_y_continuous(limits = c(0.4,1))+
  labs(x="",y=expression("Virus specificity"~italic("d'"))) -> phylum_d



#Major comment 13
library(cowplot)
load("../../02.Data/02.Rdata/iso_host_vc_R0_d.Rdata")
load("../../02.Data/02.Rdata/2.dIndex.Rdata")
iso_host_stat_cont %>%
  filter(!is.na(host_vc_R2)) -> tmp
data_r0_d <- data.frame(d=iso_host_vc_R0_d$dprime,type="Whole network")
data_r2_d <- data.frame(d=tmp$host_vc_R2,type="b")

data_r0_d %>%
  ggplot(aes(x=d))+
  geom_histogram()+
  mytheme +
  labs(x=expression("Specialization"~italic("d'")),
       y="Frequency")  -> p.r0

data_r2_d %>%
  ggplot(aes(x=d))+
  geom_histogram()+
  mytheme +
  labs(x=expression("Specialization"~italic("d'")),
       y="Frequency")  -> p.r2

plot_grid(p.r0,p.r2,labels = c("A","B")) -> p.d_dis
ggsave("../../04.Figure/01.Raw/p.d_dis.tiff",
       plot=p.d_dis,device = "tiff",compression="lzw",
       width = 6.6*2,height = 4.8,units = "in",scale = 0.8)
ggsave("../../04.Figure/01.Raw/p.d_dis.pdf",
       plot=p.d_dis,device = "pdf",
       width = 6.6*2,height = 4.8,
       units = "in",scale = 0.8)

iso_host_metadata.r0 <-
  merge.data.frame(iso_host_metadata_cont.metadata.phylm,
                   data_r0_d,by.x="GenomeID",
                   by.y=0)
iso_host_metadata.r0 %>%
  filter(log10(1/DoublingTime) > -4) -> r0.tmp


#completeness of singletons and doubletons
quality_summary %>% filter(contig_id %in% phage_fna_info_remove$V2) -> single_double_com
quality_summary %>% filter(contig_id %in% phage_fna_info_r2$V2) -> r2_com


iso_host_metadata_d <- merge.data.frame(
  iso_host_metadata,iso_d_host,by.x="Assembly.Accession",by.y = 0)

#adjust P values####
#Receptor: p value from table S5
Rep.p.val <- scan()
Rep.p.val <- na.omit(Rep.p.val)
p.adjust(Rep.p.val, method = "BH", n = length(Rep.p.val))
