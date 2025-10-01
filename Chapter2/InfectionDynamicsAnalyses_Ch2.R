###############################################
##### BAT - Analysis - Infection Dynamics #####
###############################################

setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/5.Thesis/Project/7.RADseq/Project_RAD2024/Bats/Data_Alyzee/overall/Results/InfDyn")

load("./LAST_VERSION_STATS.RData") ## Load last R history


################################################################################
# 1. Download packages and data ####
################################################################################
## 1.1 Download packages ####

library(hierfstat)
library(ade4)
library(reshape2)
library(lme4)
library(car)
library(ggplot2)
library(viridis)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(viridisLite)
library(adegenet)
library(mgcv)
library(glmmTMB)
library(sjPlot)
library(emmeans)
library(performance)
library(randomForest)
library(see)
library(DHARMa)
library(broom)
library(SNPRelate)

## 1.2 Download data ####

BAT_vcf <- read.VCF("../../vcf_after_filtering.recode.vcf")
dim(BAT_vcf) #359 35919


################################################################################
# 2. Creation of the datasets ####
################################################################################
## 2.1 Genomic dataset ####
### 2.1.1 Including host, sites, sex and host sex ####

popmap <- read.csv("../../pop_map_2.txt", header =TRUE, sep="\t")
poplist <- popmap$Site
IDlist <- popmap$ID

### 2.1.2 SNPs dataset ####

genotype_matrix <- as.matrix(BAT_vcf)
BAT <- data.frame(ID=IDlist, pop=poplist, genotype_matrix)

#Sort by population
BAT_popsorted <- BAT[order(BAT$pop),]
BAT_popsorted_matrix <- as.matrix(BAT_popsorted[,3:ncol(BAT_popsorted)])
BAT_popsorted_matrix <- apply(BAT_popsorted_matrix, 2, as.numeric)

## 2.2 Infection dataset ####
### 2.2.1 Infection all individuals ####

BAT_infection_ALL <- read.csv("./Bat_infection_ALL.csv", header=TRUE, sep=";")

### 2.2.2 Dataset infection with 2023 individuals ####

BAT_inf <- BAT_infection_ALL[BAT_infection_ALL$ID %in% BAT_popsorted$ID,]
str(BAT_inf)

### 2.2.3 Clean the dataset ####

#Replace P/NP by NA for ectoparasites
BAT_inf$Spx_tot <- ifelse(BAT_inf$Spx_tot == "P" | BAT_inf$Spx_tot == "NP", NA, BAT_inf$Spx_tot)
BAT_inf$Nyct_tot <- ifelse(BAT_inf$Nyct_tot == "P" | BAT_inf$Nyct_tot == "NP", NA, BAT_inf$Nyct_tot)

#Change factors
BAT_inf$Spx_tot <- as.numeric(BAT_inf$Spx_tot)
BAT_inf$Nyct_tot <- as.numeric(BAT_inf$Nyct_tot)
BAT_inf$AB <- as.numeric(BAT_inf$AB)
BAT_inf$Age <- factor(BAT_inf$Age, levels = c("AD", "SAD", "JUV"))
BAT_inf$Sex <- factor(BAT_inf$Sex)
BAT_inf$meanQuantity_Poly <- as.numeric(gsub(",", ".", BAT_inf$meanQuantity_Poly))

str(BAT_inf)


### 2.2.4 Polchromophilus infections ####

BAT_inf$INF_POLY <- BAT_inf$meanQuantity_Poly/(BAT_inf$meanQuantity_Poly+BAT_inf$meanQuantity_Bat)

BAT_inf$PREV_POLY <- ifelse(BAT_inf$INF_POLY > 0,1,0)


################################################################################
# 3. Genetic Diversity Measures ####
################################################################################
## 3.1 FUNI
### 3.1.1 Function ####

get.funiwn<-function(dos){
  if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
  p<-dos@p
  het<-2*p*(1-p)
  res<-apply(gaston::as.matrix(dos),1,function(x) {
    nas<-which(is.na(x));
    if(length(nas)>0){
      xs<-x[-nas];
      ps<-p[-nas];
      hets<-het[-nas];
      sum(xs^2-(1+2*ps)*xs+2*ps^2)/sum(hets);
    }
    else sum(x^2-(1+2*p)*x+2*p^2)/sum(het);
  }
  )
  return(list(Funi=unlist(res),het=sum(het)))
}


### 3.1.2 FUNI per individuals ####

funi<-get.funiwn(BAT_vcf)
Funi<-funi$Funi

### 3.1.3 Add FUNI in dataset ####

BAT_inf$FUNI <- Funi

### 3.1.4 FUNI per populations ####

Funi_pop <- data.frame(ID=IDlist, pop=poplist, Funi)
boxplot(Funi_pop$Funi~Funi_pop$pop)

Funi_kw <- kruskal.test(Funi_pop$Funi~Funi_pop$pop)#Kruskal-Wallis chi-squared = 21.784, df = 17, p-value = 0.1932

mean_Funi_pop <- aggregate(Funi ~ pop, data = Funi_pop, FUN = mean)

### 3.1.5 Plot ####

svg("./FUNI_NYCT.svg")
plot(BAT_inf$FUNI~BAT_inf$Nyct_tot, col="#E69F00", pch=19,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
axis (side=1, at=c(0, 5, 10, 15, 20))
title(xlab="Abundance bat flies", ylab="FUNI")
dev.off()

svg("./FUNI_SPX.svg")
plot(BAT_inf$FUNI~BAT_inf$Spx_tot, col="#0072B2", pch=19,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
axis (side=1, at=c(0, 10, 20, 30, 40, 50))
title(xlab="Abundance wing mites", ylab="FUNI")
dev.off()

svg("./FUNI_POLY.svg")
plot(BAT_inf$FUNI~BAT_inf$INF_POLY, col="#009E73", pch=19,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
axis (side=1, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
title(xlab="Relative abundance Polychromophilus", ylab="FUNI")
dev.off()


#Who is the outlier?
BAT_inf$ID[BAT_inf$FUNI > 0.2] #515

plot(BAT_vcf@ped$N1,BAT_vcf@ped$N2, pch=20)
which(BAT_vcf@ped$id=="515") #134
points(BAT_vcf@ped$N1[134],BAT_vcf@ped$N2[134],col="#33CCCC",pch=17)
table(as.matrix(BAT_vcf)[134,]) # 0 (26100); 1 (4469); 2 (1773)



## 3.2 FAS ####
### 3.2.1 Estimate FAS ####

kinship <- beta.dosage(genotype_matrix)
FAS <- diag(kinship)

### 3.2.2 Add FAS in dataset ####

BAT_inf$FAS <- FAS

### 3.2.3 Plot ####

plot(BAT_inf$FAS~BAT_inf$Nyct_tot)
plot(BAT_inf$FAS~BAT_inf$Spx_tot)
plot(BAT_inf$FAS~BAT_inf$INF_POLY)

#Who is the outlier?
BAT_inf$ID[BAT_inf$FAS > 0.2] #515
# 
# plot(BAT_vcf@ped$N1,BAT_vcf@ped$N2, pch=20)
# which(BAT_vcf@ped$id=="392") #24
# points(BAT_vcf@ped$N1[24],BAT_vcf@ped$N2[24],col="#33CCCC",pch=17)
# table(as.matrix(BAT_vcf)[24,]) # 0 (27338); 1 (6560); 2 (1118)



## 3.3 Heterozygosity ####
### 3.3.1 Estimate HZ ####

DATA = snpgdsVCF2GDS("../../vcf_after_filtering.recode.vcf", "vcf_after_filtering_BAT.gds")
DATA = snpgdsOpen("vcf_after_filtering_BAT.gds")
gen_mat <- snpgdsGetGeno('vcf_after_filtering_BAT.gds')
IndHet = rowSums(gen_mat==1, na.rm = TRUE) / rowSums(! is.na(gen_mat))

### 3.3.2 Add HZ to the dataset ####

BAT_inf$HZ <- IndHet


### 3.3.3 Plot ####

plot(BAT_inf$HZ~BAT_inf$Nyct_tot)
plot(BAT_inf$HZ~BAT_inf$Spx_tot)
plot(BAT_inf$HZ~BAT_inf$INF_POLY)

#Who is the outlier?
BAT_inf$ID[BAT_inf$HZ < 0.15] #515
#same than FUNI


################################################################################
# 4. Find best models ####
################################################################################

##I did FUNI then tested for FAS and heterozygosity


## 4.1 Bat flies ####
### 4.1.1 Prevalence & Distribution ####

#Prevalence

PrevNyct <- sum(BAT_inf$Nyct_tot>=1)/nrow(BAT_inf)


#Distribution --> Count so poisson family
svg("./Distrib_NYCT_sqrt.svg")
hist(
  sqrt(BAT_inf$Nyct_tot),
  breaks=100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(sqrt(BAT_inf$Nyct_tot), na.rm = TRUE), max(sqrt(BAT_inf$Nyct_tot), na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1)
dev.off()


### 4.1.2 Family selection ####

model_Nnb1 <- glmmTMB(Nyct_tot ~ HZ + Age * Sex + (1|Site), data=BAT_inf, 
                 family=nbinom1)
model_Nnb2 <- glmmTMB(Nyct_tot ~ HZ + Age * Sex + (1|Site), data=BAT_inf,
                  family=nbinom2)
model_Np <- glmmTMB(Nyct_tot ~ HZ + Age * Sex + (1|Site), data=BAT_inf,
                   family=poisson)
model_Ng <- glmmTMB(log(Nyct_tot+1) ~ HZ + Age * Sex + (1|Site), data=BAT_inf,
                    family=gaussian)
model_Ngsqrt <- glmmTMB(sqrt(Nyct_tot) ~ HZ + Age * Sex + (1|Site), data=BAT_inf,
                    family=gaussian)
AIC(model_Nnb1, model_Nnb2, model_Np, model_Ng, model_Ngsqrt)
#FUNI, FAS, HZ: Gaussian + log = best model


## Check normality and homoscedasticity
hist(residuals(model_Ng)) #Ok for normality
plot(predict(model_Ng, type = "response"),resid(model_Ng, type = "pearson")) #OK for homoscedasticity
check_model(model_Ng)
#FUNI, FAs, HZ: ok


### 4.1.3 Model selection ####

##Removing site
model_Ng_site <- glmmTMB(log(Nyct_tot+1) ~ HZ + Age * Sex, data=BAT_inf,
                    family=gaussian)
anova(model_Ng, model_Ng_site)
#FUNI, FAS, HZ: Site significant


##Removing interaction Age * Sex
model_Ng_inter <- glmmTMB(log(Nyct_tot+1) ~ HZ + Age + Sex + (1|Site), data=BAT_inf,
                         family=gaussian)
anova(model_Ng, model_Ng_inter) 
#FUNI, FAS, HZ: Interaction not significant


##Removing Sex
model_Ng_sex <- glmmTMB(log(Nyct_tot+1) ~ HZ + Age + (1|Site), data=BAT_inf,
                          family=gaussian)
anova(model_Ng_inter, model_Ng_sex) 
#FUNI, FAS, HZ: Sex significant


##Removing Age
model_Ng_age <- glmmTMB(log(Nyct_tot+1) ~ HZ + Sex + (1|Site), data=BAT_inf,
                          family=gaussian)
anova(model_Ng_inter, model_Ng_age) 
#FUNI, FAS, HZ: Age not significant


##Removing FUNI, FAS, HZ
model_Ng_FUNI <- glmmTMB(log(Nyct_tot+1) ~ Sex + (1|Site), data=BAT_inf,
                        family=gaussian)
anova(model_Ng_age, model_Ng_FUNI) 
#FUNI, FAS, HZ not significant



### 4.1.4 Best selection ####

model_Ng_FUNI <- glmmTMB(log(Nyct_tot+1) ~ Sex + (1|Site), data=BAT_inf,
                        family=gaussian)
summary(model_Ng_FUNI)#SexM
pairs(emmeans(model_Ng_FUNI,"Sex"), adjust="tukey") #F more parasitized than M

plot_model(model_Ng_FUNI, type = "emm", terms = "Sex", show.data = F)


## 4.2 Wing mites ####
### 4.2.1 Prevalence & Distribution ####


#Prevalence
PrevSpx <- sum(BAT_inf$Spx_tot>=1)/nrow(BAT_inf)


#Distribution --> Count so poisson family
svg("./Distrib_SPX_log.svg")
hist(
  log(BAT_inf$Spx_tot+1),
  breaks=100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(log(BAT_inf$Spx_tot+1), na.rm = TRUE), max(log(BAT_inf$Spx_tot+1), na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1)
dev.off()


### 4.2.2 Family selection ####

model_Snb1 <- glmmTMB(Spx_tot ~ HZ + Age + Sex + Age:Sex + (1|Site), data=BAT_inf, 
                      family=nbinom1)
model_Snb2 <- glmmTMB(Spx_tot ~ HZ + Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                      family=nbinom2)
model_Sp <- glmmTMB(Spx_tot ~ HZ +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                    family=poisson)
model_Sg <- glmmTMB(log(Spx_tot+1) ~ HZ +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                    family=gaussian)
model_Sgsqrt <- glmmTMB(sqrt(Spx_tot) ~ HZ +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                    family=gaussian)

AIC(model_Snb1, model_Snb2, model_Sp, model_Sg, model_Sgsqrt)
#FUNI, FAS, HZ: Gaussian + log = best model


## Check normality and homoscedasticity
hist(residuals(model_Sg)) #Ok for normality
plot(predict(model_Sg, type = "response"),resid(model_Sg, type = "pearson")) #OK for homoscedasticity
check_model(model_Sg)
#FUNI, FAS, HZ: ok


### 4.2.3 Model selection ####

##Removing site
model_Sg_site <- glmmTMB(log(Spx_tot+1) ~ HZ + Age + Sex + Age:Sex, data=BAT_inf,
                         family=gaussian)
anova(model_Sg, model_Sg_site) 
#FUNI, FAS, HZ: Site significant


##Removing interaction Age * Sex
model_Sg_inter <- glmmTMB(log(Spx_tot+1) ~ HZ + Age + Sex + (1|Site), data=BAT_inf,
                          family=gaussian)
anova(model_Sg, model_Sg_inter) 
#FUNI, FAS, HZ: Interaction significant


##Removing FUNI, FAS, HZ
model_Sg_FUNI <- glmmTMB(log(Spx_tot+1) ~ Sex*Age + (1|Site), data=BAT_inf,
                         family=gaussian)
anova(model_Sg, model_Sg_FUNI)
#FUNI, FAS, HZ not significant



### 4.2.4 Best selection ####
#### 4.2.4.1 FUNI & HZ ####
model_Sg_FUNI <- glmmTMB(log(Spx_tot+1) ~ Sex * Age + (1|Site), data=BAT_inf,
                         family=gaussian)
summary(model_Sg_FUNI)
#FUNI: SexM + SexM:AgeJUV

pairs(emmeans(model_Sg_FUNI,"Sex"), adjust="tukey")#F more parasitized than M
pairs(emmeans(model_Sg_FUNI,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD + no difference between AD & SAD
pairs(emmeans(model_Sg_FUNI,pairwise~ Age*Sex), type="est", adjust="tukey", tran="log")
#ADF > ADM + SADM (but no diff w/in F)
#SADF > ADM + SADM
#JUVF > ADM + SADM (but no diff between JUVF & M)
#JUVM > ADM + SADM


plot_model(model_Sg_FUNI, type="est")
plot_model(model_Sg_FUNI, type="emm", terms=c("Age", "Sex"), show.data = F)+
  theme_classic()


#### 4.2.4.2 FAS ####
model_Sg <- glmmTMB(log(Spx_tot+1) ~ FAS + Sex * Age + (1|Site), data=BAT_inf,
                    family=gaussian)
summary(model_Sg)
#FAS: FAS + SexM + SexM:AgeJUV

pairs(emmeans(model_Sg,"FAS"), adjust="tukey")#FAS is continuous
pairs(emmeans(model_Sg,"Sex"), adjust="tukey")#F more parasitized than M
pairs(emmeans(model_Sg,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD + no difference between AD & SAD
pairs(emmeans(model_Sg,pairwise~ Age*Sex), type="est", adjust="tukey", tran="log")
#ADF > ADM + SADM (but no diff w/in F)
#SADF > ADM and SADM 0.054
#JUVF > ADM + SADM (but no diff between JUVF & M)
#JUVM > ADM + SADM

plot_model(model_Sg, type="est")
plot_model(model_Sg, type="emm", terms=c("Age", "Sex"), show.data = F)

##HOW TO DO WITH CONTINUOUS DATA


## 4.3 Polychromophilus abundance ####
### 4.3.1 Prevalence & Distribution ####

#Prevalence
PrevPoly <- sum(BAT_inf$INF_POLY>0, na.rm=TRUE)/nrow(BAT_inf)


#Ditribution POLY --> Log(ratio): ???
svg("./Distrib_POLY_sqrt.svg")
hist(
  BAT_inf$INF_POLY,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(BAT_inf$INF_POLY, na.rm = TRUE), max(BAT_inf$INF_POLY, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1)
dev.off()

### 4.3.2 Family selection ####

BAT_inf_clean <- BAT_inf %>% filter(!is.na(BAT_inf$INF_POLY))
hist(BAT_inf_clean$INF_POLY, breaks=50)

PrevPoly2 <- sum(BAT_inf_clean$INF_POLY>0)/nrow(BAT_inf_clean)#0.518

model_PziB <- glmmTMB(INF_POLY ~ HZ + Age * Sex + (1|Site), data=BAT_inf_clean,
                     family=beta_family, ziformula = ~ FUNI + Age * Sex)
model_Pt <- glmmTMB(INF_POLY ~ HZ + Age * Sex + (1|Site), data=BAT_inf_clean,
                      family=tweedie)


plot_model(model_PziB, type = "pred", terms = c("Age","Sex"), show.data = F)
plot_model(model_Pt, type = "pred", terms = c("Age","Sex"), show.data = F)

check_model(model_Pt)
check_distribution(model_Pt)
mod_sim=simulateResiduals(model_Pt)
testZeroInflation(mod_sim)
plot(mod_sim)
summary(model_Pt)

AIC(model_PziB, model_Pt)
#FUNI, FAS, HZ: Tweedie = best model


## Check normality and homoscedasticity
hist(residuals(model_Pt)) #NO
plot(predict(model_PziB, type = "response"),resid(model_PziB, type = "pearson")) #NO but better Beta than Tweedie


### 4.3.3 Model selection ####

##Removing site
model_Pt_site <- glmmTMB(INF_POLY ~ HZ + Age * Sex, data=BAT_inf_clean,
                         family=tweedie)
anova(model_Pt, model_Pt_site) 
#FUNI, FAS, HZ: Site significant


##Removing interaction Age * Sex
model_Pt_inter <- glmmTMB(INF_POLY ~ HZ + Age + Sex + (1|Site), data=BAT_inf_clean,
                            family=tweedie)
anova(model_Pt, model_Pt_inter) 
#FUNI, FAS, HZ: Interaction not significant


##Removing Age
model_Pt_age <- glmmTMB(INF_POLY ~ HZ + Sex + (1|Site), data=BAT_inf_clean,
                          family=tweedie)
anova(model_Pt_inter, model_Pt_age) 
#FUNI, FAS, HZ: Age significant


##Removing Sex
model_Pt_sex <- glmmTMB(INF_POLY ~ HZ + Age + (1|Site), data=BAT_inf_clean,
                          family=tweedie)
anova(model_Pt_inter, model_Pt_sex) 
#FUNI, FAS, HZ: Sex significant


##Removing FUNI
model_Pt_FUNI <- glmmTMB(INF_POLY ~ Age + Sex + (1|Site), data=BAT_inf_clean,
                           family=tweedie)
anova(model_Pt_inter, model_Pt_FUNI) 
#FUNI, FAS, HZ not significant


### 4.3.4 Best selection ####

model_Pt_FUNI <- glmmTMB(INF_POLY ~ Sex + Age + (1|Site), data=BAT_inf_clean,
                         family=tweedie)
summary(model_Pt_FUNI)#SexM + AgeSAD + AgeJUV

pairs(emmeans(model_Pt_FUNI,"Sex"), adjust="tukey") #F more parasitized than M
pairs(emmeans(model_Pt_FUNI,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD + SAD more than AD


## 4.4 Polychromophilus presence ####
### 4.4.1 Model selection ####

model_Ppb <- glmmTMB(PREV_POLY ~ FAS + Age * Sex + (1|Site), data=BAT_inf_clean,
                    family=binomial)

#Removing site
model_Ppb_site <- glmmTMB(PREV_POLY ~ FAS + Age * Sex, data=BAT_inf_clean,
                     family=binomial)
anova(model_Ppb, model_Ppb_site) #not significant

#Removing interaction
model_Ppb_inter <- glmmTMB(PREV_POLY ~ FAS + Age + Sex, data=BAT_inf_clean,
                          family=binomial)
anova(model_Ppb_site, model_Ppb_inter) #not significant p=0.722

#Removing age
model_Ppb_age <- glmmTMB(PREV_POLY ~ FAS + Sex, data=BAT_inf_clean,
                           family=binomial)
anova(model_Ppb_inter, model_Ppb_age) #significant p=9.958e-10 ***; chisq = 41.455; chidf = 2

#Removing sex
model_Ppb_sex <- glmmTMB(PREV_POLY ~ FAS + Age, data=BAT_inf_clean,
                         family=binomial)
anova(model_Ppb_inter, model_Ppb_sex) #significant p=0.00244 **; chisq = 9.1849; chidf = 1

#Removing gendiv
model_Ppb_gendiv <- glmmTMB(PREV_POLY ~ Sex + Age, data=BAT_inf_clean,
                         family=binomial)
anova(model_Ppb_inter, model_Ppb_gendiv) 
# HZ: not significant p=0.4234; chisq = 0.6408; Chidf = 1
# FUNI: not significant p=0.3495; chisq = 0.8751; Chidf = 1
# FAS: not significant p=0.966; chisq = 0.0018; Chidf = 1


### 4.4.2 Best selection ####

model_Ppb_gendiv <- glmmTMB(PREV_POLY ~ Sex + Age, data=BAT_inf_clean,
                            family=binomial)
summary(model_Ppb_gendiv)#SexM + AgeSAD + AgeJUV

pairs(emmeans(model_Ppb_gendiv,"Sex"), adjust="tukey") #F more parasitized than M
pairs(emmeans(model_Ppb_gendiv,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD + SAD more than AD


## 4.5 Polychromophilus intensity ####
### 4.5.1 New dataset ####

BAT_inf_clean2 <- BAT_inf_clean[BAT_inf_clean$PREV_POLY == 1,]

### 4.5.2 Model selection ####

model_Pig <- glmmTMB(log(INF_POLY) ~ FAS + Age * Sex + (1|Site), data=BAT_inf_clean2,
                     family=gaussian)

#Removing site
model_Pig_site <- glmmTMB(log(INF_POLY) ~ FAS + Age * Sex, data=BAT_inf_clean2,
                          family=gaussian)
anova(model_Pig, model_Pig_site) #significant

#Removing interaction
model_Pig_inter <- glmmTMB(log(INF_POLY) ~ FAS + Age + Sex + (1|Site), data=BAT_inf_clean2,
                           family=gaussian)
anova(model_Pig_site, model_Pig_inter) #not significant p=1

#Removing age
model_Pig_age <- glmmTMB(log(INF_POLY) ~ FAS + Sex + (1|Site), data=BAT_inf_clean2,
                         family=gaussian)
anova(model_Pig_inter, model_Pig_age) #significant p=0.0002701 ***; chisq = 16.434; chidf = 2

#Removing sex
model_Pig_sex <- glmmTMB(log(INF_POLY) ~ FAS + Age + (1|Site), data=BAT_inf_clean2,
                         family=gaussian)
anova(model_Pig_inter, model_Pig_sex) #significant p=0.04091 *; chisq = 4.1797; chidf = 1

#Removing gendiv
model_Pig_gendiv <- glmmTMB(log(INF_POLY) ~ Sex + Age + (1|Site), data=BAT_inf_clean2,
                            family=gaussian)
anova(model_Pig_inter, model_Pig_gendiv) 
# HZ: not significant p=0.576; chisq = 0.3127; Chidf = 1
# FUNI: not significant p=0.9024; chisq = 0.015; Chidf = 1
# FAS: not significant p=0.6123; chisq = 0.2568; Chidf = 1

### 4.5.3 Best selection ####

model_Pig_gendiv <- glmmTMB(log(INF_POLY) ~ Sex + Age + (1|Site), data=BAT_inf_clean2,
                            family=gaussian)
summary(model_Pig_gendiv)#SexM + AgeJUV

pairs(emmeans(model_Pig_gendiv,"Sex"), adjust="tukey") #F more parasitized than M
pairs(emmeans(model_Pig_gendiv,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD



################################################################################
# 5. Visualization of the data ####
################################################################################
## 5.1 Figures ####
### 5.1.1 Genetic diversity and parasitism ####
#### 5.1.1.1 FUNI ####

library(ggplot2)
library(dplyr)
library(tidyr)

scale_ratio <- mean(c(BAT_inf$Spx_tot, BAT_inf$Nyct_tot), na.rm = TRUE) / mean(BAT_inf$INF_POLY, na.rm = TRUE)


# Reshape Spx/Nyct into long format
BAT_long <- BAT_inf %>%
  pivot_longer(cols = c(Spx_tot, Nyct_tot),
               names_to = "Parasite",
               values_to = "Abundance")

# Add INF_POLY as a separate category to include in the legend
BAT_poly <- BAT_inf %>%
  mutate(Parasite = "INF_POLY",
         Abundance = INF_POLY * scale_ratio) %>%
  select(FUNI, Parasite, Abundance)

# Combine data
BAT_combined <- bind_rows(
  BAT_long %>% mutate(Abundance = Abundance),
  BAT_poly
)

# Scientific color palette
cbScientific <- c("Spx_tot" = "#332288", "Nyct_tot" = "#E69F00", "INF_POLY" = "#6DCD4A")

svg("./FUNI_PARA.svg")
# Plot
ggplot(BAT_combined, aes(x = FUNI, y = Abundance, color = Parasite)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = cbScientific) +
  scale_y_continuous(
    name = "Ectoparasite abundance (bat flies & wing mites)",
    sec.axis = sec_axis(~./scale_ratio,
                        name = "P. murinus relative abundance",
                        breaks = pretty(BAT_inf$INF_POLY))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    axis.title.y.right = element_text(color = "black")  # ← changed from green to black
  )
dev.off()


#### 5.1.1.2 HeteroZ ####

library(ggplot2)
library(dplyr)
library(tidyr)

# Reshape Spx/Nyct into long format
BAT_long <- BAT_inf %>%
  pivot_longer(cols = c(Spx_tot, Nyct_tot),
               names_to = "Parasite",
               values_to = "Abundance")

# Add INF_POLY as a separate category to include in the legend
BAT_poly <- BAT_inf %>%
  mutate(Parasite = "INF_POLY",
         Abundance = INF_POLY * scale_ratio) %>%
  select(HZ, Parasite, Abundance)

# Combine data
BAT_combined <- bind_rows(
  BAT_long %>% mutate(Abundance = Abundance),
  BAT_poly
)

# Scientific color palette
cbScientific <- c("Spx_tot" = "#332288", "Nyct_tot" = "#E69F00", "INF_POLY" = "#6DCD4A")

svg("./HZ_PARA.svg")
# Plot
ggplot(BAT_combined, aes(x = HZ, y = Abundance, color = Parasite)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = cbScientific) +
  scale_y_continuous(
    name = "Ectoparasite abundance (bat flies & wing mites)",
    sec.axis = sec_axis(~./scale_ratio,
                        name = "P. murinus relative abundance",
                        breaks = pretty(BAT_inf$INF_POLY))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    axis.title.y.right = element_text(color = "black")  # ← changed from green to black
  )
dev.off()




#### 5.1.1.2 FAS ####

library(ggplot2)
library(dplyr)
library(tidyr)

# Reshape Spx/Nyct into long format
BAT_long <- BAT_inf %>%
  pivot_longer(cols = c(Spx_tot, Nyct_tot),
               names_to = "Parasite",
               values_to = "Abundance")

# Add INF_POLY as a separate category to include in the legend
BAT_poly <- BAT_inf %>%
  mutate(Parasite = "INF_POLY",
         Abundance = INF_POLY * scale_ratio) %>%
  select(FAS, Parasite, Abundance)

# Combine data
BAT_combined <- bind_rows(
  BAT_long %>% mutate(Abundance = Abundance),
  BAT_poly
)

# Scientific color palette
cbScientific <- c("Spx_tot" = "#332288", "Nyct_tot" = "#E69F00", "INF_POLY" = "#6DCD4A")

svg("./FAS_PARA.svg")
# Plot
ggplot(BAT_combined, aes(x = FAS, y = Abundance, color = Parasite)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = cbScientific) +
  scale_y_continuous(
    name = "Ectoparasite abundance (bat flies & wing mites)",
    sec.axis = sec_axis(~./scale_ratio,
                        name = "P. murinus relative abundance",
                        breaks = pretty(BAT_inf$INF_POLY))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    axis.title.y.right = element_text(color = "black")  # ← changed from green to black
  )
dev.off()


### 5.1.2 Bat flies ####

svg("./NYCT_Sex.svg")
ggplot(BAT_inf, aes(x = Sex, y = Nyct_tot)) +
  geom_boxplot(fill = "grey", color = "black", outlier.shape = 16, outlier.size = 2) +
  labs(
    x = "Sex",
    y = "Nyct_tot (Parasite abundance)",
    title = "Nyct_tot Abundance by Sex"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  ) +
  geom_signif(comparisons = list(c("F", "M")),
              annotations = "*",
              y_position = 23)
dev.off()


### 5.1.3 Wing mites ####

svg("./SPX_SexAge.svg")
ggplot(BAT_inf, aes(x = Sex, y = log(Spx_tot+1), fill=Age)) +
  geom_boxplot()+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  ) +
  scale_fill_brewer(palette="Greys")
dev.off()


### 5.1.4 Polychromophilus ####

svg("./POLY_SexAge.svg")
ggplot(BAT_inf, aes(x = Sex, y = INF_POLY, fill=Age)) +
  geom_boxplot()+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  ) +
  scale_fill_brewer(palette="Greys") +
  geom_signif(comparisons = list(c("F", "M")),
              annotations = "**",
              y_position = 0.6)
dev.off()


library(ggsignif)
svg("./POLY_Sex.svg")
ggplot(BAT_inf, aes(x = Sex, y = INF_POLY)) +
  geom_boxplot(fill="grey")+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  ) +
  geom_signif(comparisons = list(c("F", "M")),
              annotations = "**",
              y_position = 0.6)
dev.off()


svg("./POLY_Age.svg")
ggplot(BAT_inf, aes(x = Age, y = INF_POLY, fill=Age)) +
  geom_boxplot()+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  ) +
  scale_fill_brewer(palette="Greys")
dev.off()


################################################################################
# 6. MHC region ####
################################################################################
## 6.1 How many SNPs in the MHC region? ####

vcf_tot <- read.vcfR("../../vcf_after_filtering.recode.vcf")
pos<-vcf_tot@fix
pos<-as.data.frame(pos)

##Define loci
#class II MHC transactivator
cond1 <- pos$CHROM == "NC_081843.1" & pos$POS > 21900000 & pos$POS < 22100000
#class II antigens
cond2 <- pos$CHROM == "NC_081845.1" & pos$POS > 77000000 & pos$POS < 77400000
cond3 <- pos$CHROM == "NC_081845.1" & pos$POS > 94800000 & pos$POS < 97600000
#class I antigens
cond4 <- pos$CHROM == "NC_081854.1" & pos$POS > 104000 & pos$POS < 107000
cond5 <- pos$CHROM == "NC_081854.1" & pos$POS > 1570000 & pos$POS < 2140000
cond6 <- pos$CHROM == "NC_081854.1" & pos$POS > 9490000 & pos$POS < 9500000

# Combine all conditions
keepSNPsMHC <- cond1 | cond2 | cond3 | cond4 | cond5 | cond6
MHC_snp <- pos$ID[keepSNPsMHC] #244 SNPs in the MHC


## 6.2 Relative heterozygosity of this region ####
### 6.2.1 Calculate HetZ_MHC ####

gen_mat <- snpgdsGetGeno('vcf_after_filtering_BAT.gds')

IndHet2 = rowSums(gen_mat==1, na.rm = TRUE) / rowSums(! is.na(gen_mat))
SNPHet2 = colSums(gen_mat==1, na.rm = TRUE) / colSums(! is.na(gen_mat))
#plot(SNPHet)

gen_matMHC = gen_mat[,keepSNPsMHC]
IndHetMHC = rowSums(gen_matMHC==1, na.rm = TRUE) / rowSums(! is.na(gen_matMHC))


### 6.2.2 Compare with overall dataset ####

layout(matrix(c(1,2)))
hist(IndHet2, xlim = c(0,0.5), breaks = seq(0,0.5, 0.01))
hist(IndHetMHC, xlim = c(0,0.5), breaks = seq(0,0.5, 0.01))
dev.off()


# bootstrap 244 random

keepSNPsRand244 = sample(1:ncol(gen_mat), 244, replace = FALSE)
gen_matRand244 = gen_mat[,keepSNPsRand244]
IndHetRand244 = rowSums(gen_matRand244==1, na.rm = TRUE) / rowSums(! is.na(gen_matRand244))

#Allow to add a color to plot
add.alpha <- function(col, alpha=1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  rgb_vals <- sapply(col, col2rgb) / 255
  result <- apply(rgb_vals, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
  return(result)
}

svg("./Comparison_HetZ_MHC-overall.svg")
layout(matrix(c(1,2)))
plot(0,0, type = "n", xlim = c(0,0.5), ylim = c(0,100),
     ylab="Count", xlab="Relative heterozygosity - Overall", axes=F)
axis(side=2, at=c(0, 20, 40, 60, 80, 100))
axis (side=1, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
for (i in 1:100){
  keepSNPsRand244 = sample(1:ncol(gen_mat), 29, replace = FALSE)
  gen_matRand244 = gen_mat[,keepSNPsRand244]
  IndHetRand244 = rowSums(gen_matRand244==1, na.rm = TRUE) / rowSums(! is.na(gen_matRand244))
  
  hist(IndHetRand244, col = add.alpha("black", 0.01), border = F, breaks = seq(0,0.6, 0.025), add = TRUE)
}
# hist(IndHet, xlim = c(0,0.5), breaks = seq(0,0.5, 0.01))
hist(IndHetMHC, xlim = c(0,0.5), breaks = seq(0,0.5, 0.025), border=F,
     main="", xlab="Relative heterozygosity - MHC Region", ylab="Count")
dev.off()

### 6.2.3 Add to dataset ####

BAT_inf$HZ_MHC <- IndHetMHC


## 6.3 Models ####
### 6.3.1 Explore ####

plot(BAT_inf$HZ_MHC~log(BAT_inf$Nyct_tot+1))
plot(BAT_inf$HZ_MHC~log(BAT_inf$Spx_tot+1))
plot(BAT_inf$HZ_MHC~log(BAT_inf$INF_POLY+1))


### 6.3.2 Bat flies ####
#### 6.3.2.1 Family selection ####

model_Nnb1_MHC <- glmmTMB(Nyct_tot ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf, 
                      family=nbinom1)
model_Nnb2_MHC <- glmmTMB(Nyct_tot ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf,
                      family=nbinom2)
model_Np_MHC <- glmmTMB(Nyct_tot ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf,
                    family=poisson)
model_Ng_MHC <- glmmTMB(log(Nyct_tot+1) ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf,
                    family=gaussian)
model_Ngsqrt_MHC <- glmmTMB(sqrt(Nyct_tot) ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf,
                        family=gaussian)
AIC(model_Nnb1_MHC, model_Nnb2_MHC, model_Np_MHC, model_Ng_MHC, model_Ngsqrt_MHC)
#Gaussian + log = best model


## Check normality and homoscedasticity
hist(residuals(model_Ng_MHC)) #Ok for normality
plot(predict(model_Ng_MHC, type = "response"),resid(model_Ng_MHC, type = "pearson")) #OK for homoscedasticity
check_model(model_Ng_MHC)#ok


#### 6.3.2.2 Model selection ####

##Removing site
model_Ng_MHC_site <- glmmTMB(log(Nyct_tot+1) ~ HZ_MHC + Age * Sex, data=BAT_inf,
                         family=gaussian)
anova(model_Ng_MHC, model_Ng_MHC_site)
#Site significant


##Removing interaction Age * Sex
model_Ng_MHC_inter <- glmmTMB(log(Nyct_tot+1) ~ HZ_MHC + Age + Sex + (1|Site), data=BAT_inf,
                          family=gaussian)
anova(model_Ng_MHC, model_Ng_MHC_inter) 
#Interaction not significant


##Removing Sex
model_Ng_MHC_sex <- glmmTMB(log(Nyct_tot+1) ~ HZ_MHC + Age + (1|Site), data=BAT_inf,
                        family=gaussian)
anova(model_Ng_MHC_inter, model_Ng_MHC_sex) 
#Sex significant (Chi2 = 5.3104, Chi_df = 1, p = 0.0212*)


##Removing Age
model_Ng_MHC_age <- glmmTMB(log(Nyct_tot+1) ~ HZ_MHC + Sex + (1|Site), data=BAT_inf,
                        family=gaussian)
anova(model_Ng_MHC_inter, model_Ng_MHC_age) 
#Age not significant (Chi2 = 4.3624, Chi_df = 2, p = 0.1129)


##Removing HZ_MHC
model_Ng_MHC_HZ <- glmmTMB(log(Nyct_tot+1) ~ Sex + (1|Site), data=BAT_inf,
                         family=gaussian)
anova(model_Ng_MHC_age, model_Ng_MHC_HZ) 
#HZ not significant (Chi2 = 0.4283, Chi_df = 1, p = 0.5128)



#### 6.3.2.3 Best selection ####

model_Ng_MHC_HZ <- glmmTMB(log(Nyct_tot+1) ~ Sex + (1|Site), data=BAT_inf,
                         family=gaussian)
summary(model_Ng_MHC_HZ)#SexM
pairs(emmeans(model_Ng_MHC_HZ,"Sex"), adjust="tukey") #F more parasitized than M (est = 0.157, df = 355, t = 2.328, p = 0.0205)

plot_model(model_Ng_MHC_HZ, type = "emm", terms = "Sex", show.data = F)


### 6.3.3 Wing mites ####
#### 6.3.3.1 Family selection ####

model_Snb1_MHC <- glmmTMB(Spx_tot ~ HZ_MHC + Age + Sex + Age:Sex + (1|Site), data=BAT_inf, 
                      family=nbinom1)
model_Snb2_MHC <- glmmTMB(Spx_tot ~ HZ_MHC + Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                      family=nbinom2)
model_Sp_MHC <- glmmTMB(Spx_tot ~ HZ +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                    family=poisson)
model_Sg_MHC <- glmmTMB(log(Spx_tot+1) ~ HZ_MHC +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                    family=gaussian)
model_Sgsqrt_MHC <- glmmTMB(sqrt(Spx_tot) ~ HZ_MHC +  Age + Sex + Age:Sex + (1|Site), data=BAT_inf,
                        family=gaussian)

AIC(model_Snb1_MHC, model_Snb2_MHC, model_Sp_MHC, model_Sg_MHC, model_Sgsqrt_MHC)
#Gaussian + log = best model


## Check normality and homoscedasticity
hist(residuals(model_Sg_MHC)) #Ok for normality
plot(predict(model_Sg_MHC, type = "response"),resid(model_Sg_MHC, type = "pearson")) #OK for homoscedasticity
check_model(model_Sg_MHC)#ok


#### 6.3.3.2 Model selection ####

##Removing site
model_Sg_MHC_site <- glmmTMB(log(Spx_tot+1) ~ HZ_MHC + Age + Sex + Age:Sex, data=BAT_inf,
                         family=gaussian)
anova(model_Sg_MHC, model_Sg_MHC_site) 
#Site significant


##Removing interaction Age * Sex
model_Sg_MHC_inter <- glmmTMB(log(Spx_tot+1) ~ HZ_MHC + Age + Sex + (1|Site), data=BAT_inf,
                          family=gaussian)
anova(model_Sg_MHC, model_Sg_MHC_inter) 
#Interaction significant (Chi2 = 17.622, Chi_df = 2, p = 0.0001491***)


##Removing HZ_MHC
model_Sg_MHC_HZ <- glmmTMB(log(Spx_tot+1) ~ Sex*Age + (1|Site), data=BAT_inf,
                         family=gaussian)
anova(model_Sg_MHC, model_Sg_MHC_HZ)
#HZ not significant (Chi2 = 0.0032, Chi_df = 1, p = 0.9552)



#### 6.3.3.3 Best selection ####

model_Sg_MHC_HZ <- glmmTMB(log(Spx_tot+1) ~ Sex*Age + (1|Site), data=BAT_inf,
                           family=gaussian)
summary(model_Sg_MHC_HZ)#SexM + SexM:AgeJUV
pairs(emmeans(model_Sg_MHC_HZ,"Sex"), adjust="tukey") #F more parasitized than M (est = 0.378, df = 351, t = 4.386, p = <0.0001)
pairs(emmeans(model_Sg_MHC_HZ,"Age"), adjust="tukey") #AD less parasitized than JUV (est = -0.5257, df = 351, t = -5.380, p = <0.0001)
                                                      #SAD less parasitized than JUV (est = -0.4763, df = 351, t = -4.162, p = 0.0001)
                                                      #No diff AD-SAD (est = -0.0494, df = 351, t = -0.525, p = 0.8589)
pairs(emmeans(model_Sg_MHC_HZ,~Age*Sex), adjust="tukey")

plot_model(model_Sg_MHC_HZ, type = "emm", terms = c("Sex", "Age"), show.data = F)


### 6.3.4 Polychromophilus ####
#### 6.3.4.1 Family selection ####

BAT_inf_clean <- BAT_inf %>% filter(!is.na(BAT_inf$INF_POLY))
hist(BAT_inf_clean$INF_POLY, breaks=50)

PrevPoly2 <- sum(BAT_inf_clean$INF_POLY>0)/nrow(BAT_inf_clean)#0.518

model_PziB_MHC <- glmmTMB(INF_POLY ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf_clean,
                      family=beta_family, ziformula = ~ FUNI + Age * Sex)
model_Pt_MHC <- glmmTMB(INF_POLY ~ HZ_MHC + Age * Sex + (1|Site), data=BAT_inf_clean,
                    family=tweedie)


plot_model(model_PziB_MHC, type = "pred", terms = c("Age","Sex"), show.data = F)
plot_model(model_Pt_MHC, type = "pred", terms = c("Age","Sex"), show.data = F)

check_model(model_Pt_MHC)
check_distribution(model_Pt_MHC)
mod_sim_MHC=simulateResiduals(model_Pt_MHC)
testZeroInflation(mod_sim_MHC)
plot(mod_sim_MHC)
summary(model_Pt_MHC)

AIC(model_PziB_MHC, model_Pt_MHC)
#Tweedie = best model


## Check normality and homoscedasticity
hist(residuals(model_Pt_MHC)) #NO
plot(predict(model_Pt_MHC, type = "response"),resid(model_Pt_MHC, type = "pearson")) #NO but better Beta than Tweedie


#### 6.3.4.2 Model selection ####

##Removing site
model_Pt_MHC_site <- glmmTMB(INF_POLY ~ HZ_MHC + Age * Sex, data=BAT_inf_clean,
                         family=tweedie)
anova(model_Pt_MHC, model_Pt_MHC_site) 
#Site significant


##Removing interaction Age * Sex
model_Pt_MHC_inter <- glmmTMB(INF_POLY ~ HZ_MHC + Age + Sex + (1|Site), data=BAT_inf_clean,
                          family=tweedie)
anova(model_Pt_MHC, model_Pt_MHC_inter) 
#Interaction not significant (Chi2 = 0.4453, Chi_df = 2, p = 0.8004)


##Removing Age
model_Pt_MHC_age <- glmmTMB(INF_POLY ~ HZ_MHC + Sex + (1|Site), data=BAT_inf_clean,
                        family=tweedie)
anova(model_Pt_MHC_inter, model_Pt_MHC_age) 
#Age significant (Chi2 = 36.93, Chi_df = 2, p = 9.567e-09***)


##Removing Sex
model_Pt_MHC_sex <- glmmTMB(INF_POLY ~ HZ_MHC + Age + (1|Site), data=BAT_inf_clean,
                        family=tweedie)
anova(model_Pt_MHC_inter, model_Pt_MHC_sex) 
#Sex significant  (Chi2 = 8.3217, Chi_df = 1, p = 0.003917**)


##Removing HZ_MHC
model_Pt_MHC_HZ <- glmmTMB(INF_POLY ~ Age + Sex + (1|Site), data=BAT_inf_clean,
                         family=tweedie)
anova(model_Pt_MHC_inter, model_Pt_MHC_HZ) 
#HZ not significant  (Chi2 = 0.2634, Chi_df = 1, p = 0.6078)


#### 6.3.4.3 Best selection ####

model_Pt_MHC_HZ <- glmmTMB(INF_POLY ~ Age + Sex + (1|Site), data=BAT_inf_clean,
                           family=tweedie)
summary(model_Pt_MHC_HZ)#SexM + AgeSAD + AgeJUV

pairs(emmeans(model_Pt_MHC_HZ,"Sex"), adjust="tukey") #F more parasitized than M 
pairs(emmeans(model_Pt_MHC_HZ,"Age"), adjust="tukey") #JUV more parasitized than AD & SAD + SAD more than AD


## 6.4 Graphical representation ####

library(ggplot2)
library(dplyr)
library(tidyr)

scale_ratio <- mean(c(BAT_inf$Spx_tot, BAT_inf$Nyct_tot), na.rm = TRUE) / mean(BAT_inf$INF_POLY, na.rm = TRUE)

# Reshape Spx/Nyct into long format
BAT_long2 <- BAT_inf %>%
  pivot_longer(cols = c(Spx_tot, Nyct_tot),
               names_to = "Parasite",
               values_to = "Abundance")

# Add INF_POLY as a separate category to include in the legend
BAT_poly2 <- BAT_inf %>%
  mutate(Parasite = "INF_POLY",
         Abundance = INF_POLY * scale_ratio) %>%
  select(HZ_MHC, Parasite, Abundance)

# Combine data
BAT_combined2 <- bind_rows(
  BAT_long2 %>% mutate(Abundance = Abundance),
  BAT_poly2
)

# Scientific color palette
cbScientific <- c("Spx_tot" = "#332288", "Nyct_tot" = "#E69F00", "INF_POLY" = "#6DCD4A")

svg("./HZ_MHC_PARA.svg")
# Plot
ggplot(BAT_combined2, aes(x = HZ_MHC, y = Abundance, color = Parasite)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = cbScientific) +
  scale_y_continuous(
    name = "Ectoparasite abundance (bat flies & wing mites)",
    sec.axis = sec_axis(~./scale_ratio,
                        name = "P. murinus relative abundance",
                        breaks = pretty(BAT_inf$INF_POLY))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    axis.title.y.right = element_text(color = "black")  # ← changed from green to black
  )
dev.off()



