#####################################################
##### Chapter 4 - Analysis - Infection Dynamics #####
#####################################################

setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/5.Thesis/Project/8.Infection dynamic/Analysis")

load("./InfDyn_lastversion.RData") ## Load last R history


################################################################################
# 1. Download packages and data ####
################################################################################
## 1.1 Download packages ####

library(mgcv)
library(ggplot2)
library(glmmTMB)
library(sjPlot)
library(emmeans)
library(performance)
library(piecewiseSEM)
library(cowplot)
library(dplyr)
library(tidyr)

## 1.2 Download data ####
### 1.2.1 Dataset ####

Dodo_Parasito <- read.csv("./Capture_dodo_ALL.csv", h = T, sep = ";")

### 1.2.2 Create prevalences & quantities ####

Dodo_Parasito$INF_POLY <- Dodo_Parasito$meanQuantity_Poly/(Dodo_Parasito$meanQuantity_Poly+Dodo_Parasito$meanQuantity_Bat)
Dodo_Parasito$PREV_POLY <- ifelse(Dodo_Parasito$INF_POLY > 0,1,0)


### 1.2.3 Clean dataset ####

## Convert Date
Dodo_Parasito$Date <- as.Date(Dodo_Parasito$Date, format = "%d/%m/%Y")


## Clean dataset for models
Dodo_Parasito1 <- Dodo_Parasito[!is.na(Dodo_Parasito$Date)
                                & !is.na(Dodo_Parasito$Site)
                                & !is.na(Dodo_Parasito$Sex)
                                & !is.na(Dodo_Parasito$Repro)
                                & !is.na(Dodo_Parasito$Age)
                                & !is.na(Dodo_Parasito$Spx_tot)
                                & !is.na(Dodo_Parasito$Nyct_tot),]

Dodo_Parasito1$Date_numeric <- as.numeric(Dodo_Parasito1$Date)
Dodo_Parasito1$Year <- format(Dodo_Parasito1$Date, "%Y")
#Dodo_Parasito1$Year <- factor(Dodo_Parasito1$Year)
#Dodo_Parasito1$Site <- factor(Dodo_Parasito1$Site)
Dodo_Parasito1$MonthDay <- format(Dodo_Parasito1$Date, "%m/%d")


################################################################################
# 2. Temporal analyses - only on 2023 - DATE ####
################################################################################
## 2.1 Datasets ####

## 2023 only
Dodo_Parasito2023 <- Dodo_Parasito1[Dodo_Parasito1$Year == "2023",]
table(Dodo_Parasito2023$Site)

## Dorigny and MR only
Dodo_Parasito23_site <- Dodo_Parasito1[Dodo_Parasito1$Year == "2023"&Dodo_Parasito1$Site == c("Dorigny", "MaisonRiviere"),]

### 2.2 Test effect date across all sites ####
#Assumption: there is the same dynamic of parasitism across all sites
#### 2.2.1 SPX ####

overall1 <- gam(Spx_tot ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023, family = poisson())
overall1_gauss <- gam(log(Spx_tot+1) ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023, family = gaussian())
AIC(overall1, overall1_gauss)#normal
summary(overall1_gauss) #Significant effect of date
plot(overall1_gauss,pages=1,residuals=TRUE, xlab="Date", ylab="log(SPX+1)")

#### 2.2.3 NYCT ####

overall11 <- gam(Nyct_tot ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023, family = poisson())
overall11_gauss <- gam(log(Nyct_tot+1) ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023, family = gaussian())
AIC(overall11, overall11_gauss)#normal
summary(overall11_gauss) #No significant effect of date
plot(overall11_gauss,pages=1,residuals=TRUE, xlab="Date", ylab="log(NYCT+1)")
length(Dodo_Parasito2023$INF_POLY[Dodo_Parasito2023$INF_POLY > 0])


#### 2.2.4 POLY ####

#Dodo_Parasito2023b <- Dodo_Parasito2023[!is.na(Dodo_Parasito2023$PREV_POLY),]

overall111 <- gam(PREV_POLY ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023, family = binomial())
summary(overall111) #No significant effect of date
plot(overall111,pages=1,residuals=TRUE, xlab="Date", ylab="Poly_prevalence")

Dodo_Parasito2023_int <- Dodo_Parasito2023[Dodo_Parasito2023$PREV_POLY == 1,]
table(Dodo_Parasito2023_int$Site)
overall111_int <- gam(log(INF_POLY) ~ s(Date_numeric, k = 5), data = Dodo_Parasito2023_int, family = gaussian())
summary(overall111_int) #Significant effect of date
plot(overall111_int,pages=1,residuals=TRUE, xlab="Date", ylab="Poly_intensity")


### 2.3 Plots effect of date over all sites ####

# SPX
p <- ggplot(data = Dodo_Parasito2023, mapping = aes(x = Date_numeric, y = log(Spx_tot+1))) +
  geom_point(size = 0.5, color = "grey") + theme_bw() + ylab(expression("LOG(Wing mite abundance + 1)")) +
  geom_smooth(method="gam", formula = y ~ s(x, k = 5),color = "black",method.args = list(family = "gaussian")) +
  scale_x_continuous(
    name = "Date",
    breaks = Dodo_Parasito2023$Date_numeric[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")],
    labels = format(Dodo_Parasito2023$Date[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")], "%d-%b")  # ou "%Y-%m-%d"
  ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y =element_text(size=8))
p

# NYCT
q <- ggplot(data = Dodo_Parasito2023, mapping = aes(x = Date_numeric, y = log(Nyct_tot+1))) +
  geom_point(size = 0.5, color = "grey") + theme_bw() + ylab(expression("LOG(Bat fly abundance + 1)")) +
  geom_smooth(method="gam", formula = y ~ s(x, k = 5),color = "black",method.args = list(family = "gaussian")) +
  scale_x_continuous(
    name = "Date",
    breaks = Dodo_Parasito2023$Date_numeric[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")],
    labels = format(Dodo_Parasito2023$Date[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")], "%d-%b")  # ou "%Y-%m-%d"
  ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y =element_text(size=8))
q

# POLY intensity
r <- ggplot(data = Dodo_Parasito2023_int, mapping = aes(x = Date_numeric, y = log(INF_POLY))) +
  geom_point(size = 0.5, color = "grey") + theme_bw() + ylab(expression("LOG("*italic(P.murinus)*" intensity)")) +
  geom_smooth(method="gam", formula = y ~ s(x, k = 5),color = "black",method.args = list(family = "gaussian")) +
  scale_x_continuous(
    name = "Date",
    breaks = Dodo_Parasito2023$Date_numeric[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")],
    labels = format(Dodo_Parasito2023$Date[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")], "%d-%b")  # ou "%Y-%m-%d"
  )+ theme(axis.title.y =element_text(size=8))
r

# POLY prev
s <- ggplot(data = Dodo_Parasito2023, mapping = aes(x = Date_numeric, y = PREV_POLY)) +
  geom_point(size = 0.5, color = "grey") + theme_bw() +
  geom_smooth(method="gam", formula = y ~ s(x, k = 5),color = "black",method.args = list(family = "binomial")) +
  scale_y_continuous(
    name = expression(italic(P.murinus)*" presence"),
    breaks =c(0,1),
    labels = c("Absent", "Present")
  ) +
  scale_x_continuous(
    "Date",
    breaks = Dodo_Parasito2023$Date_numeric[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")],
    labels = format(Dodo_Parasito2023$Date[Dodo_Parasito2023$Date %in% c("2023-06-07", "2023-07-06", "2023-08-07", "2023-09-06", "2023-10-04")], "%d-%b")  # ou "%Y-%m-%d"
  )  + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y =element_text(size=8))
s



figure <- plot_grid(p, q, s, r, ncol = 1, align = "v", rel_heights = c(1, 1, 1, 1),labels = c("A", "B", "C", "D"), label_x = 0.92, label_y = 0.95)
figure


ggsave("Fig_temp_dyn.pdf",plot = figure,dpi = "retina",width = 4, height = 8)


table(is.na(Dodo_Parasito2023$PREV_POLY))


### 2.4 Build an date effect per site ####
#Est-ce que effectivement y a la meme dynamique across all sites?

Dodo_Parasito23_site$Site1 <- factor(Dodo_Parasito23_site$Site)
#On ne prend que D et MR car il nous fallait assez de sampling points

#### 2.4.1 SPX ####

overall2 <- gam(Spx_tot ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site, family = poisson())
overall2_gauss <- gam(log(Spx_tot+1) ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site, family = gaussian())
AIC(overall2, overall2_gauss)#gaussian
summary(overall2_gauss) #Significant effect of date and sites
plot(overall2_gauss,pages=1,residuals=TRUE, xlab="Date", ylab="log(SPX+1)")


#### 2.4.2 NYCT ####

overall22 <- gam(Nyct_tot ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site, family = poisson())
overall22_gauss <- gam(log(Nyct_tot+1) ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site, family = gaussian())
AIC(overall22, overall22_gauss)#gaussian
summary(overall22_gauss) #Significant effect of date and sites
plot(overall22_gauss,pages=1,residuals=TRUE, xlab="Date", ylab="log(NYCT+1)")

#### 2.4.3 Poly ####
##### 2.4.3.1 Prevalence ####

overall222 <- gam(PREV_POLY ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site, family = binomial())
summary(overall222) #No significant effect of date
plot(overall222,pages=1,residuals=TRUE, xlab="Date", ylab="Poly_prevalence")


##### 2.4.3.2 Intensity ####

Dodo_Parasito23_site_int <- Dodo_Parasito23_site[Dodo_Parasito23_site$PREV_POLY == 1,]
overall222_int <- gam(log(INF_POLY) ~ Site1 + s(Date_numeric, by = Site1, k = 5), data = Dodo_Parasito23_site_int, family = gaussian())
summary(overall222_int) #Significant effect of date
plot(overall222_int,pages=1,residuals=TRUE, xlab="Date", ylab="Poly_int")


################################################################################
# 3. Temporal analyses - only on 2023 - REPRO TIME ####
################################################################################
## 3.1 New dataset ####

#New time variable based on biology and our data 
Dodo_Parasito1$Repro_time <- ifelse(Dodo_Parasito1$MonthDay < "06/27", "P",
                                    ifelse(Dodo_Parasito1$MonthDay > "07/28", "PL", "L"))

#Remove juveniles to test Repro
Dodo_Parasito11 <- Dodo_Parasito1[Dodo_Parasito1$Age != "JUV",]
table(Dodo_Parasito11$Site,Dodo_Parasito11$Year)

## 3.2 SPX ####
### 3.2.1 Family selection ####

model <- glmmTMB(Spx_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=nbinom1)
model1 <- glmmTMB(Spx_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=nbinom2)
model11 <- glmmTMB(Spx_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=poisson)
model111 <- glmmTMB(log(Spx_tot + 1) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                   family=gaussian)
AIC(model,model1,model11,model111)
#Best model: gaussian 

### 3.2.3 Check model assumption ####

hist(residuals(model111)) # normalité ok
plot(predict(model111, type = "response"),resid(model111, type = "pearson")) # homogénéité des variances ok


### 3.2.4 Model selection ####

model2 <- glmmTMB(log(Spx_tot + 1) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model111,model2) # Site significatif

model2 <- glmmTMB(log(Spx_tot + 1) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model111,model2) # Année significatif

model2 <- glmmTMB(log(Spx_tot + 1) ~Age + Sex + Sex:Repro + Age:Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model111,model2) # Age:Repro non-significatif

model3 <- glmmTMB(log(Spx_tot + 1) ~Age + Sex + Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model2,model3) # Age:Sex non-significatif

model4 <- glmmTMB(log(Spx_tot + 1) ~ Age + Sex + Sex:Repro + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model3,model4) # Repro_time significatif

model4 <- glmmTMB(log(Spx_tot + 1) ~ Age + Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model3,model4) # Repro significatif -> Donc Sex significatif 

model4 <- glmmTMB(log(Spx_tot + 1) ~ Sex + Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model3,model4) # Age non-significatif


### 3.2.5 Best selection ####

pairs(emmeans(model4,"Sex"), adjust="tukey") # Femelles plus parasité que males Mais attention parce que Repo est aussi signif !
pairs(emmeans(model4, ~ Repro | Sex), adjust="tukey") # pas de différence chez les males - Les femelles lactantes sont plus parasité que les pregnant ou les autres. Les post-lactantes ont un niveau intermediaire 
pairs(emmeans(model4, ~ Repro), adjust="tukey") # pas de différence chez les males - Les femelles lactantes sont plus parasité que les pregnant ou les autres. Les post-lactantes ont un niveau intermediaire 
pairs(emmeans(model4,"Repro_time"), adjust="tukey") # courbe en cloche avec plus de parasite en periode de lactation 

## 3.3 NYCT ####
### 3.3.1 Family selection ####

model_Nnb1 <- glmmTMB(Nyct_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=nbinom1)
model_Nnb2 <- glmmTMB(Spx_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=nbinom2)
model_Np <- glmmTMB(Nyct_tot ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                   family=poisson)
model_Ng <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                    family=gaussian)
AIC(model_Nnb1,model_Nnb2,model_Np,model_Ng)
#Best model: gaussian 


### 3.3.2 Check model assumptions ####

hist(residuals(model_Ng)) # normalité ok
plot(predict(model_Ng, type = "response"),resid(model_Ng, type = "pearson")) # homogénéité ok
check_model(model_Ng)


### 3.3.3 Model selection ####

model_Ng_site <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng,model_Ng_site) #Site significant

model_Ng_year <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng,model_Ng_year) #Year significant

model_Ng_AR <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Sex:Repro + Age:Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng,model_Ng_AR) #Age:Repro not significant

model_Ng_AS <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng_AR,model_Ng_AS) #Age:Sex not significant

model_Ng_Repro <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng_AS,model_Ng_Repro) #Repro not significant

model_Ng_Rtime <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng_Repro,model_Ng_Rtime) #Repro_time significant

model_Ng_Sex <- glmmTMB(log(Nyct_tot + 1) ~ Age + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng_Repro,model_Ng_Sex) #Sex significant

model_Ng_Age <- glmmTMB(log(Nyct_tot + 1) ~ Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=gaussian)
anova(model_Ng_Repro,model_Ng_Age) #Age not significant


### 3.3.4 Best model ####


pairs(emmeans(model_Ng_Age,"Sex"), adjust="tukey") #Females more parasitized than males
pairs(emmeans(model_Ng_Age,"Repro_time"), adjust="tukey") #PL>L>P

svg("./Nyct_Sex.svg")
plot_model(model_Ng_Age, type="emm", terms="Sex", show.data = F)+
  theme_classic()+
  xlab("Sex")+
  ylab("log(NYCT+1)")
dev.off()

svg("./Nyct_ReproTime.svg")
plot_model(model_Ng_Age, type="emm", terms="Repro_time", show.data = F)+
  theme_classic()+
  xlab("Reproduction period")+
  ylab("log(NYCT+1)")
dev.off()


## 3.4 POLY prevalence ####
### 3.4.1 Model selection ####
model <- glmmTMB(PREV_POLY ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=binomial)
model1 <- glmmTMB(PREV_POLY ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site), data=Dodo_Parasito11, 
                 family=binomial)
anova(model,model1) # Year significant 
model1 <- glmmTMB(PREV_POLY ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model,model1) # Site non-significatif mais on laisse dans le modèle pour controle de la non-indépendance 

model1 <- glmmTMB(PREV_POLY ~ Age + Sex + Sex:Repro + Age:Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                 family=binomial)
anova(model,model1) # Age:Repro not significant

model2 <- glmmTMB(PREV_POLY ~ Age + Sex + Sex:Repro + Repro_time  + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model1,model2) # Age:Sex not significant

model3 <- glmmTMB(PREV_POLY ~ Age + Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model2,model3) # Repro not significant

model4 <- glmmTMB(PREV_POLY ~ Age + Sex + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model3,model4) # Repro_time significant

model4 <- glmmTMB(PREV_POLY ~ Age + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model3,model4) # Sex significant

model4 <- glmmTMB(PREV_POLY ~ Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11, 
                  family=binomial)
anova(model3,model4) # Age significant


### 3.4.2 Best model ####

pairs(emmeans(model3,"Sex"), adjust="tukey") #Females more parasitized than males
pairs(emmeans(model3,"Repro_time"), adjust="tukey") #PL < L et P
pairs(emmeans(model3,"Age"), adjust="tukey") # SAD > AD

## 3.5 POLY intensity ####
### 3.5.1 Clean dataset ####

# Dodo_Parasito11b <- Dodo_Parasito11_int[!is.na(Dodo_Parasito11$PREV_POLY),]
# table(Dodo_Parasito11b$Site,Dodo_Parasito11b$Year)

Dodo_Parasito11_int <- Dodo_Parasito11[Dodo_Parasito11$PREV_POLY == 1,]
Dodo_Parasito11_int <- Dodo_Parasito11_int[!is.na(Dodo_Parasito11_int$PREV_POLY),]
table(Dodo_Parasito11_int$Site,Dodo_Parasito11_int$Year)

### 3.5.2 Model selection ####

model <- glmmTMB(log(INF_POLY) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                 family=gaussian)
model1 <- glmmTMB(log(INF_POLY) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Site), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model,model1) # Year non-significatif mais on laisse dans le modèle pour controle de la non-indépendance  
model1 <- glmmTMB(log(INF_POLY) ~ Age + Sex + Sex:Repro + Age:Sex + Age:Sex:Repro + Repro_time + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model,model1) # Site non-significatif mais on laisse dans le modèle pour controle de la non-indépendance 

model1 <- glmmTMB(log(INF_POLY) ~ Age + Sex + Sex:Repro + Age:Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model,model1) # Age:Repro not significant

model2 <- glmmTMB(log(INF_POLY)  ~ Age + Sex + Sex:Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model1,model2) # Age:Sex not significant

model3 <- glmmTMB(log(INF_POLY)  ~ Age + Sex + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model2,model3) # Repro not significant

model4 <- glmmTMB(log(INF_POLY)  ~ Age + Sex + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model3,model4) # Repro_time not significant

model5 <- glmmTMB(log(INF_POLY)  ~ Age + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model4,model5) # Sex not significant

model5 <- glmmTMB(log(INF_POLY)  ~ Sex + (1|Site) + (1|Year), data=Dodo_Parasito11_int, 
                  family=gaussian)
anova(model4,model5) # Age not significant

### Rien de significatif pour l'intensité 




## 3.6 Vizualization of the data ####
### 3.6.1 Create datasets ####
#### 3.6.1.1 Dataset Sex ####

abondance_sex <- Dodo_Parasito11 %>%
  group_by(Sex) %>%
  summarise(
    mean_Spx = mean(log(Spx_tot + 1), na.rm = TRUE),
    se_Spx = sd(log(Spx_tot + 1), na.rm = TRUE) / sqrt(sum(!is.na(Spx_tot))),    
    mean_Nyct = mean(log(Nyct_tot + 1), na.rm = TRUE),
    se_Nyct = sd(log(Nyct_tot + 1), na.rm = TRUE) / sqrt(sum(!is.na(Nyct_tot))),
    mean_ppoly = mean(PREV_POLY, na.rm = TRUE)*100,
    se_ppoly = sd(PREV_POLY, na.rm = TRUE) *100 / sqrt(sum(!is.na(PREV_POLY))),
    n = n()
  )
abondance_sex1 <- Dodo_Parasito11_int %>%
  group_by(Sex) %>%
  summarise(
    mean_ipoly = mean(log(INF_POLY), na.rm = TRUE),
    se_ipoly = sd(log(INF_POLY), na.rm = TRUE) / sqrt(sum(!is.na(INF_POLY))),    
  )
abondance_combined <- left_join(abondance_sex, abondance_sex1, by = "Sex")

abondance_long <- abondance_combined %>%
  pivot_longer(
    cols = c(mean_Spx, mean_Nyct, mean_ppoly, mean_ipoly,
             se_Spx, se_Nyct, se_ppoly, se_ipoly),
    names_to = c(".value", "Variable"),  # crée deux colonnes : "mean" et "se"
    names_pattern = "(mean|se)_(.*)"
  )
abondance_long <- abondance_long %>%
  mutate(Sex = recode(Sex, "M" = "Male", "F" = "Female"))

abondance_long <- abondance_long %>%
  mutate(Variable = factor(Variable, levels = c("ipoly","Nyct", "ppoly","Spx"),
                           labels = c(
                             "LOG(italic(Polychromophilus)~intensity)",
                             "LOG(italic(Nycteribiidae) + 1)",
                             "italic(Polychromophilus)~prevalence~('%')",
                             "LOG(italic(Spinturnix) + 1)"
                           )))


#### 3.6.1.2 Dataset Repro status ####

abondance_repro <- Dodo_Parasito11 %>%
  group_by(Repro) %>%
  summarise(
    mean_Spx = mean(log(Spx_tot + 1), na.rm = TRUE),
    se_Spx = sd(log(Spx_tot + 1), na.rm = TRUE) / sqrt(sum(!is.na(Spx_tot))),
    n = n()
  )
abondance_repro <- abondance_repro %>%
  mutate(Sex = ifelse(Repro %in% c("R", "NR"), "Male", "Female"))

abondance_repro <- abondance_repro %>%
  mutate(Repro = factor(Repro, levels = c("NR","R", "O", "P","L", "PL"), labels = c("Non-Reproductive","Reproductive","Other","Pregnant","Lactating", "Post-Lactating")))


#### 3.6.1.3 Dataset Age ####

abondance_age <- Dodo_Parasito11 %>%
  group_by(Age) %>%
  summarise(
    mean_PolyP = mean(log(PREV_POLY + 1), na.rm = TRUE),
    se_PolyP = sd(log(PREV_POLY + 1), na.rm = TRUE) / sqrt(sum(!is.na(PREV_POLY))),
    n = n()
  )
abondance_age <- abondance_age %>%
  mutate(Age = recode(Age, "AD" = "Adult", "SAD" = "Subadult", "JUV" = "Juvenile"))




### 3.6.2 Effect of sex for each parasites ####

t <- ggplot(abondance_long, aes(x = Sex, y = mean, fill = Variable)) +
  geom_pointrange(aes(ymin = mean - 1.94 * se, ymax = mean + 1.94 * se)) + 
  facet_wrap(~ Variable, scales = "free_y",     labeller = label_parsed) +
  labs(y = "Mean (± 95% CI)", x = "Sex") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("Fig_sex_effect.pdf",plot = t,dpi = "retina",width = 6, height = 6)


### 3.6.3 Effect of reproductive status on SPX ####

u <- ggplot(abondance_repro, aes(x = Repro, y = mean_Spx, color = Sex)) +
  geom_pointrange(aes(ymin = mean_Spx - 1.94 * se_Spx, ymax = mean_Spx + 1.94 * se_Spx)) + 
  labs(x = "Reproductive status", 
       y = expression("LOG("*italic("Spinturnix")*" + 1)"),
       color = "Sex") +
  scale_color_manual(values = c("Male" = "black", "Female" = "grey")) +
  theme_bw()
u

ggsave("Fig_repro_effect.pdf",plot = u,dpi = "retina",width = 7.5, height = 3)


### 3.6.4 Effect of age for each parasites on Poly prev ####

v <- ggplot(abondance_age, aes(x = Age, y = mean_PolyP, color = Age)) +
  geom_pointrange(aes(ymin = mean_PolyP - 1.94 * se_PolyP, ymax = mean_PolyP + 1.94 * se_PolyP)) + 
  labs(x = "Age", 
       y = expression(""*italic("P. murinus")*" presence"),
       color = "Age") +
  scale_color_grey() +
  theme_bw()
v

ggsave("Fig_age_effect.pdf",plot = v,dpi = "retina",width = 5, height = 3)



#######################
### SEM
######################



Dodo_Parasito_sem <- Dodo_Parasito11[!is.na(Dodo_Parasito11$PREV_POLY),]
Dodo_Parasito_sem$LOG_Nyct_tot <- log(Dodo_Parasito_sem$Nyct_tot + 1)
Dodo_Parasito_sem$LOG_Spx_tot <- log(Dodo_Parasito_sem$Spx_tot + 1)
Dodo_Parasito_sem$LOG_INF_POLY<- log(Dodo_Parasito_sem$INF_POLY + 1)

### AVEC PREVALENCE
# modèle de base: on compare tous les modèles entre eux
sem1 = psem(
  glmmTMB(PREV_POLY ~ Age + Repro + Repro_time + LOG_Nyct_tot + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=binomial),
#  glmmTMB(LOG_INF_POLY ~ Age + Repro + Repro_time + LOG_Nyct_tot + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  glmmTMB(LOG_Nyct_tot ~ Age + Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  glmmTMB(LOG_Spx_tot ~ Age + Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  LOG_Spx_tot %~~% LOG_Nyct_tot
)
summary(sem1)
dSep(sem1)


# sem1 dit que y a un lien entre spx~prevPOLY donc on reconstruit
sem2 = psem(
  glmmTMB(PREV_POLY ~ Age + Repro + Repro_time + LOG_Nyct_tot + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=binomial),
  #  glmmTMB(LOG_INF_POLY ~ Age + Repro + Repro_time + LOG_Nyct_tot + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  glmmTMB(LOG_Nyct_tot ~ Age + Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  glmmTMB(LOG_Spx_tot ~ Age + Repro + Repro_time + PREV_POLY + (1|Site) + (1|Year), data=Dodo_Parasito_sem, family=gaussian),
  LOG_Spx_tot %~~% LOG_Nyct_tot
)
summary(sem2)

plot(sem1)


### MODELE AVEC INTENSITY
Dodo_Parasito_sem_int <- Dodo_Parasito11_int[!is.na(Dodo_Parasito11_int$INF_POLY),]
Dodo_Parasito_sem_int$LOG_Nyct_tot <- log(Dodo_Parasito_sem_int$Nyct_tot + 1)
Dodo_Parasito_sem_int$LOG_Spx_tot <- log(Dodo_Parasito_sem_int$Spx_tot + 1)
Dodo_Parasito_sem_int$LOG_INF_POLY<- log(Dodo_Parasito_sem_int$INF_POLY + 1)

sem2 = psem(
  glmmTMB(LOG_INF_POLY ~ Age + Repro + Repro_time + LOG_Nyct_tot + (1|Site) + (1|Year), data=Dodo_Parasito_sem_int, family=gaussian),
  glmmTMB(LOG_Nyct_tot ~ Age + Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito_sem_int, family=gaussian),
  glmmTMB(LOG_Spx_tot ~ Age + Repro + Repro_time + (1|Site) + (1|Year), data=Dodo_Parasito_sem_int, family=gaussian),
  LOG_Spx_tot %~~% LOG_Nyct_tot
)
summary(sem2)



################################################################################
# 4. Bat effects - NYCT ####
################################################################################
## 4.1 All ages ####
### 4.1.1 Family selection ####

model_Np1 <- glmmTMB(Nyct_tot ~ Age + Sex + Age:Sex + (1|Site) + (1|Year), data=Dodo_Parasito1, 
                    family=poisson)
model_Ng1 <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Age:Sex + (1|Site) + (1|Year), data=Dodo_Parasito1, 
                    family=gaussian)
AIC(model_Np1,model_Ng1)
#Best model: gaussian 


### 4.1.2 Check model assumptions ####

hist(residuals(model_Ng1)) # normalité ok
plot(predict(model_Ng1, type = "response"),resid(model_Ng1, type = "pearson")) # homogénéité ok
check_model(model_Ng1)


### 4.1.3 Model selection ####

model_Ng1_site <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Age:Sex + (1|Year), data=Dodo_Parasito1, 
                         family=gaussian)
anova(model_Ng1,model_Ng1_site) #Site significant

model_Ng1_year <- glmmTMB(log(Nyct_tot + 1) ~Age + Sex + Age:Sex + (1|Site), data=Dodo_Parasito1, 
                         family=gaussian)
anova(model_Ng1,model_Ng1_year) #Year significant

model_Ng1_AS <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + (1|Site) + (1|Year), data=Dodo_Parasito1, 
                     family=gaussian)
anova(model_Ng1,model_Ng1_AS) #Interaction not significant

model_Ng1_Age<- glmmTMB(log(Nyct_tot + 1) ~ Sex + (1|Site) + (1|Year), data=Dodo_Parasito1, 
                        family=gaussian)
anova(model_Ng1_AS, model_Ng1_Age) #Age significant

model_Ng1_Sex<- glmmTMB(log(Nyct_tot + 1) ~ Age + (1|Site) + (1|Year), data=Dodo_Parasito1, 
                        family=gaussian)
anova(model_Ng1_AS, model_Ng1_Sex) #Sex significant



### 4.1.4 Vizualisation of the data ####


pairs(emmeans(model_Ng1_AS,"Sex"), adjust="tukey") #Females more parasitized than males
pairs(emmeans(model_Ng1_AS,"Age"), adjust="tukey") #Juv>SAD, Juv>AD, AD=SAD

svg("./Nyct_Sex_allAge.svg")
plot_model(model_Ng1_AS, type="emm", terms="Sex", show.data = F)+
  theme_classic()+
  xlab("Sex")+
  ylab("log(NYCT+1)")
dev.off()

svg("./Nyct_Age_allAge.svg")
plot_model(model_Ng1_AS, type="emm", terms="Age", show.data = F)+
  theme_classic()+
  xlab("Age")+
  ylab("log(NYCT+1)")
dev.off()


## 4.2 Only adults ####
### 4.2.1 Dataset only adults ####

#Remove juveniles to test Repro
Dodo_ParasitoAD <- Dodo_Parasito1[Dodo_Parasito1$Age != "JUV",]


### 4.2.2 Family selection ####

model_Np2 <- glmmTMB(Nyct_tot ~ Age + Sex + Repro + Age:Sex + Age:Repro + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                     family=poisson)
model_Ng2 <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Repro + Age:Sex + Age:Repro + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                     family=gaussian)
AIC(model_Np2,model_Ng2)
#Best model: gaussian 


### 4.2.3 Check model assumptions ####

hist(residuals(model_Ng2)) # normalité ok
plot(predict(model_Ng2, type = "response"),resid(model_Ng2, type = "pearson")) # homogénéité ok
check_model(model_Ng2)


### 4.2.4 Model selection ####

model_Ng2_site <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Repro + Age:Sex + Age:Repro + (1|Year), data=Dodo_ParasitoAD, 
                          family=gaussian)
anova(model_Ng2,model_Ng2_site) #Site significant

model_Ng2_year <- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Repro + Age:Sex + Age:Repro + (1|Site), data=Dodo_ParasitoAD, 
                          family=gaussian)
anova(model_Ng2,model_Ng2_year) #Year significant

model_Ng2_AS <- glmmTMB(log(Nyct_tot + 1) ~  Age + Sex + Repro + Age:Repro + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                        family=gaussian)
anova(model_Ng2,model_Ng2_AS) #Interaction age and sex significant

model_Ng2_AR <- glmmTMB(log(Nyct_tot + 1) ~  Age + Sex + Repro + Age:Sex + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                        family=gaussian)
anova(model_Ng2,model_Ng2_AR) #Interaction age and repro not significant

model_Ng2_Age<- glmmTMB(log(Nyct_tot + 1) ~ Sex + Repro + Age:Sex + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                        family=gaussian)
anova(model_Ng2_AR, model_Ng2_Age) #Age significant

model_Ng2_Sex<- glmmTMB(log(Nyct_tot + 1) ~ Age + Repro + Age:Sex + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                        family=gaussian)
anova(model_Ng2_AR, model_Ng2_Sex) #Sex significant

model_Ng2_repro<- glmmTMB(log(Nyct_tot + 1) ~ Age + Sex + Age:Sex + (1|Site) + (1|Year), data=Dodo_ParasitoAD, 
                        family=gaussian)
anova(model_Ng2_AR, model_Ng2_repro) #Repro significant


### 4.1.4 Vizualisation of the data ####

# 
# pairs(emmeans(model_Ng2_AR,"Sex"), adjust="tukey") #Females more parasitized than males
# pairs(emmeans(model_Ng2_AR,"Age"), adjust="tukey") #Juv>SAD, Juv>AD, AD=SAD
# pairs(emmeans(model_Ng2_AR,"Age"), adjust="tukey") #Juv>SAD, Juv>AD, AD=SAD
# pairs(emmeans(model_Ng2_AR,"Repro"), adjust="tukey") #Juv>SAD, Juv>AD, AD=SAD
# 
# svg("./Nyct_Sex_allAge.svg")
# plot_model(model_Ng1_AS, type="emm", terms="Sex", show.data = F)+
#   theme_classic()+
#   xlab("Sex")+
#   ylab("log(NYCT+1)")
# dev.off()
# 
# svg("./Nyct_Age_allAge.svg")
# plot_model(model_Ng1_AS, type="emm", terms="Age", show.data = F)+
#   theme_classic()+
#   xlab("Age")+
#   ylab("log(NYCT+1)")
# dev.off()
