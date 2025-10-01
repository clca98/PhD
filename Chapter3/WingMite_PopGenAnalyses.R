#######################################################
##### Wing mites - Analyses - Population Genetics #####
#######################################################

setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/5.Thesis/Project/7.RADseq/Project_RAD2024/Parasites/Spx/New_sequencing/5.Analysis/minDP10/Nomisdata/")


################################################################################
# 1. Download packages and data ####
################################################################################
## 1.1 Download packages ####

library(hierfstat)
library(gaston)
library(vcfR)
library(ade4)
library(reshape2)
library(lme4)
library(car)
library(ggplot2)
library(pheatmap)
library(viridis)
library(fields)
library(vegan)
library(Matrix)
library(plyr)
library(dplyr)
library(tidyr)
library(gplots)
library(calibrate)
library(Rmisc)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(SNPRelate)

## 1.2 Download data ####

SPX_m3M5_minDP10 <- read.VCF("./SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.recode.vcf")
dim(SPX_m3M5_minDP10) #171 9525


################################################################################
# 2. Exploring data ####
################################################################################
## 2.1 Distribution of the genotypes ####

plot(SPX_m3M5_minDP10@ped$N1,SPX_m3M5_minDP10@ped$N2, pch=20)

## 2.2 Heterozygosity per individuals ####

Hetz <- read.table("./SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.het", header=TRUE, sep="\t")
Hetz$HETZ <- (Hetz$N_SITES-Hetz$O.HOM.)/Hetz$N_SITES
hist(Hetz$HETZ, breaks=100)

## 2.3 Coverage per individuals ####

Cov <- read.table("./SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.idepth", header=TRUE, sep="\t")
hist(Cov$MEAN_DEPTH, breaks=100)
mean(Cov$MEAN_DEPTH)

## 2.4 Heterozity Coverage correlation ####

plot(Hetz$HETZ~Cov$MEAN_DEPTH)
mod1 <- lm(Hetz$HETZ~Cov$MEAN_DEPTH)
summary(mod1)# R2 = 0.009999 ; p = 0.1011
abline(mod1, col="grey45", lwd = 2) #Control

cor.test(Hetz$HETZ,Cov$MEAN_DEPTH)#R2 = 0.1257872; p = 0.1011


################################################################################
# 3. Creating datasets ####
################################################################################
## 3.1 Including host, sites, sex and host sex ####

popmap <- read.csv("../../Popmap_sexes.csv", header =TRUE, sep=";")

new_popmap <- popmap[popmap$ID %in% SPX_m3M5_minDP10@ped$id,]
poplist <- new_popmap$Site
spxsex <- new_popmap$Spx_sex
batsex <- new_popmap$Bat_sex
IDlist <- new_popmap$ID
Library <- new_popmap$Library 

table(new_popmap$Site)
table(new_popmap$Library)


## 3.2 Remove kin individuals k = 0.125 ####

SPX_prune125 <- read.table("./betaSPX_prune0125.table")

## 3.3 SNPs dataset ####

genotype_matrix <- as.matrix(SPX_m3M5_minDP10)
SPX <- data.frame(ID=IDlist, pop=poplist, spxsex=spxsex, batsex=batsex, genotype_matrix)

#Sort by population
SPX_popsorted <- SPX[order(SPX$pop),]
SPX_popsorted_matrix <- as.matrix(SPX_popsorted[,5:ncol(SPX_popsorted)])
SPX_popsorted_matrix <- apply(SPX_popsorted_matrix, 2, as.numeric)

#Without kin individuals
SPX_nokin <- SPX_popsorted[SPX_popsorted$ID %in% rownames(SPX_prune125),]
SPX_nokin_matrix <- as.matrix(SPX_nokin[,5:ncol(SPX_nokin)])
SPX_nokin_matrix <- apply(SPX_nokin_matrix, 2, as.numeric)


## 3.4 Genotypes dataset ####

Genotype_SPX <- SPX[,-c(1:4)]
Genotype_SPX[Genotype_SPX==0] <- 11
Genotype_SPX[Genotype_SPX==1] <- 12
Genotype_SPX[Genotype_SPX==2] <- 22
Genotype_SPX <- data.frame(ID=SPX[,1],pop=SPX[,2],spxsex=SPX[,3], batsex=SPX[,4], Genotype_SPX)

#Without kin individuals
SPX_nokin_geno <- SPX_nokin[,-c(1:4)]
SPX_nokin_geno[SPX_nokin_geno==0] <- 11
SPX_nokin_geno[SPX_nokin_geno==1] <- 12
SPX_nokin_geno[SPX_nokin_geno==2] <- 22
SPX_nokin_geno <- data.frame(ID=SPX_nokin[,1],pop=SPX_nokin[,2],nyctsex=SPX_nokin[,3], batsex=SPX_nokin[,4], SPX_nokin_geno)
SPX_nokin_geno_matrix <- apply(as.matrix(SPX_nokin_geno[,5:ncol(SPX_nokin_geno)]), 2, as.numeric)


################################################################################
# 4. Genetic Analyses - overall dataset ####
################################################################################
## 4.1 Genetic diversity ####

# per pop
HetZ_pop <- data.frame(ID=IDlist, pop=poplist, spxsex=spxsex, batsex=batsex, Hetz)
boxplot(HetZ_pop$HETZ~HetZ_pop$pop)

test <- kruskal.test(HetZ_pop$HETZ~HetZ_pop$pop)
result <- TukeyHSD(test)
which(result$`HetZ_pop$pop`[,4] < 0.05)#Moudon-LaCascade --> effet taille de pop


mean_het_per_pop <- aggregate(HETZ ~ pop, data = HetZ_pop, FUN = mean)

# Heterozygosity per individuals
DATA = snpgdsVCF2GDS("./SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.recode.vcf", "SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.gds")
DATA = snpgdsOpen("SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.gds")
gen_mat <- snpgdsGetGeno('SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.gds')
IndHet = rowSums(gen_mat==1, na.rm = TRUE) / rowSums(! is.na(gen_mat))

HetZ_pop2 <- data.frame(ID=IDlist, pop=poplist, IndHet)
boxplot(HetZ_pop2$IndHet~HetZ_pop2$pop)

test <- kruskal.test(HetZ_pop2$IndHet~HetZ_pop2$pop)#Kruskal-Wallis chi-squared = 23.4, df = 13, p-value = 0.03712
result <- TukeyHSD(test)
which(result$`HetZ_pop$pop`[,4] < 0.05)#Moudon-LaCascade --> effet taille de pop


#Without kin individuals
SPXHetZ_nokin <- HetZ_pop2[HetZ_pop2$ID %in% rownames(SPX_prune125),]
test2 <- kruskal.test(SPXHetZ_nokin$IndHet~SPXHetZ_nokin$pop)#Kruskal-Wallis chi-squared = 20.98, df = 13, p-value = 0.07333
mean_het_per_pop <- aggregate(IndHet ~ pop, data = SPXHetZ_nokin, FUN = mean)


## 4.2 PCA per site ####
### 4.2.1 On the entire dataset ####

pca <- indpca(SPX_nokin_geno[,c(1,5:ncol(SPX_nokin_geno))])
plot(pca)

# % of genetic variance that explain each axis 
pca_eig <- 100*pca$ipca$eig/sum(pca$ipca$eig) #PC1 = 1.0322726% & PC2 = 0.8582099% | nokin: PC1=0.8908002 ; PC2=0.8849327
barplot(pca_eig[1:10])
sum(pca_eig[1:2])#1.890482
sum(pca_eig[1:10])#8.496456


pca_ALL <- cbind(SPX_nokin_geno[,1:4], pca$ipca$li[, 1:2])
pca_sortedpop <- pca_ALL[order(pca_ALL$pop),]

# Obtenir les couleurs pour chaque population
pop_levels <- unique(pca_sortedpop$pop)
col_pop <- viridis(length(pop_levels))
names(col_pop) <- pop_levels  # Associe chaque couleur à une population
pch_pop <- 1:length(pop_levels)
names(pch_pop) <- pop_levels

# Appliquer les couleurs en fonction de chaque individu
ind_colors <- col_pop[as.character(pca_sortedpop$pop)]
ind_pch <- pch_pop[as.character(pca_sortedpop$pop)]

svg("PCA_site_SPX_minDP10_nomissingdata_nokin.svg", bg = "transparent")
plot(pca_sortedpop[,5], pca_sortedpop[,6], col=ind_colors, pch=ind_pch, lwd=1,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-8, -4, 0, 4, 8))
axis (side=1, at=c(-8, -4, 0, 4, 8))
title(xlab="PC1 (0.891% explained variation)", 
      ylab="PC2 (0.885% explained variation)", cex.lab=1.5)
dev.off()
svg("PCA_site_SPX_minDP10_nomissingdata_nokin_legend.svg", bg = "transparent")
plot.new()
legend("topright", legend=1:14, col=col_pop, pch = pch_pop, title = "Pop")
dev.off()

### 4.2.2 Without outliers ####

# 2 outliers: 535S2 and 544S2
# Hetz = 0.1724680 & 0.1704507 (mean Het = 0.1569139)
# New PCA w/out them

pca_sortedpop2 <- pca_sortedpop[!pca_sortedpop$ID %in% c("535S2_mergedALL", "544S2_mergedALL"),]
pop_levels2 <- unique(pca_sortedpop2$pop)
col_pop2 <- viridis(length(pop_levels2))
names(col_pop2) <- pop_levels2  # Associe chaque couleur à une population
pch_pop2 <- 1:length(pop_levels2)
names(pch_pop2) <- pop_levels2

# Appliquer les couleurs en fonction de chaque individu
ind_colors2 <- col_pop2[as.character(pca_sortedpop2$pop)]
ind_pch2 <- pch_pop2[as.character(pca_sortedpop2$pop)]

svg("PCA_site_SPX_minDP10_nomissingdata2.svg", bg = "transparent")
plot(pca_sortedpop2[,5], pca_sortedpop2[,6], col=ind_colors2, pch=ind_pch2, lwd=1,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-6, -4, -2, 0, 2, 4, 6, 8, 10))
axis (side=1, at=c(-2, 0, 2, 4, 6))
title(xlab="PC1 (1.032% explained variation)", 
      ylab="PC2 (0.858% explained variation)", cex.lab=1.5)
legend("bottomright", legend=1:14, col=col_pop2, pch = pch_pop2, title = "Pop")
dev.off()


## 4.3 Distribution of the rare alleles per sampling sites ####

SPX_alleles <- SPX[,5:ncol(SPX_popsorted)]
mac <- colSums(SPX_alleles, na.rm = TRUE)
mac_rare <- mac[mac<=10]
loci_Rare <- attributes(mac_rare)$names

Present_in_Npop <- c() #vector that tells number of pop the allele is present in
for (i in 1:length(loci_Rare)){
  locuscol <- match(loci_Rare[i],names(SPX))
  Which_pops <- unique(SPX$pop[SPX[,locuscol]>0])
  Which_pops <- Which_pops[!is.na(Which_pops)]
  Present_in_Npop <- c(Present_in_Npop, length(Which_pops))
}

hist(
  Present_in_Npop,
  breaks=seq(0.5, max(Present_in_Npop)+0.5, by=1),
  col = "lightblue",
  border = "white",
  main = "Rare alleles",
  xlab = "Number of sites",
  ylab = "Count",
  xlim = c(min(Present_in_Npop, na.rm = TRUE), max(Present_in_Npop, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1,
  xaxt="n"
)
axis(1, at=1:10, labels=1:10)


## 4.4 Hierarchical Analysis of Variance - RARE ALLELES ####
### 4.4.1 Fs.dosage on rare alleles ####

SPX_rareall <- cbind(SPX[,c(1,2)],SPX[,names(mac_rare)])

HierVar_rareall <- fs.dosage(SPX_rareall[,c(3:3446)], pop=SPX_rareall[,2])

svg("./Fsdosage_rareAlleles_SPX.svg")
plot(HierVar_rareall)
dev.off()

par(mfrow=c(1,1))


### 4.4.2 Plot Pairwise kinship per pop - FsM - rare alleles ####

kinship_rareall_SPX <- HierVar_rareall$FsM
kinship_rareall_SPX[lower.tri(kinship_rareall_SPX)] <- NA

colFsM <- colorRampPalette(brewer.pal(9, "PuBu"))(100)

svg("./Kinship_RareAll_SPX.svg")
heatmap(kinship_rareall_SPX, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFsM, scale = "none", na.rm = TRUE, zlim = range(kinship_rareall_SPX, na.rm = TRUE))
dev.off()
svg("./Kinship_RareAll_SPX_scale.svg")
image.plot(legend.only = TRUE, zlim = range(kinship_rareall_SPX, na.rm = TRUE), 
           col = colFsM, legend.line = 2)
dev.off()



### 4.4.3 Plot Pairwise Fst per pop - Fs2x2 - rare alleles ####


Pfst_RareAll_SPX <- HierVar_rareall$Fst2x2
Pfst_RareAll_SPX[upper.tri(Pfst_RareAll_SPX)] <- NA
diag(Pfst_RareAll_SPX) <- NA

colFs2x2 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

svg("./Fst_RareAll_SPX.svg")
heatmap(Pfst_RareAll_SPX, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFs2x2, scale = "none", na.rm = TRUE, zlim = range(Pfst_RareAll_SPX, na.rm = TRUE))
dev.off()
svg("./Fst_RareAll_SPX_scale.svg")
image.plot(legend.only = TRUE, zlim = range(Pfst_RareAll_SPX, na.rm = TRUE), 
           col = colFs2x2, legend.line = 2)
dev.off()




### 4.4.4 Bootstrapping For INBREEDING COEFF & FST PER POP - rare alleles ####


# Number of bootstraps
n_bootstraps <- 100
n_loci <- ncol(SPX_rareall[, 3:ncol(SPX_rareall)])

# Get unique population names
populations <- unique(SPX_rareall[, 2])
nb_pop <- length(populations)

# Initialize empty lists to store bootstrap results
results_Fs_RareAll <- data.frame(
  bootstrap = integer(0),
  site = character(0),
  Fis_value = numeric(0),
  Fst_value = numeric(0),
  stringsAsFactors = FALSE
)

for (i in 1:n_bootstraps) {
  sampled_loci <- sample(3:ncol(SPX_rareall), n_loci, replace = TRUE)
  vcf_sample <- SPX_rareall[, c(2, sampled_loci)]
  
  HierVar_boot <- fs.dosage(vcf_sample[, -1], pop = vcf_sample[, 1])
  
  if (is.null(HierVar_boot)) next

  for (pop_name in names(HierVar_boot$Fs[1, ])) {
    results_Fs_RareAll <- rbind(results_Fs_RareAll, data.frame(
      bootstrap = i,
      site = pop_name,
      Fis_value = HierVar_boot$Fs[1, pop_name],
      Fst_value = HierVar_boot$Fs[2, pop_name],
      stringsAsFactors = FALSE
    ))
  }
}

# Boxplot visualization

svg("./SPX_FISbootstrap_RareAll.svg")
boxplot(Fis_value ~ site, data = results_Fs_RareAll,
        main = "Mean FIS per population in Rare alleles",
        xlab = "", 
        ylab = "Mean FIS",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()

svg("./SPX_FSTbootstrap_RareAll.svg")
boxplot(Fst_value ~ site, data = results_Fs_RareAll,
        main = "Mean FST per population in Rare alleles",
        xlab = "Population", 
        ylab = "Mean FST",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()


## 4.5 Hierarchical Analysis of Variance - OVERALL dataset ####
### 4.5.1 Fs.dosage ####

SPX_popsorted_noout <- SPX_popsorted[SPX_popsorted$ID != "498S1_mergedALL" & SPX_popsorted$ID != "792S1_mergedALL",]


HierVar_SPX <- fs.dosage(SPX_popsorted[,c(5:ncol(SPX_popsorted))], pop=SPX_popsorted[,2])

svg("./Fsdosage_noCascade_SPX.svg")
plot(HierVar_SPX)
dev.off()

par(mfrow=c(1,1))


#Without kin individuals
HierVar_SPXnokin <- fs.dosage(SPX_nokin[,c(5:ncol(SPX_nokin))], pop=SPX_nokin[,2])



## Fs.dosage on NYCT to have the same ranges

global_zlim_FsM <- range(HierVar_SPXnokin$FsM, HierVar_Nyctnokin$FsM, na.rm = TRUE)
global_zlim_Fst2x2 <- range(HierVar_SPXnokin$Fst2x2, HierVar_Nyctnokin$Fst2x2, na.rm = TRUE)


### 4.5.2 Plot Pairwise kinship per pop - FsM ####

kinship_mx_SPXnokin <- HierVar_SPXnokin$FsM
kinship_mx_SPXnokin[lower.tri(kinship_mx_SPXnokin)] <- NA

colFsM <- colorRampPalette(brewer.pal(9, "PuBu"))(100)

svg("./Kinship_Blue_SPXnokin.svg")
heatmap(kinship_mx_SPXnokin, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFsM, scale = "none", na.rm = TRUE, zlim = global_zlim_FsM)
dev.off()
svg("./Kinship_blue_SPXnokin_scale.svg")
image.plot(legend.only = TRUE, zlim = global_zlim_FsM, 
           col = colFsM, legend.line = 2)
dev.off()


### 4.5.3 Plot Pairwise Fst per pop - Fs2x2 ####

Pfst_mx_SPXnokin <- HierVar_SPXnokin$Fst2x2
Pfst_mx_SPXnokin[upper.tri(Pfst_mx_SPXnokin)] <- NA
diag(Pfst_mx_SPXnokin) <- NA

colFs2x2 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

svg("./Fst_red_SPXnokin.svg")
heatmap(Pfst_mx_SPXnokin, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFs2x2, scale = "none", na.rm = TRUE, zlim = global_zlim_Fst2x2)
dev.off()
svg("./Fst_red_SPXnokin_scale.svg")
image.plot(legend.only = TRUE, zlim = global_zlim_Fst2x2, 
           col = colFs2x2, legend.line = 2)
dev.off()


### 4.5.4 Bootstrapping For INBREEDING COEFF & FST PER POP ####

# Number of bootstraps
n_bootstraps <- 100
n_loci <- ncol(SPX_nokin[, 5:ncol(SPX_nokin)])

# Get unique population names
populations <- unique(SPX_nokin[, 2])
nb_pop <- length(populations)

# Initialize empty lists to store bootstrap results
results_Fs_SPXnokin <- data.frame(
  bootstrap = integer(0),
  site = character(0),
  Fis_value = numeric(0),
  Fst_value = numeric(0),
  stringsAsFactors = FALSE
)

# Store hierarchical variance components
FsMSPXnokin_bootstrap <- list()
Fs2x2SPXnokin_bootstrap <- list()


for (i in 1:n_bootstraps) {
  sampled_loci <- sample(5:ncol(SPX_nokin), n_loci, replace = TRUE)
  vcf_sample <- SPX_nokin[, c(2, sampled_loci)]
  
  HierVar_bootSPXnokin <- fs.dosage(vcf_sample[, -1], pop = vcf_sample[, 1])
  
  if (is.null(HierVar_bootSPXnokin)) next
  
  FsMSPXnokin_bootstrap[[i]] <- HierVar_bootSPXnokin$FsM
  Fs2x2SPXnokin_bootstrap[[i]] <- HierVar_bootSPXnokin$Fst2x2
  
  for (pop_name in names(HierVar_bootSPXnokin$Fs[1, ])) {
    results_Fs_SPXnokin <- rbind(results_Fs_SPXnokin, data.frame(
      bootstrap = i,
      site = pop_name,
      Fis_value = HierVar_bootSPXnokin$Fs[1, pop_name],
      Fst_value = HierVar_bootSPXnokin$Fs[2, pop_name],
      stringsAsFactors = FALSE
    ))
  }
}

# Compute confidence intervals
ci_95 <- function(values) {
  quantile(values, probs = c(0.025, 0.975), na.rm = TRUE)
}

# Compute CI for Fis and Fst values per population
Fis_CISPXnokin <- aggregate(Fis_value ~ site, data = results_Fs_SPXnokin, ci_95)
Fst_CISPXnokin <- aggregate(Fst_value ~ site, data = results_Fs_SPXnokin, ci_95)

# Compute CI for each entry in FsM and Fs2x2 matrices
FsM_matrixSPXnokin <- simplify2array(FsMSPXnokin_bootstrap)
Fs2x2_matrixSPXnokin <- simplify2array(Fs2x2SPXnokin_bootstrap)

FsM_CI_matrixSPXnokin <- apply(FsM_matrixSPXnokin, c(1,2), ci_95)
Fs2x2_CI_matrixSPXnokin <- apply(Fs2x2_matrixSPXnokin, c(1,2), ci_95)

# Boxplot visualization

svg("./SPX_FISbootstrap_nokin.svg")
boxplot(Fis_value ~ site, data = results_Fs_SPXnokin,
        main = "Mean FIS per population",
        xlab = "", 
        ylab = "Mean FIS",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()

svg("./SPX_FSTbootstrap_nokin.svg")
boxplot(Fst_value ~ site, data = results_Fs_SPXnokin,
        main = "Mean FST per population",
        xlab = "Population", 
        ylab = "Mean FST",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()

# Output confidence intervals
write.csv(Fis_CISPXnokin, "SPX_Fis_CI.csv", row.names = FALSE)
write.csv(Fst_CISPXnokin, "SPX_Fst_CI.csv", row.names = FALSE)


## 4.6 Pairwise kinship per individuals ####
### 4.6.1 Kinship matrix ####

kinship <- beta.dosage(SPX_popsorted_matrix)
diag(kinship) = NA

#Without kin individuals
kinship_SPXnokin <- beta.dosage(SPX_nokin_matrix)
diag(kinship_SPXnokin) = NA



svg("Distrib_kinship_SPXnokin.svg")
hist(
  kinship_SPXnokin,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(kinship_SPXnokin, na.rm = TRUE), max(kinship_SPXnokin, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()


kinship[upper.tri(kinship)] <- NA

svg("./Hist_kinship.svg")
hist(kinship, breaks=100)
dev.off()

image(kinship)

popcount <- table(SPX$pop)
pop_cumul <- cumsum(popcount)

image(1:nrow(kinship), 1:ncol(kinship), kinship, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(kinship, na.rm = TRUE))
abline(v = pop_cumul + 0.5, lwd = 2, col = "white")
plot.new()
image.plot(legend.only = TRUE, zlim = range(kinship, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)

### 4.6.2 Exploration of kinship links ####

## LINKS > 0.2
kinship[upper.tri(kinship)] <- NA
pairs20_SPX <- which(kinship > 0.2, arr.ind = TRUE) #276
pairs25_SPX <- which(kinship > 0.25, arr.ind = TRUE) #30
pairs40_SPX <- which(kinship > 0.4, arr.ind = TRUE) #2

# Extract individual names for the pairs
ind1 <- SPX_popsorted$ID[pairs20_SPX[, 1]]  # First individual in each pair
ind2 <- SPX_popsorted$ID[pairs20_SPX[, 2]]  # Second individual in each pair
kinship_table <- data.frame(Individual1 = ind1, Individual2 = ind2, Kinship = kinship[pairs20_SPX])
kinship_table_order <- kinship_table[order(kinship_table$Individual1),]

# Count pairs per individuals
pair_counts <- table(kinship_table$Individual2)
pair_counts_df <- data.frame(Individual = names(pair_counts), Num_Pairs = as.numeric(pair_counts))

# Extract individual names for the pairs > 0.4
ind1.40 <- SPX_popsorted$ID[pairs40_SPX[, 1]]  # First individual in each pair
ind2.40 <- SPX_popsorted$ID[pairs40_SPX[, 2]]  # Second individual in each pair
kinship_table_40 <- data.frame(Individual1 = ind1.40, Individual2 = ind2.40, Kinship = kinship[pairs40_SPX])


### 4.6.3 Distribution of the genotypes ####

plot(SPX_m3M5_minDP10@ped$N1,SPX_m3M5_minDP10@ped$N2, pch=20)

#Find outliers with a lot of links
which(SPX_m3M5_minDP10@ped$id=="498S1_mergedALL") #67
table(as.matrix(SPX_m3M5_minDP10)[67,]) # 0 (6039); 1 (866); 2 (63)
points(SPX_m3M5_minDP10@ped$N1[67],SPX_m3M5_minDP10@ped$N2[67],col="#33CCCC",pch=17)

which(SPX_m3M5_minDP10@ped$id=="792S2_mergedALL") #167
table(as.matrix(SPX_m3M5_minDP10)[167,]) # 0 (6458); 1 (955); 2 (58)
points(SPX_m3M5_minDP10@ped$N1[167],SPX_m3M5_minDP10@ped$N2[167],col="#ff7f0e",pch=17)


#Find the twins
which(SPX_m3M5_minDP10@ped$id=="477S2_mergedALL") #55
table(as.matrix(SPX_m3M5_minDP10)[55,]) # 0 (6084); 1 (1130); 2 (251)
points(SPX_m3M5_minDP10@ped$N1[55],SPX_m3M5_minDP10@ped$N2[55],col="red",pch=20)

which(SPX_m3M5_minDP10@ped$id=="486S1_mergedALL") #58
table(as.matrix(SPX_m3M5_minDP10)[58,]) # 0 (6510); 1 (1096); 2 (155)
points(SPX_m3M5_minDP10@ped$N1[58],SPX_m3M5_minDP10@ped$N2[58],col="orange",pch=20)

which(SPX_m3M5_minDP10@ped$id=="792S1_mergedALL") #166
table(as.matrix(SPX_m3M5_minDP10)[166,]) # 0 (6773); 1 (1141); 2 (109)
points(SPX_m3M5_minDP10@ped$N1[166],SPX_m3M5_minDP10@ped$N2[166],col="green3",pch=20)

which(SPX_m3M5_minDP10@ped$id=="792S2_mergedALL") #167
table(as.matrix(SPX_m3M5_minDP10)[167,]) # 0 (6458); 1 (955); 2 (58)
points(SPX_m3M5_minDP10@ped$N1[167],SPX_m3M5_minDP10@ped$N2[167],col="lightblue3",pch=20)


### 4.6.4 OUtliers from the PCA ####

which(SPX_popsorted$ID == "535S2_mergedALL") #123
which(SPX_popsorted$ID == "544S2_mergedALL") #128
kinship[128,123] #0.264406

which(SPX_m3M5_minDP10@ped$id=="535S2_mergedALL") #101
table(as.matrix(SPX_m3M5_minDP10)[101,]) # 0 (7124); 1 (1536); 2 (246)
points(SPX_m3M5_minDP10@ped$N1[101],SPX_m3M5_minDP10@ped$N2[101],col="red",pch=17)

which(SPX_m3M5_minDP10@ped$id=="544S2_mergedALL") #106
table(as.matrix(SPX_m3M5_minDP10)[106,]) # 0 (7155); 1 (1524); 2 (162)
points(SPX_m3M5_minDP10@ped$N1[106],SPX_m3M5_minDP10@ped$N2[106],col="green3",pch=17)


## 4.7 Sex-biased dispersal ####

SPX_mAIc <- sexbias.test(Genotype_SPX[,c(2,5:ncol(Genotype_SPX))], Genotype_SPX[,3], test="mAIc") # t = -1.154827; p = 0.2492973


## 4.8 Distribution of the alleles ####

svg("Allelic_freq_SPX.svg")
par(bg = "#ece8dc")
hist(
  SPX_m3M5_minDP10@p,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(SPX_m3M5_minDP10@p, na.rm = TRUE), max(SPX_m3M5_minDP10@p, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()


## 4.9 Tajima's D ####

tajima_d <- read.table("./SPX_tajima_results_meanDP10_nomissdata.Tajima.D", header = TRUE)
mean_TD <- mean(tajima_d$TajimaD) #-0.5388354
sd_TD <- sd(tajima_d$TajimaD) #0.6801651

hist(
  tajima_d$TajimaD,
  breaks = 100,
  col = "lightblue",
  border = "white",
  main = "Distribution of Tajima'D per SNPs",
  xlab = "Tajima's D",
  ylab = "# SNPs")
abline(v=0, col="red", lty=2, lwd=1.5)

#Significatively less than 0?
wilcox.test(tajima_d$TajimaD, mu = 0, alternative = "less") #V = 54992098, p-value < 2.2e-16

tajima_d$Z_value <- (tajima_d$TajimaD-mean_TD)/sd_TD
hist(tajima_d$Z_value, breaks=100)
abline(v=1.96, col="red", lty=2, lwd=1.5)


## 4.10 Isolation by distance ####

## GEOGRAPHIC DISTANCES
Coord_pop <- read.csv("./Coord_pops.csv", header=TRUE, sep=";")
distgeo_coord <- Coord_pop[, 4:5]
row.names(distgeo_coord) <- Coord_pop$Site
distgeo <- dist(distgeo_coord)


plot(distgeo, as.dist(Pfst_mx_SPXnokin))
#plot(log(distgeo), as.dist(Pfst_mx_CONCAT_nohwe/(1-Pfst_mx_CONCAT_nohwe)))

mantel(distgeo, as.dist(Pfst_mx_SPXnokin), permutations=1000) #r = -0.1827 ; p = 0.83516

svg("./SPX_MantelIBD.svg")
plot(distgeo, as.dist(Pfst_mx_SPXnokin), type="p", pch=19, cex=1.5,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-0.005, 0, 0.005))
axis(side=1, at=c(20000, 40000, 60000, 80000))
title(xlab="Geographical distance (m)", 
      ylab="Genetic distances (FST)", cex.lab=1.5)
dev.off()


################################################################################
# 5. Genetic Analyses - no outliers ####
################################################################################
## 5.1 Dataset ####

SPX_noout <- SPX[!(SPX$ID %in% c("498S1_mergedALL", "792S2_mergedALL")), ]


#Sort by population
SPX_noout_popsorted <- SPX_noout[order(SPX_noout$pop),]
SPX_noout_popsorted_matrix <- as.matrix(SPX_noout_popsorted[,5:ncol(SPX_noout_popsorted)])
SPX_noout_popsorted_matrix <- apply(SPX_noout_popsorted_matrix, 2, as.numeric)


## 5.2 Hierarchical Analysis of Variance ####
### 5.2.1 Using Fs.dosage ####

HierVar_SPX_noout <- fs.dosage(SPX_noout[,c(5:ncol(SPX_noout))], pop=SPX_noout[,2])

plot(HierVar_SPX_noout)

par(mfrow=c(1,1))


### 5.2.2 Plot Pairwise kinship per pop - FsM ####

kinship_mx_noout <- HierVar_SPX_noout$FsM
kinship_mx_noout[upper.tri(kinship_mx_noout)] <- NA
diag(kinship_mx_noout) <- NA

svg("FsM_kinship_SPX_noout.svg", bg = "transparent")
heatmap(kinship_mx_noout, Rowv = NA, Colv = NA, symm = TRUE, 
        col = viridis(100), scale = "none", na.rm = TRUE)
dev.off()
svg("FsM_kinship_SPX_noout_scale.png", bg = "transparent")
image.plot(legend.only = TRUE, zlim = range(HierVar_SPX_noout$FsM, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)
dev.off()


### 5.2.3 Plot Pairwise Fst per pop - Fs2x2 ####

Pfst_mx_noout <- HierVar_SPX_noout$Fst2x2
Pfst_mx_noout[upper.tri(Pfst_mx_noout)] <- NA

heatmap(Pfst_mx_noout, Rowv = NA, Colv = NA, symm = TRUE, 
        col = magma(100), scale = "none", na.rm = TRUE)
image.plot(legend.only = TRUE, zlim = range(Pfst_mx_noout, na.rm = TRUE), 
           col = magma(100), legend.line = 2)


### 5.2.4 Bootstrapping For INBREEDING COEFF & FST PER POP ####

n_bootstraps <- 100
n_loci <- ncol(SPX_noout[, 5:ncol(SPX_noout)])

# Get unique population names
populations <- unique(SPX_noout[, 2])
nb_pop <- length(populations)

# Initialize an empty data frame with specified columns
results_Fs <- data.frame(
  bootstrap = integer(0),
  site = character(0),
  Fis_values = numeric(0),
  #Kinship_values = numeric(0),
  #PairwiseFst_values = numeric(0),
  Fst_values = numeric(0),
  stringsAsFactors = FALSE
)

for (i in 1:n_bootstraps) {
  # Sample loci with replacement
  sampled_loci <- sample(5:ncol(SPX_noout), n_loci, replace = TRUE)
  
  # Subset the data for these loci
  vcf_sample <- SPX_noout[, c(2, sampled_loci)]
  
  # Calculate statistics
  HierVar_boot <- fs.dosage(vcf_sample[, -1], pop = vcf_sample[, 1])
  
  # If fs.dosage returns NULL or an error, skip this iteration
  if (is.null(HierVar_boot)) next
  
  # Extract Fs values for each population
  for (pop_name in names(HierVar_boot$Fs[1,])){
    Fis_values <- HierVar_boot$Fs[1, pop_name] #inbreeding coeff per pop
    Fst_values <- HierVar_boot$Fs[2, pop_name] #Fst coeff per pop 
    
    results_Fs <- rbind(results_Fs, data.frame(
      bootstrap = i,
      site = pop_name,
      Fis_value = Fis_values,
      Fst_value = Fst_values,
      stringsAsFactors = FALSE
    ))
  }
}

str(results_Fs)

boxplot(Fis_value ~ site, data = results_Fs,
        main = "Mean Fis per pop",
        xlab = "", 
        ylab = "Mean Fis",
        col = "grey",
        las=2)
abline(h = 0, col = "#F4A6A6", lwd = 2, lty=2)

boxplot(Fst_value ~ site, data = results_Fs,
        main = "Mean Fst per pop",
        xlab = "", 
        ylab = "Mean Fst",
        col = "grey",
        las=2)
abline(h = 0, col = "#F4A6A6", lwd = 2, lty=2)



## 5.3 Pairwise kinship per individuals ####

kinship_noout <- beta.dosage(SPX_noout_popsorted_matrix)
diag(kinship_noout) = NA
image(kinship_noout)


hist(kinship_noout, breaks = 100, col="black")
svg("Distrib_kinship_noout.svg")
hist(
  kinship_noout,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(kinship_noout, na.rm = TRUE), max(kinship_noout, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()


image(1:nrow(kinship_noout), 1:ncol(kinship_noout), kinship_noout, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(kinship_noout, na.rm = TRUE))
plot.new()
image.plot(legend.only = TRUE, zlim = range(kinship_noout, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)


## WHICH IND KINSHIP > 0.25 (parent-offsprings OR full sib) ??
kinship_noout[upper.tri(kinship_noout)] <- NA
pairs20_noout <- which(kinship_noout > 0.20, arr.ind = TRUE)

# Extract individual names for the pairs
ind1.noout <- SPX_noout_popsorted$ID[pairs20_noout[, 1]]  # First individual in each pair
ind2.noout <- SPX_noout_popsorted$ID[pairs20_noout[, 2]]  # Second individual in each pair
kinship_tablenoout <- data.frame(Individual1 = ind1.noout, Individual2 = ind2.noout, Kinship = kinship_noout[pairs20_noout])
kinship_tablenoout_order <- kinship_tablenoout[order(kinship_tablenoout$Individual1),]





