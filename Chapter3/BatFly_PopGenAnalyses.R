####################################################
#### Bat Flies - Analysis - Population Genetics ####
####################################################


setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/5.Thesis/Project/7.RADseq/Project_RAD2024/Parasites/Nyct/V2_Results/M3m4/withHWE/mac/no_outliers")

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


## 1.2 Download data ####

Nyct_M3m4 <- read.VCF("./min_cov_5_clean_80miss_cov10_65_hweq_mac.recode.vcf")
dim(Nyct_M3m4) #426 1544


################################################################################
# 2. Exploring outliers ####
################################################################################
## 2.1 Kinship matrix ####

k <- beta.dosage(Nyct_M3m4)

## 2.2 Explore Inbreeding and kinship ####
links_count <- apply(k, 1, function(x) sum(x > 0.125, na.rm = TRUE))
plot(diag(k), log(links_count), 
     xlab = "Inbreeding Coefficient", 
     ylab = "Number of Links > 0.125", 
     main = "Inbreeding vs. Relatedness Links", pch=20, axes=F)
axis(1)  # Default x-axis
custom_ticks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
axis(2, at = log(custom_ticks), labels = custom_ticks)  # Custom Y-axis labels
points(diag(k)[links_count == 363], log(links_count[links_count == 363]),
       col = "orange", pch = 19)
points(diag(k)[links_count == 421], log(links_count[links_count == 421]),
       col = "red", pch = 19)


## 2.3 Pruning for k>0.125 ####

links_count2 <- apply(k, 1, function(x) sum(x > 0.0625, na.rm = TRUE))

Prune0125 <- read.table("./beta_prune0125.table", header=TRUE)
Prune0125_matrix <- as.matrix(Prune0125)
Prune0125_matrix <- apply(Prune0125, 2, as.numeric)
image(1:nrow(Prune0125_matrix), 1:ncol(Prune0125_matrix), Prune0125_matrix, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(Prune0125_matrix, na.rm = TRUE))
plot.new()
image.plot(legend.only = TRUE, zlim = range(Prune0125_matrix, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)


## 2.4 Distribution of genotypes - plot ####

point_colors <- ifelse(Nyct_M3m4@ped$id %in% rownames(Prune0125), "black", "blue")
plot(Nyct_M3m4@ped$N1,Nyct_M3m4@ped$N2, col=point_colors, pch=20)
points(Nyct_M3m4@ped$N1[245],Nyct_M3m4@ped$N2[245],col="red",pch=20)
points(Nyct_M3m4@ped$N1[167],Nyct_M3m4@ped$N2[167],col="orange",pch=20)


Prune00625 <- read.table("./beta_prune00625.table", header=TRUE)
point_colors <- ifelse(Nyct_M3m4@ped$id %in% rownames(Prune00625), "black", "purple")
plot(Nyct_M3m4@ped$N1,Nyct_M3m4@ped$N2, col=point_colors, pch=20)
points(Nyct_M3m4@ped$N1[245],Nyct_M3m4@ped$N2[245],col="red",pch=20)
points(Nyct_M3m4@ped$N1[167],Nyct_M3m4@ped$N2[167],col="orange",pch=20)


################################################################################
# 3. Creating datasets ####
################################################################################
## 3.1 Including host, sites, sex and host sex ####

popmap <- read.csv("./Site_sexs_M3m4_noout.csv", header =TRUE, sep=";")
poplist <- popmap$Site
nyctsex <- popmap$Sex_Nyct
batsex <- popmap$Sex_Bat
IDlist <- popmap$ID

## 3.2 Add relative heterozygosity ####

# Per individuals
HetZ_indiv <- read.table("./min_cov_5_clean_80miss_cov10_65_hweq_mac.het", header=TRUE, sep="\t")
HetZ_indiv$HETZ <- (HetZ_indiv$N_SITES-HetZ_indiv$O.HOM.)/HetZ_indiv$N_SITES

#Per population
HetZ_pop <- data.frame(ID=IDlist, pop=poplist, nyctsex=nyctsex, batsex=batsex, HetZ_indiv)
write.table(HetZ_pop, "./HetZ_pop_noout.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## 3.3 SNPs dataset ####

genotype_matrix <- as.matrix(Nyct_M3m4)
Nyct <- data.frame(ID=IDlist, pop=poplist, nyctsex=nyctsex, batsex=batsex, genotype_matrix)
Nyct_noout <- Nyct[!(Nyct$ID %in% c("508N1.merged", "553N1.merged")), ]


#Sort by population
Nyct_noout_popsorted <- Nyct_noout[order(Nyct_noout$pop),]
Nyct_noout_popsorted_matrix <- as.matrix(Nyct_noout_popsorted[,5:ncol(Nyct_noout_popsorted)])
Nyct_noout_popsorted_matrix <- apply(Nyct_noout_popsorted_matrix, 2, as.numeric)


## 3.4 Genotype dataset ####

Genotype_Nyct_noout <- Nyct_noout[,-c(1:4)]
Genotype_Nyct_noout[Genotype_Nyct_noout==0] <- 11
Genotype_Nyct_noout[Genotype_Nyct_noout==1] <- 12
Genotype_Nyct_noout[Genotype_Nyct_noout==2] <- 22
Genotype_Nyct_noout <- data.frame(ID=Nyct_noout[,1],pop=Nyct_noout[,2],nyctsex=Nyct_noout[,3], batsex=Nyct_noout[,4], Genotype_Nyct_noout)

Genotype_Nyct_noout_matrix <- apply(as.matrix(Genotype_Nyct_noout[,5:1548]), 2, as.numeric)
heatmap(Genotype_Nyct_noout_matrix)


################################################################################
# 4. Genetic Analyses ####
################################################################################
## 4.1 Heterozygosity per sampling sites and individuals ####

boxplot(HetZ_pop$HETZ~HetZ_pop$pop)

hist(HetZ_pop$HETZ, breaks=100)
abline(v=0.06229773, col="red", lwd=2)#min (696N1 from Yvonand)
abline(v=0.06303116, col="purple", lwd=2)#508N1
abline(v=0.06945482, col="blue", lwd=2)#553N1


## 4.2 Distribution of the rare alleles per sampling sites ####

Nyct_noout_alleles <- Nyct_noout[,5:1548]
mac <- colSums(Nyct_noout_alleles, na.rm = TRUE)
mac_rare <- mac[mac==10]
loci_Rare <- attributes(mac_rare)$names

Present_in_Npop <- c() #vector that tells number of pop the allele is present in
for (i in 1:length(loci_Rare)){
  locuscol <- match(loci_Rare[i],names(Nyct_noout))
  Which_pops <- unique(Nyct_noout$pop[Nyct_noout[,locuscol]>0])
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


## 4.3 Allelic frequency distribution ####

svg("Allelic_freq.svg")
par(bg = "#ece8dc")
hist(
  Nyct_M3m4@p,
  breaks = 100,
  col = "#69b3a2",
  border = "#ece8dc",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(Nyct_M3m4@p, na.rm = TRUE), max(Nyct_M3m4@p, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()


## 4.4 Tajima's D ####
#TESTING SIGNIFICANCE TAJIMA'S D

tajima_d <- read.table("../nomac/hwe_nomac_tajima_results.Tajima.D", header = TRUE)
mean_TD <- mean(tajima_d$TajimaD) #-0.7084649
sd_TD <- sd(tajima_d$TajimaD) #0.3554116

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
wilcox.test(tajima_d$TajimaD, mu = 0, alternative = "less") #V = 1004951, p-value < 2.2e-16

tajima_d$Z_value <- (tajima_d$TajimaD-mean_TD)/sd_TD
hist(tajima_d$Z_value, breaks=100)
abline(v=1.96, col="red", lty=2, lwd=1.5)


## 4.5 Hierarchical Analysis of Variance ####
### 4.5.1 Fs.dosage matrices ####

HierVar_M3m4_noout <- fs.dosage(Nyct_noout[,c(5:1548)], pop=Nyct_noout[,2])

plot(HierVar_M3m4_noout)

par(mfrow=c(1,1))


## Fs.dosage on NYCT to have the same ranges

global_zlim_FsM <- range(HierVar_SPX_noout$FsM, HierVar_M3m4_noout$FsM, na.rm = TRUE)
global_zlim_Fst2x2 <- range(HierVar_SPX_noout$Fst2x2, HierVar_M3m4_noout$Fst2x2, na.rm = TRUE)


### 4.5.2 Pairwise kinship per pop - FsM ####

kinship_mx_noout <- HierVar_M3m4_noout$FsM
kinship_mx_noout[lower.tri(kinship_mx_noout)] <- NA

colFsM_noout <- colorRampPalette(brewer.pal(9, "PuBu"))(100)

svg("FsM_kinship_noout_NYCT.svg", bg = "transparent")
heatmap(kinship_mx_noout, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFsM_noout, scale = "none", na.rm = TRUE, zlim=global_zlim_FsM)
dev.off()
svg("FsM_kinship_noout_scale_NYCT.png", bg = "transparent")
image.plot(legend.only = TRUE, zlim = range(HierVar_M3m4_noout$FsM, na.rm = TRUE), 
           col = colFsM_noout, legend.line = 2, zlim=global_zlim_FsM)
dev.off()


### 4.5.3 Plot Pairwise Fst per pop - Fs2x2 ####

Pfst_mx_noout <- HierVar_M3m4_noout$Fst2x2
Pfst_mx_noout[upper.tri(Pfst_mx_noout)] <- NA
diag(Pfst_mx) <- NA

colFs2x2 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

svg("./Fst2x2_noout_NYCT.svg")
heatmap(Pfst_mx_noout, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFs2x2, scale = "none", na.rm = TRUE, zlim=global_zlim_Fst2x2)
dev.off()
svg("./Fst2x2_noout_NYCT_scale.svg")
image.plot(legend.only = TRUE, zlim = range(Pfst_mx_noout, na.rm = TRUE), 
           col = colFs2x2, legend.line = 2, zlim=global_zlim_Fst2x2)
dev.off()


### 4.5.4 Bootstrapping Pop-specific FIS & FST ####


n_bootstraps <- 100
n_loci <- ncol(Nyct_noout[, 5:1548])

# Get unique population names
populations <- unique(Nyct_noout[, 2])
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
  sampled_loci <- sample(5:1548, n_loci, replace = TRUE)
  
  # Subset the data for these loci
  vcf_sample <- Nyct_noout[, c(2, sampled_loci)]
  
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


## 4.6 Pairwise kinship per individuals ####
### 4.6.1 Kinship matrix ####

kinship_noout <- beta.dosage(Nyct_noout_popsorted_matrix)
diag(kinship_noout) = NA


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



image(kinship_noout)


image(1:nrow(kinship_noout), 1:ncol(kinship_noout), kinship_noout, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(kinship_noout, na.rm = TRUE))
plot.new()
image.plot(legend.only = TRUE, zlim = range(kinship_noout, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)


### 4.6.2 Exploring kinship among individuals ####

## WHICH IND KINSHIP > 0.25 (parent-offsprings OR full sib) ??
kinship_noout[upper.tri(kinship_noout)] <- NA
pairs25_noout <- which(kinship_noout > 0.25, arr.ind = TRUE)


## WHICH IND KINSHIP > 0 ??
pairs0_noout <- which(kinship_noout > 0, arr.ind = TRUE)

# Extract individual names for the pairs
ind1.1 <- Nyct_noout_popsorted$ID[pairs0_noout[, 1]]  # First individual in each pair
ind2.1 <- Nyct_noout_popsorted$ID[pairs0_noout[, 2]]  # Second individual in each pair
kinship_table2 <- data.frame(Individual1 = ind1.1, Individual2 = ind2.1, Kinship = kinship_noout[pairs0_noout])

## WHICH IND KINSHIP > 0.125 (Half-sib) ??
pairs0125_noout <- which(kinship_noout > 0.125 & kinship_noout <0.25, arr.ind = TRUE) #980 pairs

# Extract individual names for the pairs
ind1.halfSib <- Nyct_noout_popsorted$ID[pairs0125_noout[, 1]]  # First individual in each pair
ind2.halfSib <- Nyct_noout_popsorted$ID[pairs0125_noout[, 2]]  # Second individual in each pair
kinship_table_halfSib <- data.frame(Individual1 = ind1.halfSib, Individual2 = ind2.halfSib, Kinship = kinship_noout[pairs0125_noout])

hist(kinship_table2$Kinship, breaks=100)
abline(v=0.125, col="red", lwd=2)#halfsib


## 4.7 PCA per site ####

pca <- indpca(Genotype_Nyct_noout[,c(1,5:ncol(Genotype_Nyct_noout))])
plot(pca)

# % of genetic variance that explain each axis 
pca_eig <- 100*pca$ipca$eig/sum(pca$ipca$eig) #PC1 = 0.95186044% & PC2 = 0.93598838%
barplot(pca_eig[1:10])
sum(pca_eig[1:2])#1.887849
sum(pca_eig[1:10])#8.666291


pca_ALL <- cbind(Nyct_noout[,1:4], pca$ipca$li[, 1:2])
pca_sortedpop <- pca_ALL[order(pca_ALL$pop),]

# Obtenir les couleurs pour chaque population
pop_levels <- unique(pca_sortedpop$pop)
col_pop <- viridis(length(pop_levels))
names(col_pop) <- pop_levels  # Associe chaque couleur à une population


# Appliquer les couleurs en fonction de chaque individu
ind_colors <- col_pop[as.character(pca_sortedpop$pop)]

svg("PCA_site.svg", bg = "transparent")
plot(pca_sortedpop[,5], pca_sortedpop[,6], col=ind_colors, pch=19, lwd=4,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))
axis (side=1, at=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))
title(xlab="PC1 (0.949% explained variation)", 
      ylab="PC2 (0.934% explained variation)", cex.lab=1.5)
dev.off()
plot.new()
svg("PCA_site_legend.svg", bg = "transparent")
legend("topright", legend=1:14, col=col_pop, pch = 19, title = "Pop")
dev.off()


## 4.8 Isolation by distance ####

str(popmap)
str(Pfst_mx)


## GEOGRAPHIC DISTANCES
Coord_pop <- read.csv("./Coord_pops.csv", header=TRUE, sep=";")
distgeo_coord <- Coord_pop[, 4:5]
row.names(distgeo_coord) <- Coord_pop$Site
distgeo <- dist(distgeo_coord)


plot(distgeo, as.dist(Pfst_mx_noout))
#plot(log(distgeo), as.dist(Pfst_mx_CONCAT_nohwe/(1-Pfst_mx_CONCAT_nohwe)))

mantel(distgeo, as.dist(Pfst_mx_noout), permutations=1000) # -r = 0.0002689 ; p = 0.5025
model <- lm(as.vector(as.dist(Pfst_mx_noout))~as.vector(distgeo))
summary(model)

plot(distgeo, as.dist(Pfst_mx_noout), type="p", pch=19, cex=1.5,
     xlab="", ylab="", axes =F, col="#4CB8B1")
axis(side=2, at=c(-0.005, -0.004, -0.0020, 0,0.0020, 0.0040))
axis(side=1, at=c(0, 20000, 40000, 60000, 80000, 100000))
title(xlab="Geographical distance (m)", 
      ylab="Genetic distances (FST)", cex.lab=1.5)
abline(model, col="red", lty=2)


## 4.9 Sex-biased dispersal ####

### Calculate corrected Assignments Index

AIc_ALL_mAIc <- sexbias.test(Genotype_Nyct_noout[,c(2,5:ncol(Genotype_Nyct_noout))], Genotype_Nyct_noout[,3], test="mAIc") # t = -1.06869; p = 0.2862803
AIc_ALL_vAIc <- sexbias.test(Genotype_Nyct_noout[,c(2,5:ncol(Genotype_Nyct_noout))], Genotype_Nyct_noout[,3], test="vAIc") # t = 0.8477027; p = 0.256687

