################################################
##### BAT - Analysis - Population Genetics #####
################################################

setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/5.Thesis/Project/7.RADseq/Project_RAD2024/Bats/Data_Alyzee/overall/Results")

load("./LAST_VERSION.RData") ## Load last R history


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
library(viridisLite)
library(gdsfmt)
library(SNPRelate)
library(LEA)
library(ggridges)
library(adegenet)

## 1.2 Download data ####

BAT_vcf <- read.VCF("../vcf_after_filtering.recode.vcf")
dim(BAT_vcf) #359 35919


################################################################################
# 2. Exploration data ####
################################################################################
## 2.1 Distribution of the genotypes ####
plot(BAT_vcf@ped$N1,BAT_vcf@ped$N2, pch=20)

## 2.2 Heterozygosity per individuals ####
Hetz <- read.table("./min_cov_5_clean_85miss_meanDP_hwe.het", header=TRUE, sep="\t")
Hetz$HETZ <- (Hetz$N_SITES-Hetz$O.HOM.)/Hetz$N_SITES
hist(Hetz$HETZ, breaks=100)

## 2.3 Coverage per individuals ####
Reads <- read.table("../nbreads_INDV_BAT.txt", header=TRUE, sep="\t")
hist(Reads$Reads, breaks=100)
mean(Reads$Reads)#1695454

reads_BAT_order <- Reads[order(Reads$Reads),]
barplot(reads_BAT_order$Reads)

svg("./NbReads.svg")
barplot(reads_BAT_order$Reads, border="white",space=0,
        xlab="", ylab="", axes =F,  col="black")
axis(2, at=seq(0, 3000000, by=1000000), labels=seq(0, 3000000, by=1000000))
title(xlab="Bat individuals", 
      ylab="Number of reads", cex.lab=1.5)
dev.off()


## 2.4 Reads per individuals ####
Cov <- read.table("./min_cov_5_clean_85miss_meanDP_hwe.idepth", header=TRUE, sep="\t")
hist(Cov$MEAN_DEPTH, breaks=100)
mean(Cov$MEAN_DEPTH)#13.88X

## 2.5 Heterozygosity Coverage correlation - minDP5 ####
plot(Hetz$HETZ[Hetz$HETZ>0.18]~Cov$MEAN_DEPTH[Hetz$HETZ>0.18])
mod1 <- lm(Hetz$HETZ[Hetz$HETZ>0.18]~Cov$MEAN_DEPTH[Hetz$HETZ>0.18])
summary(mod1)# R2 = 0.02483  ; p = 0.001598
abline(mod1, col="grey45", lwd = 2) #Control

cor.test(Hetz$HETZ,Cov$MEAN_DEPTH)#cor = 0.1660009; p = 0.001598


## 2.6 Distribution of alleles ####

svg("Allelic_freq_BAT.svg")
hist(
  BAT_vcf@p,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(BAT_vcf@p, na.rm = TRUE), max(BAT_vcf@p, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()


## 2.7 Tajima's D significance #####

tajima_d <- read.table("./Tajima_resultat.Tajima.D", header = TRUE)
mean_TD <- mean(tajima_d$TajimaD, na.rm=TRUE)#-0.08111
sd_TD <- sd(tajima_d$TajimaD, na.rm=TRUE)#0.929733

hist(
  tajima_d$TajimaD,
  breaks = 100,
  col = "black",
  border = "white",
  main = "Distribution of Tajima'D per SNPs",
  xlab = "Tajima's D",
  ylab = "# SNPs")
abline(v=0, col="red", lty=2, lwd=1.5)

#Significatively less than 0?
wilcox.test(tajima_d$TajimaD, mu = 0, alternative = "less") #V = 1014198358, p-value < 2.2e-16


## 2.8 Distribution SNPs per chromosomes ####

vcf_tot <- read.vcfR("../vcf_after_filtering.recode.vcf")
pos<-vcf_tot@fix
pos<-as.data.frame(pos)
chrom<-as.factor(pos$CHROM)
chomo<-levels(chrom)

table(factor(chrom, levels = chomo))

##in MHC??

sum(pos$CHROM=="NC_081845.1" & pos$POS>94750000 & pos$POS < 96250000)
MHC_snp <- pos$ID[pos$CHROM=="NC_081845.1" & pos$POS>94750000 & pos$POS < 96250000]



################################################################################
# 3. Creation of the datasets ####
################################################################################
## 3.1 Including host, sites, sex and host sex ####

popmap <- read.csv("./pop_map_2.txt", header =TRUE, sep="\t")
poplist <- popmap$Site
IDlist <- popmap$ID

## 3.2 SNPs dataset ####

genotype_matrix <- as.matrix(BAT_vcf)
BAT <- data.frame(ID=IDlist, pop=poplist, genotype_matrix)

#Sort by population
BAT_popsorted <- BAT[order(BAT$pop),]
BAT_popsorted_matrix <- as.matrix(BAT_popsorted[,3:ncol(BAT_popsorted)])
BAT_popsorted_matrix <- apply(BAT_popsorted_matrix, 2, as.numeric)


## 3.3 Genotypes dataset ####

Genotype_BAT <- BAT[,-c(1:2)]
Genotype_BAT[Genotype_BAT==0] <- 11
Genotype_BAT[Genotype_BAT==1] <- 12
Genotype_BAT[Genotype_BAT==2] <- 22
Genotype_BAT <- data.frame(ID=BAT[,1],pop=BAT[,2], Genotype_BAT)


## 3.4 Datasets without kin individuals ####
### 3.4.1 Remove kin individuals ####

#Cf Pruning_functions_V2.r

#Pruning & new beta.dosage
kinship_nokin <- pruner(b, 0.12) #for b cf script Pruning_functions_V2.r


### 3.4.1 SNPs dataset without kin individuals ####

BAT_nokin <- BAT_popsorted[BAT_popsorted$ID %in% rownames(kinship_nokin),]
BAT_nokin_matrix <- as.matrix(BAT_nokin[,3:ncol(BAT_nokin)])
BAT_nokin_matrix <- apply(BAT_nokin_matrix, 2, as.numeric)


### 3.4.2 Genotype dataset without kin individuals ####

BAT_nokin_geno <- BAT_nokin[,-c(1:2)]
BAT_nokin_geno[BAT_nokin_geno==0] <- 11
BAT_nokin_geno[BAT_nokin_geno==1] <- 12
BAT_nokin_geno[BAT_nokin_geno==2] <- 22
BAT_nokin_geno <- data.frame(ID=BAT_nokin[,1],pop=BAT_nokin[,2], BAT_nokin_geno)
BAT_nokin_geno_matrix <- apply(as.matrix(BAT_nokin_geno[,3:ncol(BAT_nokin_geno)]), 2, as.numeric)



################################################################################
# 4. Genetic Analysis ####
################################################################################
## 4.1 Genetic diversity ####
### 4.1.1 Heterozygosity ####

DATA = snpgdsVCF2GDS("../vcf_after_filtering.recode.vcf", "vcf_after_filtering_BAT.gds")
DATA = snpgdsOpen("vcf_after_filtering_BAT.gds")
gen_mat <- snpgdsGetGeno('vcf_after_filtering_BAT.gds')
IndHet = rowSums(gen_mat==1, na.rm = TRUE) / rowSums(! is.na(gen_mat))

HetZ_pop <- data.frame(ID=IDlist, pop=poplist, IndHet)
boxplot(HetZ_pop$IndHet~HetZ_pop$pop)

Het_kw <- kruskal.test(HetZ_pop$IndHet~HetZ_pop$pop)#Kruskal-Wallis chi-squared = 19.103, df = 17, p-value = 0.3226

mean_het_per_pop <- aggregate(IndHet ~ pop, data = HetZ_pop, FUN = mean)


### 4.1.2 FUNI ####

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

funi<-get.funiwn(BAT_vcf)
Funi<-funi$Funi

Funi_pop <- data.frame(ID=IDlist, pop=poplist, Funi)
boxplot(Funi_pop$Funi~Funi_pop$pop)

Funi_kw <- kruskal.test(Funi_pop$Funi~Funi_pop$pop)#Kruskal-Wallis chi-squared = 21.784, df = 17, p-value = 0.1932

mean_Funi_pop <- aggregate(Funi ~ pop, data = Funi_pop, FUN = mean)


## 4.2 PCA per site ####
### 4.2.1 On the entire dataset ####

pca <- indpca(Genotype_BAT[,c(1,3:ncol(Genotype_BAT))])
plot(pca)

# % of genetic variance that explain each axis 
pca_eig <- 100*pca$ipca$eig/sum(pca$ipca$eig) #PC1 = 0.600351482% & PC2 = 0.587500395%
barplot(pca_eig[1:10])
sum(pca_eig[1:2])#1.187852
sum(pca_eig[1:10])#5.197666

pca_ALL <- cbind(BAT[,1:2], pca$ipca$li[, 1:2])
pca_sortedpop <- pca_ALL[order(pca_ALL$pop),]

# Identification of the outliers
pca$ipca$li[,1]>10 #1 & 244
Genotype_BAT$ID[244] #315 & 666
pca$ipca$li[,1]<c(-50)#12 & 247
Genotype_BAT$ID[247]#335 & 669
pca$ipca$li[,2]<c(-30)#13 & 236
Genotype_BAT$ID[236]#336 & 657


### 4.2.2 without the recaptures ####

pca_sortedpop$ID <- as.character(pca_sortedpop$ID)
pca_sortedpop2 <- pca_sortedpop[!(pca_sortedpop$ID %in% c("315", "335", "336","666", "669", "657")),]
pop_levels2 <- unique(pca_sortedpop2$pop)
col_pop2 <- viridis(length(pop_levels2))
names(col_pop2) <- pop_levels2  # Associe chaque couleur à une population
pch_pop2 <- 1:length(pop_levels2)
names(pch_pop2) <- pop_levels2

# Appliquer les couleurs en fonction de chaque individu
ind_colors2 <- col_pop2[as.character(pca_sortedpop2$pop)]
indpch2 <- pch_pop2[as.character(pca_sortedpop2$pop)]

svg("PCA_site_BAT_norecapture.svg", bg = "transparent")
plot(pca_sortedpop2[,3], pca_sortedpop2[,4], col=ind_colors2, pch=indpch2, lwd=1,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-8, -6, -4, -2, 0, 2, 4))
axis (side=1, at=c(-2, 0, 2, 4))
title(xlab="PC1 (0.600% explained variation)", 
      ylab="PC2 (0.588% explained variation)", cex.lab=1.5)
legend("bottomright", legend=1:18, col=col_pop2, pch = pch_pop2, title = "Pop")
dev.off()

### 4.2.3 without the kin individuals ####

pca_nokin <- indpca(BAT_nokin_geno[,c(1,3:ncol(BAT_nokin_geno))])
plot(pca_nokin)

# % of genetic variance that explain each axis 
pca_nokin_eig <- 100*pca_nokin$ipca$eig/sum(pca_nokin$ipca$eig) #PC1=0.4858247 ; PC2=0.4325732
barplot(pca_nokin_eig[1:10])
sum(pca_nokin_eig[1:2])#0.9183978
sum(pca_nokin_eig[1:10])#4.236481


pca_nokin_ALL <- cbind(BAT_nokin_geno[,1:2], pca_nokin$ipca$li[, 1:2])
pca_nokin_sortedpop <- pca_nokin_ALL[order(pca_nokin_ALL$pop),]

# Obtenir les couleurs pour chaque population
pop_nokin_levels <- unique(pca_nokin_sortedpop$pop)
col_nokin_pop <- viridis(length(pop_nokin_levels))
names(col_nokin_pop) <- pop_nokin_levels  # Associe chaque couleur à une population
pch_nokin_pop <- 1:length(pop_nokin_levels)
names(pch_nokin_pop) <- pop_nokin_levels

# Appliquer les couleurs en fonction de chaque individu
ind_nokin_colors <- col_nokin_pop[as.character(pca_nokin_sortedpop$pop)]
ind_nokin_pch <- pch_nokin_pop[as.character(pca_nokin_sortedpop$pop)]

svg("PCA_site_BAT_nokin.svg", bg = "transparent")
plot(pca_nokin_sortedpop[,3], pca_nokin_sortedpop[,4], col=ind_nokin_colors, pch=ind_nokin_pch, lwd=1,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-10, 0, 10, 20, 30))
axis (side=1, at=c(-5, 0, 5, 10))
title(xlab="PC1 (0.486% explained variation)", 
      ylab="PC2 (0.433% explained variation)", cex.lab=1.5)
dev.off()
svg("PCA_site_BAT_nokin_legend.svg", bg = "transparent")
plot.new()
legend("topright", legend=1:18, col=col_nokin_pop, pch = pch_nokin_pop, title = "Pop")
dev.off()


### 4.2.4 Parallel Coordinate plot without kin individuals ####

group_colors <- as.factor(BAT_nokin_geno$pop)
color_palette <- viridis(length(unique(group_colors)))  # Assign colors

svg("./Parallele_coord_plot.svg")
matplot(t(as.matrix(pca_nokin$ipca$li[,1:2])), type="l",
        lty=1, col = color_palette[group_colors],
        xlab="", ylab="", axes=F)
axis(side=2, at=c(-10, 0, 10, 20, 30))
axis (side=1, at=c("PC1", "PC2"))
title(xlab="Principal Components", 
      ylab="Normalized value", cex.lab=1.5)
dev.off()
plot.new()
svg("./Parallele_coord_legend.svg")
legend("topright", legend = 1:18, col = color_palette, lty = 1, cex = 0.7)
dev.off()


## 4.3 Rare alleles ####
### 4.3.1 Distribution of the rare alleles per sampling sites ####

BAT_alleles <- BAT[,3:ncol(BAT_popsorted)]
mac <- colSums(BAT_alleles, na.rm = TRUE)
mac_rare <- mac[mac<=10]
loci_Rare <- attributes(mac_rare)$names

Present_in_Npop <- c() #vector that tells number of pop the allele is present in
for (i in 1:length(loci_Rare)){
  locuscol <- match(loci_Rare[i],names(BAT))
  Which_pops <- unique(BAT$pop[BAT[,locuscol]>0])
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


### 4.3.2 Hierarchical analysis on rare alleles ####

BAT_rareall <- cbind(BAT[,c(1,2)],BAT[,names(mac_rare)])

HierVar_BATrareall <- fs.dosage(BAT_rareall[,c(3:ncol(BAT_rareall))], pop=BAT_rareall[,2])


#### 4.3.2.1 Pairwise kinship per pop - FsM - RARE ALLELES ####

kinship_rareall_BAT <- HierVar_BATrareall$FsM
kinship_rareall_BAT[lower.tri(kinship_rareall_BAT)] <- NA

colFsM_rareall <- colorRampPalette(brewer.pal(9, "PuBu"))(100)

svg("./Kinship_RareAll_BAT.svg")
heatmap(kinship_rareall_BAT, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFsM_rareall, scale = "none", na.rm = TRUE, zlim = range(kinship_rareall_BAT, na.rm = TRUE))
dev.off()
svg("./Kinship_RareAll_BAT_scale.svg")
image.plot(legend.only = TRUE, zlim = range(kinship_rareall_BAT, na.rm = TRUE), 
           col = colFsM_rareall, legend.line = 2)
dev.off()


#### 4.3.2.2 Pairwise Fst per pop - Fs2x2 - RARE ALLELES ####

Pfst_RareAll_BAT <- HierVar_BATrareall$Fst2x2
Pfst_RareAll_BAT[upper.tri(Pfst_RareAll_BAT)] <- NA
diag(Pfst_RareAll_BAT) <- NA

colFs2x2_rareall <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

svg("./Fst_RareAll_BAT.svg")
heatmap(Pfst_RareAll_BAT, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFs2x2_rareall, scale = "none", na.rm = TRUE, zlim = range(Pfst_RareAll_BAT, na.rm = TRUE))
dev.off()
svg("./Fst_RareAll_BAT_scale.svg")
image.plot(legend.only = TRUE, zlim = range(Pfst_RareAll_BAT, na.rm = TRUE), 
           col = colFs2x2_rareall, legend.line = 2)
dev.off()


#### 4.3.2.3 Population-specific FIS & FST - RARE ALLELES ####

n_bootstraps <- 100
n_loci_rareall <- ncol(BAT_rareall[, 3:ncol(BAT_rareall)])

# Get unique population names
populations <- unique(BAT_rareall[, 2])
nb_pop <- length(populations)

# Initialize empty lists to store bootstrap results
results_Fs_RareAll <- data.frame(
  bootstrap = integer(0),
  site = character(0),
  Fis_value_rareall = numeric(0),
  Fst_value_rareall = numeric(0),
  stringsAsFactors = FALSE
)

# Store hierarchical variance components
FsMrareall_bootstrap <- list()
Fs2x2rareall_bootstrap <- list()


for (i in 1:n_bootstraps) {
  sampled_loci <- sample(3:ncol(BAT_rareall), n_loci_rareall, replace = TRUE)
  vcf_sample <- BAT_rareall[, c(2, sampled_loci)]
  
  HierVar_boot_rareall <- fs.dosage(vcf_sample[, -1], pop = vcf_sample[, 1])
  
  if (is.null(HierVar_boot_rareall)) next
  
  FsMrareall_bootstrap[[i]] <- HierVar_boot_rareall$FsM
  Fs2x2rareall_bootstrap[[i]] <- HierVar_boot_rareall$Fst2x2
  
  for (pop_name in names(HierVar_boot_rareall$Fs[1, ])) {
    results_Fs_RareAll <- rbind(results_Fs_RareAll, data.frame(
      bootstrap = i,
      site = pop_name,
      Fis_value_rareall = HierVar_boot_rareall$Fs[1, pop_name],
      Fst_value_rareall = HierVar_boot_rareall$Fs[2, pop_name],
      stringsAsFactors = FALSE
    ))
  }
}


# Compute confidence intervals
ci_95 <- function(values) {
  quantile(values, probs = c(0.025, 0.975), na.rm = TRUE)
}

# Compute CI for Fis and Fst values per population
Fis_CIrareall <- aggregate(Fis_value_rareall ~ site, data = results_Fs_RareAll, ci_95)
Fst_CIrareall <- aggregate(Fst_value_rareall ~ site, data = results_Fs_RareAll, ci_95)

# Compute CI for each entry in FsM and Fs2x2 matrices
FsM_matrix_rareall <- simplify2array(FsMrareall_bootstrap)
Fs2x2_matrix_rareall <- simplify2array(Fs2x2rareall_bootstrap)

FsM_CI_matrix_rareall <- apply(FsM_matrix_rareall, c(1,2), ci_95)
Fs2x2_CI_matrix_rareall <- apply(Fs2x2_matrix_rareall, c(1,2), ci_95)


# Boxplot visualization

svg("./BAT_FISbootstrap_RareAll.svg")
boxplot(Fis_value_rareall ~ site, data = results_Fs_RareAll,
        main = "Mean FIS per population in Rare alleles",
        xlab = "", 
        ylab = "Mean FIS",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()

svg("./BAT_FSTbootstrap_RareAll.svg")
boxplot(Fst_value_rareall ~ site, data = results_Fs_RareAll,
        main = "Mean FST per population in Rare alleles",
        xlab = "Population", 
        ylab = "Mean FST",
        col = "lightgrey",
        las = 2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
dev.off()



## 4.4 Hierarchical Analysis of Variance ####

HierVar_BAT <- fs.dosage(BAT[,c(3:ncol(BAT_popsorted))], pop=BAT[,2])

plot(HierVar_BAT)

par(mfrow=c(1,1))


### 4.4.1 Pairwise kinship per pop - FsM ####

kinship_mx <- HierVar_BAT$FsM
kinship_mx[lower.tri(kinship_mx)] <- NA

colFsM <- colorRampPalette(brewer.pal(9, "PuBu"))(100)

svg("./Kinship_Blue_BAT.svg")
heatmap(kinship_mx, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFsM, scale = "none", na.rm = TRUE)
dev.off()
plot.new()
svg("./Kinship_blue_SPX_scale.svg")
image.plot(legend.only = TRUE, zlim = range(HierVar_BAT$FsM, na.rm = TRUE), 
           col = colFsM, legend.line = 2)
dev.off()


### 4.4.2 Pairwise Fst per pop - Fs2x2 ####

Pfst_mx <- HierVar_BAT$Fst2x2
Pfst_mx[upper.tri(Pfst_mx)] <- NA
diag(Pfst_mx) <- NA

colFs2x2 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

svg("./Fst_red_BAT.svg")
heatmap(Pfst_mx, Rowv = NA, Colv = NA, symm = TRUE, 
        col = colFs2x2, scale = "none", na.rm = TRUE)
dev.off()
plot.new()
svg("./Fst_red_BAT_scale.svg")
image.plot(legend.only = TRUE, zlim = range(Pfst_mx, na.rm = TRUE), 
           col = colFs2x2, legend.line = 2)
dev.off()


### 4.4.3 Population-specific FIS & FST ####
#### 4.4.3.1 Linkage disequilibrium between loci ####

# Calculate in the cluster


#### 4.4.3.2 Bootstrapping Pop-specific FIS & FST ####

n_bootstraps <- 100
n_loci <- ncol(BAT_popsorted[, 3:ncol(BAT_popsorted)])

# Get unique population names
populations <- unique(BAT_popsorted[, 2])
nb_pop <- length(populations)

# Initialize an empty data frame with specified columns
results_Fs_BAT <- data.frame(
  bootstrap = integer(0),
  site = character(0),
  Fis_values = numeric(0),
  #Kinship_values = numeric(0),
  #PairwiseFst_values = numeric(0),
  Fst_values = numeric(0),
  stringsAsFactors = FALSE
)

# Store hierarchical variance components
FsMBAT_bootstrap <- list()
Fs2x2BAT_bootstrap <- list()


for (i in 1:n_bootstraps) {
  # Sample loci with replacement
  sampled_loci <- sample(3:ncol(BAT_popsorted), n_loci, replace = TRUE)
  
  # Subset the data for these loci
  vcf_sample <- BAT_popsorted[, c(2, sampled_loci)]
  
  # Calculate statistics
  HierVar_boot_BAT <- fs.dosage(vcf_sample[, -1], pop = vcf_sample[, 1])
  
  # If fs.dosage returns NULL or an error, skip this iteration
  if (is.null(HierVar_boot_BAT)) next
  
  FsMBAT_bootstrap[[i]] <- HierVar_boot_BAT$FsM
  Fs2x2BAT_bootstrap[[i]] <- HierVar_boot_BAT$Fst2x2
  
  # Extract Fs values for each population
  for (pop_name in names(HierVar_boot_BAT$Fs[1,])){
    Fis_values_BAT <- HierVar_boot_BAT$Fs[1, pop_name] #inbreeding coeff per pop
    Fst_values_BAT <- HierVar_boot_BAT$Fs[2, pop_name] #Fst coeff per pop 
    
    results_Fs_BAT <- rbind(results_Fs_BAT, data.frame(
      bootstrap = i,
      site = pop_name,
      Fis_value_BAT = Fis_values_BAT,
      Fst_value_BAT = Fst_values_BAT,
      stringsAsFactors = FALSE
    ))
  }
}

# Compute confidence intervals
ci_95 <- function(values) {
  quantile(values, probs = c(0.025, 0.975), na.rm = TRUE)
}

# Compute CI for Fis and Fst values per population
Fis_CIBAT <- aggregate(Fis_value_BAT ~ site, data = results_Fs_BAT, ci_95)
Fst_CIBAT <- aggregate(Fst_value_BAT ~ site, data = results_Fs_BAT, ci_95)

# Compute CI for each entry in FsM and Fs2x2 matrices
FsM_matrixBAT <- simplify2array(FsMBAT_bootstrap)
Fs2x2_matrixBAT <- simplify2array(Fs2x2BAT_bootstrap)

FsM_CI_matrixBAT <- apply(FsM_matrixBAT, c(1,2), ci_95)
Fs2x2_CI_matrixBAT <- apply(Fs2x2_matrixBAT, c(1,2), ci_95)

str(results_Fs_BAT)

svg("./BAT_popspecific_FIS.svg")
boxplot(Fis_value_BAT ~ site, data = results_Fs_BAT,
        main = "Mean Fis per pop",
        xlab = "", 
        ylab = "Mean Fis",
        col = "lightgrey",
        las=2)
abline(h = 0, col = "red", lwd = 2, lty=2)
dev.off()

svg("./BAT_popspecific_FST.svg")
boxplot(Fst_value_BAT ~ site, data = results_Fs_BAT,
        main = "Mean Fst per pop",
        xlab = "", 
        ylab = "Mean Fst",
        col = "lightgrey",
        las=2)
abline(h = 0, col = "red", lwd = 2, lty=2)
dev.off()


# Mean overall stat
mean(results_Fs_BAT$Fis_value_BAT[results_Fs_BAT$site=="All"]) #0.003568107
mean(results_Fs_BAT$Fst_value_BAT[results_Fs_BAT$site=="All"]) #0.003316113


## 4.5 Pairwise kinship across individuals ####
### 4.5.1 Kinship matrix ####

kinship <- beta.dosage(BAT_popsorted_matrix)
diag(kinship) = NA
image(kinship)


svg("Distribution_kinship_BAT.svg")
hist(
  kinship,
  breaks = 100,
  col = "black",
  border = "white",
  main="",
  xlab = "Values",
  ylab = "Frequency",
  xlim = c(min(kinship, na.rm = TRUE), max(kinship, na.rm = TRUE)),
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1,
  las = 1
)
dev.off()




popcount <- table(BAT$pop)

svg("./kinsihp_overall_indv.svg")
image(1:nrow(kinship), 1:ncol(kinship), kinship, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(kinship, na.rm = TRUE))
dev.off()
plot.new()
svg("./kinsihp_overall_indv_scale.svg")
image.plot(legend.only = TRUE, zlim = range(kinship, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)
dev.off()

## LINKS > 0.2
kinship[upper.tri(kinship)] <- NA
pairs20_BAT <- which(kinship > 0.2, arr.ind = TRUE) #27

# Extract individual names for the pairs
ind1 <- BAT_popsorted$ID[pairs20_BAT[, 1]]  # First individual in each pair
ind2 <- BAT_popsorted$ID[pairs20_BAT[, 2]]  # Second individual in each pair
kinship_table <- data.frame(Individual1 = ind1, Individual2 = ind2, Kinship = kinship[pairs20_BAT])
kinship_table_order <- kinship_table[order(kinship_table$Individual1),]

# Count pairs per individuals
pair_counts <- table(kinship_table$Individual2)
pair_counts_df <- data.frame(Individual = names(pair_counts), Num_Pairs = as.numeric(pair_counts))


### 4.5.2 Fraternity coefficient ####

snpgdsVCF2GDS('../vcf_after_filtering.recode.vcf','vcf_after_filtering.recode.gds',method=c("biallelic.only"))
gds<- snpgdsOpen("vcf_after_filtering.recode.gds",allow.duplicate=TRUE)
k <- snpgdsIBDMoM(gds,autosome.only=FALSE)
plot(k$k1~k$k0,pch=20)

colors <- ifelse((k$k0 >= 0 & k$k0 <= 0.18 & k$k1 >= 0.8 & k$k1 <= 1), "#7570b3",
                 ifelse((k$k0 >= 0.15 & k$k0 <= 0.35 & k$k1 >= 0.35 & k$k1 <= 0.65), "#d95f02",
                        ifelse((k$k0 >= 0 & k$k0 <= 0.1 & k$k1 >= 0 & k$k1 <= 0.1),"#F0E442", "grey" )))

svg("./testCoefFrat.svg")
pdf("./testCoefFrat.pdf")
tiff("./FratCoeff2.tiff")
plot(k$k0, k$k1, col = colors, pch = 20,cex=1.5, xlab = "k0", ylab = "k1")
dev.off()

### Parent-offspring relationship
data_points_bis <- k$k0 >= 0 & k$k0 <= 0.18 & k$k1 >= 0.8 & k$k1 <= 1
colnames(data_points_bis)<-k$sample.id
rownames(data_points_bis)<-k$sample.id

indices <- which(data_points_bis == TRUE, arr.ind = TRUE)
results <- data.frame(
  Individu1 = rownames(data_points_bis)[indices[, 1]],  
  Individu2 = colnames(data_points_bis)[indices[, 2]],
  data_points_bis = data_points_bis[indices],
  k0_values = k$k0[indices],
  k1_values = k$k1[indices]
)
write.csv(results,"PO_relationship.csv")


### Full Siblings relationship
data_points_FS <- k$k0 >= 0.15 & k$k0 <= 0.35 & k$k1 >= 0.35 & k$k1 <= 0.65
colnames(data_points_FS)<-k$sample.id
rownames(data_points_FS)<-k$sample.id

indices_FS <- which(data_points_FS == TRUE, arr.ind = TRUE)
results_FS <- data.frame(
  Individu1 = rownames(data_points_FS)[indices_FS[, 1]],  
  Individu2 = colnames(data_points_FS)[indices_FS[, 2]],
  data_points_FS = data_points_FS[indices_FS],  
  k0_values = k$k0[indices_FS],
  k1_values = k$k1[indices_FS]
)
write.csv(results_FS,"FS_relationship.csv")


### Identical
data_points_self <- k$k0 >= 0 & k$k0 <= 0.1 & k$k1 >= 0 & k$k1 <= 0.1
diag(data_points_self)=NA
colnames(data_points_self)<-k$sample.id
rownames(data_points_self)<-k$sample.id

indices_self <- which(data_points_self == TRUE, arr.ind = TRUE)
results_self <- data.frame(
  Individu1 = rownames(data_points_self)[indices_self[, 1]],  
  Individu2 = colnames(data_points_self)[indices_self[, 2]],
  data_points_self = data_points_self[indices_self]
)
write.csv(results_self,"self_relationship.csv")


### 4.5.3 Genotypes - kinship special cases ####

which(BAT_vcf@ped$id=="434") #62
which(BAT_vcf@ped$id=="497") #116
which(BAT_vcf@ped$id=="810") #353

svg("./Genotypes_kinshipspecialcases.svg")
plot(BAT_vcf@ped$N1,BAT_vcf@ped$N2, pch=20)
points(BAT_vcf@ped$N1[62],BAT_vcf@ped$N2[62],col="#71C837",pch=18)
points(BAT_vcf@ped$N1[116],BAT_vcf@ped$N2[116],col="#FF6600",pch=18)
points(BAT_vcf@ped$N1[353],BAT_vcf@ped$N2[353],col="#BC5FD3",pch=18)
dev.off()


### 4.5.4 Pairwise kinship without kin individuals ####

diag(kinship_nokin) = NA
image(kinship_nokin)

popcount_nikin <- table(BAT_popsorted_nokin$pop)

svg("./kinship_overall_indv_NOKIN.svg")
image(1:nrow(kinship_nokin), 1:ncol(kinship_nokin), kinship_nokin, 
      col = viridis(100), axes = FALSE, xlab = "", ylab = "", 
      zlim = range(kinship_nokin, na.rm = TRUE))
dev.off()
plot.new()
svg("./kinship_overall_indv_nokin_scale.svg")
image.plot(legend.only = TRUE, zlim = range(kinship_nokin, na.rm = TRUE), 
           col = viridis(100), legend.line = 2)
dev.off()


## 4.6 Sex-biased dispersal ####

Data_ALL_2023 <- read.table("./DATA_ALL_2023.txt", header=TRUE, sep="\t", dec=",")

# On adults only!
Genotype_BAT$sex <- NA
Genotype_BAT$age <- NA
common_ids <- Genotype_BAT$ID %in% Data_ALL_2023$ID  # identifie les ID communs
Genotype_BAT$sex[common_ids] <- Data_ALL_2023$Sex[match(Genotype_BAT$ID[common_ids], Data_ALL_2023$ID)]
Genotype_BAT$age[common_ids] <- Data_ALL_2023$Age[match(Genotype_BAT$ID[common_ids], Data_ALL_2023$ID)]
Genotype_BAT_order <- Genotype_BAT[,c(1:2, 35922, 35923, 3:35921)] #reorder columns
Genotype_BAT_ADonly <- Genotype_BAT_order[Genotype_BAT_order$age== "AD" | Genotype_BAT_order$age== "SAD",]


BAT_mAIc <- sexbias.test(Genotype_BAT_ADonly[,c(2,5:ncol(Genotype_BAT_ADonly))], 
                         Genotype_BAT_ADonly[,3], 
                         test="mAIc") # t = -1.788627; p = 0.07487271



## 4.7 Isolation by distance ####

Coord_pop <- read.csv("./Coord_pops_RADBat.csv", header=TRUE, sep=",")
distgeo_coord <- Coord_pop[, 6:7]
row.names(distgeo_coord) <- Coord_pop$Site
distgeo <- dist(distgeo_coord)

svg("./Isolation_distance.svg")
plot(distgeo, as.dist(Pfst_mx))
dev.off()
#plot(log(distgeo), as.dist(Pfst_mx_CONCAT_nohwe/(1-Pfst_mx_CONCAT_nohwe)))

mantel(distgeo, as.dist(Pfst_mx), permutations=1000) #r = 0.2356 ; p = 0.083916



## 4.8 Admixture plot ####

#between no kin individuals
BAT_popsorted_nokin_matrix <- BAT_popsorted[BAT_popsorted$ID %in% colnames(kinship_nokin),]

# Extract the genetic data (columns 3:ncol(BAT_popsorted_nokin))
geno <- t(as.matrix(BAT_popsorted_nokin_matrix[, 3:ncol(BAT_popsorted_nokin)]))
dim(geno)#35919 328
geno[is.na(geno)]=9

write.table(geno, file = "geno.geno", quote = F, col.names = F, row.names = F, sep = "")

#Perform admixture
structure <- snmf("geno.geno", K=1:10, repetitions = 10,entropy=TRUE,project="new")

#Cross entropy for each K
entropy <- sapply(1:10, function(k) cross.entropy(structure, K=k))
best_indices <- which(entropy == min(entropy), arr.ind = TRUE) #K1, run 8

#Boxplot entropy
svg("Entropy_boxplot_10runs.svg", bg = "transparent")
boxplot(entropy, col="grey", axes = F)
axis(side=2, at=c(0.495, 0.500, 0.505, 0.510, 0.515))
axis (side=1, at=c(1:10))
title(xlab="Number of clusters (K)", 
      ylab="Entropy", cex.lab=1.5)
dev.off()


#Table with assignment values for each cluster (each K)
Qtable = read.table("./geno.snmf/K1/run8/geno_r8.1.Q")
pop_labels=

#Plot cluster
AdmixPlot <- barplot(as.matrix(t(Qtable)), col = "lightblue", 
        border = NA, space = 0, names.arg=)


## 4.9 Estimate of Ne ####
### 4.9.1 NeEstimator ####

#ON THE CLUSTER


### 4.9.2 By hand - JG script ####

library(gaston)
library(hierfstat)
bed <- read.VCF("../vcf_after_filtering.recode.vcf",get.info=TRUE)

#as vcf is missing chr id, the following identify them
chrom.lim <- which(bed@snps$pos[2:35919]-bed@snps$pos[1:35918]< -1000000)
lg<-vector(length=length(bed@p))
lg[1:chrom.lim[1]]<-1
for (i in 2:length(chrom.lim)){
  lg[(chrom.lim[i-1]+1):chrom.lim[i]]<-i
}

bed@snps$chr<-lg

#5 remaining loci unassigned
bed@snps$chr[bed@snps$chr==0] <- NA

#nb inds
n<-dim(bed)[2]

#sampling correction. for a start, use all inds
r2s<-1/n+3.19/n^2


#maf filter at 40%
bedmaf40<-bed[,bed@p>=0.4]

r2<-(cor(as.matrix(bedmaf40),use="pairwise.complete.obs"))^2

mr2<-mean(mat2vec(r2))

#identify pairs of loci on different chr
sc<-outer(bedmaf40@snps$chr,bedmaf40@snps$chr,FUN<-function(x,y) x!=y)

#mean of r2 for pairs of loci between chromosomes 
mr2<-mean(r2*sc)


#r2 prime, as a first pass
#TODO: get n for each pair of loci instead of overall
r2p<-mr2-r2s

#formula for NE in the RM case

ne<-(1/3+(1/9-2.76*r2p)^.5)/2/r2p
ne

#and the NE estimate for maf 0.4:

#[1] 1543.99

#same, but for maf10
bedmaf10<-bed[,bed@p>=0.1]

#get all r2
r2<-(cor(as.matrix(bedmaf10),use="pairwise.complete.obs"))^2

#identify pairs between chromosomes
sc<-outer(bedmaf10@snps$chr,bedmaf10@snps$chr,FUN<-function(x,y) x!=y)
mr2<-mean(r2*sc,na.rm=TRUE)


r2p<-mr2-r2s
ne<-(1/3+(1/9-2.76*r2p)^.5)/2/r2p
ne

### 4.9.3 With mutation rate ####

# Theta Watterson

# #From cluster:
# module load gcc/12.3.0
# module load vcftools/0.1.16
# vcftools --vcf vcf_after_filtering.recode.vcf --site-pi --out pi_output
# grep -v "^CHROM" pi_output.sites.pi | wc -l ##35919
# awk 'BEGIN { an = 0; for (i = 1; i < 718; i++) an += 1/i; print an }' ##7.15299
# theta=$(echo "scale=8; 35919 / 7.15299" | bc) ##theta=5021.53644839

theta = 5021.53644839
u=10^-8

Ne_t=theta/(4*u) #125'538'411'210

K=35919
n=c(1:718)
an=sum(1/n-1)

theta2 = K / an #-50.52996
Ne_t2=theta2/(4*u) #-1263248976


# Heterozygozity

H=mean(HetZ_pop$IndHet)
u=10^-8

Ne_co = H/(4*u*(1-H)) #5'746'398
