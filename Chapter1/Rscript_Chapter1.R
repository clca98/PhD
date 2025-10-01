#########################################################
##### Genetic Comparison of Bed Bugs at local scale #####
#########################################################

setwd("C:/Users/ccastex/OneDrive - Université de Lausanne/6.Cimex paper")


#########################################################
#### 1. Loading packages ####
#########################################################

library("hierfstat")
library ("pegas")
library("adegenet")
library("ade4")
library("factoextra")
library("ggfortify")
library("ggplot2")
library("dplyr")
library("pals")
library("vegan")
library("corrplot")
library("ade4")
library("pheatmap")
library("RColorBrewer")
library("hrbrthemes")
library("corMLPE")
library("ape")

# remove.packages("rlang")       
# install.packages("rlang")


setwd("./6.Microsat/Res_1indallpop/")



#########################################################
#### 2. Analysis with all individuals per sites #####
#########################################################
## 2.1 Download table ####
microsat_all<-read.fstat.data("./6.Microsat/Res_1indallpop/MicrosatCimexALL_10loci_poporder.dat")

## 2.2 Calculate kinship per individuals ####

#Dosage matrix
dos_microsatall <- fstat2dos(microsat_all[,2:ncol(microsat_all)], diploid=TRUE)

#Estimate pairwise kinships
r_matrixall <- beta.dosage(dos_microsatall)
#values=pairwise kinships
#diag=inbreeding
colnames(r_matrixall) = paste(1:370, sep = "")
rownames(r_matrixall) = paste(1:370, sep = "")
kinship_part <- r_matrixall
diag(kinship_part) <- 0
lower_triangular <- lower.tri(r_matrixall,diag = F)

# Set the diagonal to zeros (keep only kinship values)
diag(lower_triangular) <- 0

# Multiply r_matrixall by lower_triangular to keep only the lower triangular part
kinship_part <- r_matrixall * lower_triangular
kinship_part[kinship_part == 0] <- NA

#Estimate kinship for bat
mean(kinship_part[1:59, 1:59], na.rm=TRUE) #0.3748651
mean(kinship_part[60:370, 60:370], na.rm=TRUE) #0.09504722
mean(kinship_part, na.rm=TRUE) #5.737926e-17


#Estimate kinship for human

jpeg("pheatmap_relatedness_allind_scale.jpg")
svg("./pheatmap_relatedness_allind_scale.svg", width = 10, height = 5)
pheatmap(r_matrixall, cluster_row=FALSE, cluster_col=FALSE, 
         color=inferno(75), show_colnames = FALSE, show_rownames = FALSE,)
dev.off()



## 2.3 Pairwise kinship per pop over all individuals ####

#Per pop
r_matrixfs <- fs.dosage(dos_microsatall, microsat_all[,1]) #meme info que relatedness mais average par pop
image(r_matrixfs$FsM)

#Per host
pop_bnh <- c(rep("bats", 59), rep("humans",311))
total_fs <- fs.dosage(dos_microsatall, pop_bnh)
R_Human <- (2*0.2956)/(1+0.5960) #0.3704261
R_Bat <- (2*0.5134)/(1+0.5960) #0.6433584


## 2.4 PCA over all individuals ####

pca <- indpca(microsat_all)

# % of genetic variance that explain each axis 
pca_eig <- 100*pca$ipca$eig/sum(pca$ipca$eig) #PC1 = 18.14% & PC2 = 13.87%
barplot(pca_eig[1:10])
sum(pca_eig[1:2])#32.01
sum(pca_eig[1:10])#73.33

plot(pca)

#Per host
Hosts <- as.data.frame(c(rep("bats",59),rep("humans", 311)))
pca_ALL <- cbind(Hosts, pca$ipca$li[, 1:2])
colHost <- c(rep("mediumpurple3",59),rep("orange2", 311))


svg("PCA_host.svg", width = 10, height = 5)
plot(pca_ALL[,2], pca_ALL[,3], col=colHost, pch=19, lwd=4,
     xlab="", ylab="", axes=F)
axis(side=2, at=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))
axis (side=1, at=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))
title(xlab="PC1 (18.14% explained variation)", 
      ylab="PC2 (13.87% explained variation)", cex.lab=1.5)
legend("topright", legend=c("Humans", "Bats"), col=c("orange2","mediumpurple3"),
       pch=19, cex=1.25, box.lty=0)
dev.off()




#########################################################
#### 3. Analysis with >5 individuals per pop #####
#########################################################
## 3.1 Download table ####

microsat_17pop<-read.fstat.data("./6.Microsat/Res_1indallpop/Microsat_17pop.dat")

## 3.2 Pairwise kinship per pop ####

#Dosage matrix
dos_microsat17pop <- fstat2dos(microsat_17pop[,2:ncol(microsat_17pop)], diploid=TRUE)
r_matrix17pop <- beta.dosage(dos_microsat17pop)
image(r_matrix17pop)

r_matrixfs_17pop <- fs.dosage(dos_microsat17pop, microsat_17pop[,1])
image(r_matrixfs_17pop$FsM) #FsM is pairwise kinship per pop and Fs2x2 is pairwise FsT


## 3.3 Genetic diversity ####
### 3.3.1 Heterozygosities and FIS ####
##Overall
Gendivall <- basic.stats(microsat_17pop)
Gendivall
write.csv("Gendiv_allind.csv", sep = ";")
wilcox.test(Gendivall$Hs[,1:2], Gendivall$Hs[,3:17], paired=TRUE)

##Bat-asso
GendivB <- basic.stats(microsat_17pop[1:59,])
GendivB
wilcox.test(GendivB$perloc[,1], GendivB$perloc[,2], paired=TRUE) #V=3, p=0.02439
boxplot(GendivB$perloc[,1], GendivB$perloc[,2])

##Human-asso
GendivH <- basic.stats(microsat_17pop[60:321,])
GendivH
wilcox.test(GendivH$perloc[,1], GendivH$perloc[,2], paired=TRUE) #V=1, p=0.01285
boxplot(GendivH$perloc[,1], GendivH$perloc[,2])

##Between hosts
boxplot(GendivB$perloc[,2], GendivH$perloc[,2])
wilcox.test(GendivB$perloc[,2], GendivH$perloc[,2], paired=TRUE) #V=50, p=0.01953


### 3.3.2 Allelic richness ####
##Overall
allele.count(microsat_17pop, diploid=TRUE)
allrich <- allelic.richness(microsat_17pop,min.n=NULL,diploid=TRUE)

##Bat-asso
allele.count(microsat_17pop[1:59,], diploid=TRUE)
allrichB <- allelic.richness(microsat_17pop[1:59,],min.n=NULL,diploid=TRUE)

##Human-asso
allele.count(microsat_17pop[60:321,], diploid=TRUE)
allrichH <- allelic.richness(microsat_17pop[60:321,],min.n=NULL,diploid=TRUE)


#Ar_test <- wilcox.test(, paired=TRUE)#V = , p-value = 


### 3.3.3 Private alleles ####
##Overall
nb.alleles(microsat_17pop,diploid=TRUE)

##PER POPULATION
POP = microsat_17pop[,1]
POP2 = c(POP,POP)

# initiate the table
TabPrivatePerPop = as.data.frame(matrix(0, nrow = length(unique(POP)), ncol = ncol(microsat_17pop)-1))
colnames(TabPrivatePerPop) = colnames(microsat_17pop)[2:ncol(microsat_17pop)]
rownames(TabPrivatePerPop) = unique(POP)

for(i in 1:ncol(TabPrivatePerPop)){
  alleles = c(substr(microsat_17pop[,i+1],1,3), substr(microsat_17pop[,i+1],4,6))
  for (Allele in levels(as.factor(alleles))){
    IsInPops = as.vector(na.omit(unique(POP2[alleles==Allele])))
    if(length(IsInPops)==1){
      TabPrivatePerPop[rownames(TabPrivatePerPop)==IsInPops,i] = TabPrivatePerPop[rownames(TabPrivatePerPop)==IsInPops,i]+1
    }
  }
}

colSums(TabPrivatePerPop)


##PER HOST

#add a host column
Host <- as.data.frame(c(rep("bats",59),rep("humans", 262)))
colnames(Host) = "Host"
microsat_host <- cbind(microsat_17pop, Host)
microsat_host <- microsat_host[,c(1,12, 2:11)]

HOST = microsat_host[,2]
HOST2 = c(HOST, HOST)

# initiate the table
TabPrivatePerHost = as.data.frame(matrix(0, nrow = length(unique(HOST)), ncol = ncol(microsat_host)-2))
colnames(TabPrivatePerHost) = colnames(microsat_host)[3:ncol(microsat_host)]
rownames(TabPrivatePerHost) = unique(HOST)

for(i in 1:ncol(TabPrivatePerHost)){
  alleles_host = c(substr(microsat_host[,i+2],1,3), substr(microsat_host[,i+2],4,6))
  for (Allele2 in levels(as.factor(alleles_host))){
    IsInHosts = as.vector(na.omit(unique(HOST2[alleles_host==Allele2])))
    if(length(IsInHosts)==1){
      TabPrivatePerHost[rownames(TabPrivatePerHost)==IsInHosts,i] = TabPrivatePerHost[rownames(TabPrivatePerHost)==IsInHosts,i]+1
    }
  }
}

rowSums(TabPrivatePerHost)

wilcox.test(as.numeric(TabPrivatePerHost[1,]), as.numeric(TabPrivatePerHost[2,]), paied=TRUE) #W=33, p=0.2078


### 3.3.4 Fis per sites ####

Fis_pop <- c()
for (pop in unique(microsat_17pop$Pop)){
  Pop=microsat_17pop[which(microsat_17pop$Pop==pop),]
  Gendiv = basic.stats(Pop)
  Fis_pop = c(Fis_pop, Gendiv$overall[9])
}
Fis_pop


#CI of FIS
boot.ppfis(microsat_17pop[1:59,], nboot=1000) #pop1&2 don't include 0
boot.ppfis(microsat_17pop[60:321,], nboot=1000) #pop2,13&14 don't include 0

boot.ppfis(microsat_host[,2:12], nboot=1000)#hosts don't include 0

FIS_human <- data.frame()
for(popH in length(microsat_17pop[60:321,]$Pop))
  

## 3.4 Pairwise FST ####

#Overall Fst
Gendiv_17pop <- basic.stats(microsat_17pop)
GendivB_17pop <- basic.stats(microsat_17pop[1:59,])
GendivH_17pop <- basic.stats(microsat_17pop[60:321,])

#Pairwise Fst
PFst_microsat17pop <- genet.dist(microsat_17pop, method="WC84") #faire une matrice pairwise FST sans NA
#Fs2x2 from fs.dosage should work too

Pfst_table <- read.csv("./6.Microsat/Res_17pop_paper/pFST_table.csv", header=T, sep=";")

#heatmap
ggplot(Pfst_table, aes(Pop1, Pop2, fill= value))+ 
  geom_tile()+
  scale_fill_viridis_c(option = "magma")


## Bonferroni significance

#to get p-values for pairwise FST from FSTAT, use test.g on all 17*16/2 population pairs. 
#You need 2720 (20*17*16/2) randomizations to achieve Bonferroni significance at 5% (100*17*16/2 at 1%)

pops<-as.integer(names(table(microsat_17pop[,1])))
pp.pval<-matrix(nrow=17,ncol=17)
for (i in 2:17) for (j in 1:(i-1)){
  print(paste0(pops[i],", ",pops[j]))
  pp.pval[i,j]<-test.g(microsat_17pop[microsat_17pop[,1]%in%c(pops[i],pops[j]),-1],
                       microsat_17pop[microsat_17pop[,1]%in%c(pops[i],pops[j]),1],
                       nperm=2720)$p.val
} 

pp.pval


## 3.5 Partial Mantel test ####

distgeo <- read.csv("./6.Microsat/Res_1indallpop/Coord_microsat_17pop.csv", header=TRUE, sep=";")

###coord suisses
distgeobis <- distgeo[,6:7]
row.names(distgeobis) <- distgeo$Pop_name
distgeobis_dist <- dist(distgeobis)

plot(distgeobis_dist, PFst_microsat17pop)
plot(log(distgeobis_dist), PFst_microsat17pop)

#### Geneva sampling
# 1-GE + 1-noGE + 0-reste
partial_GE <- read.csv("./PartialMantel_1GE.csv", header=TRUE, sep=";", row.names=1)
partial_GE <- as.matrix(partial_GE)
partial_GE_matrix <- as.dist(partial_GE, upper=F)

# # Is there IBD if I control for Geneva sampling?
# mantel.partial(PFst_microsat17pop, distgeobis_dist, partial_GE_matrix, permutations = 1000) #r=0.04727, p=0.3976
# 
# # Is there a Geneva effect if I control for the geographical distance?
# mantel.partial(PFst_microsat17pop, partial_GE_matrix, distgeobis_dist, permutations = 1000) #r=-0.1569, p=0.91808 


# List of rows & columns with GENEVA sites
GenevaRows = c(9:16)
GenevaColums = c(9:16)
as.matrix(PFst_microsat17pop)[GenevaRows,GenevaColums]
as.matrix(distgeobis_dist)[GenevaRows,GenevaColums]

# List of rows & columns without GENEVA sites
NoGenevaRows = c(1:8,17)
NoGenevaColums = c(1:8,17)
as.matrix(PFst_microsat17pop)[NoGenevaRows,NoGenevaColums]
as.matrix(distgeobis_dist)[NoGenevaRows,NoGenevaColums]

# List of rows & columns GE/noGE
GEnoGERows = c(9:16)
GEnoGEColums = c(1:8,17)
as.matrix(PFst_microsat17pop)[GEnoGERows,GEnoGEColums]
as.matrix(distgeobis_dist)[GEnoGERows,GEnoGEColums]

#### Host sampling
# 1-B + 1-H + 0-reste
partial_HB <- read.csv("./PartialMantel_1HB.csv", header=TRUE, sep=";", row.names=1)
partial_HB <- as.matrix(partial_HB)
partial_HB_matrix <- as.dist(partial_HB, upper=F)

# Is there a geographical effect if I control for the host sampling?
mantel.partial(PFst_microsat17pop, distgeobis_dist, partial_HB_matrix, permutations = 1000) #r=0.2795, p=0.006993 

#Is there a host effect if I control for the geographical distances ?
mantel.partial(PFst_microsat17pop, partial_HB_matrix, distgeobis_dist, permutations = 1000) #r=0.254, p=0.056943


##### GLS tests

datacheckGE_Host <- read.csv("./gls_test.csv", header=TRUE, sep=";")

#### Collinearity of the matrices
mantel(-partial_HB_matrix, distgeobis_dist, permutations=1000) #Geo-Host -> r = 0.6562, p = 0.00999
# correlated (bat-asso bugs are samples in the same areas and human-asso bugs also but both of them are sampled in different areas)
mantel(-partial_GE_matrix, distgeobis_dist, permutations=1000) #Geo-Site -> r = 0.5267, p = 0.000999
# correlated (bugs in GE are closer than bugs outside GE)
mantel(partial_HB_matrix, partial_GE_matrix) #Site-Host -> r = 0.00418, p = 0.49 
# not correlated (human-asso noGE does not compensate the fact that there is no bat-asso elsewhere)

###
# Test w/out geographical distances since correlated w/other matrices
###
mantel.partial(PFst_microsat17pop, partial_GE_matrix, partial_HB_matrix, permutations=1000)#r=-0.2124, p=0.998
summary(gls(GenetDist ~ scale(Geneva) + scale(Host), correlation=corMLPE(form=~ID1+ID2), data=datacheckGE_Host))

###
# Test w/ geographical distances on Human-asso only
###

# List of rows & columns with GENEVA sites
GenevaRowsH = c(7:14)
GenevaColumsH = c(7:14)
as.matrix(PFst_microsatH)[GenevaRowsH,GenevaColumsH]
as.matrix(distgeobisH_dist)[GenevaRowsH,GenevaColumsH]

# List of rows & columns without GENEVA sites
NoGenevaRowsH = c(1:6,15)
NoGenevaColumsH = c(1:6,15)
as.matrix(PFst_microsatH)[NoGenevaRowsH,NoGenevaColumsH]
as.matrix(distgeobisH_dist)[NoGenevaRowsH,NoGenevaColumsH]

# List of rows & columns GE/noGE
GEnoGERowsH = c(7:14)
GEnoGEColumsH = c(1:6,15)
as.matrix(PFst_microsatH)[GEnoGERowsH,GEnoGEColumsH]
as.matrix(distgeobisH_dist)[GEnoGERowsH,GEnoGEColumsH]


### IBD on Geneva only
plot(as.matrix(distgeobisH_dist)[GenevaRowsH,GenevaColumsH], as.matrix(PFst_microsatH)[GenevaRowsH,GenevaColumsH])
mantel(as.matrix(PFst_microsatH)[GenevaRowsH,GenevaColumsH], as.matrix(distgeobisH_dist)[GenevaRowsH,GenevaColumsH], permutations=1000)#r=-0.1151, p=0.59241 

plot(as.matrix(distgeobisH_dist)[GenevaRowsH,GenevaColumsH], as.matrix(PFst_microsatH)[GenevaRowsH,GenevaColumsH], type="p", pch=1, cex=1.5,
     xlab="", ylab="", axes =F)
points(distgeobisH_dist[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)], 
       PFst_microsatH[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)], 
       col=c("orange1"), pch=19, cex=1.5)#GE/GE - olivedrab3 - #606C38
axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(side=1, at=c(0, 1000, 2000, 3000, 4000, 5000))
title(xlab="Geographical distance (m)", 
      ylab="Genetic distances (FST)", cex.lab=1.5)


###IBD on noGeneva only
plot(as.matrix(distgeobis_dist)[NoGenevaRows,NoGenevaColums], as.matrix(PFst_microsat17pop)[NoGenevaRows,NoGenevaColums])
mantel(as.matrix(PFst_microsat17pop)[NoGenevaRows,NoGenevaColums], as.matrix(distgeobis_dist)[NoGenevaRows,NoGenevaColums], permutations=1000)#r=0.08724, p=0.29471


###
# Plot in Human-asso only
###

PFst_microsatH <- genet.dist(microsat_17pop[60:321,], method="WC84") #faire une matrice pairwise FST sans NA
distgeobisH <- distgeo[3:17,6:7]
row.names(distgeobisH) <- distgeo[3:17,]$Pop_name
distgeobisH_dist <- dist(distgeobisH)

plot(distgeobisH_dist, PFst_microsatH, type="p", pch=1, cex=1.5,
     xlab="", ylab="", axes =F)
points(distgeobisH_dist[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)], 
       PFst_microsatH[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)], 
       col=c("orange1"), pch=19, cex=1.5)#GE/GE - olivedrab3 - #606C38
points(distgeobisH_dist[c(1:5, 14:18, 27:30, 39:41, 50, 51, 60, 69)], 
       PFst_microsatH[c(1:5, 14:18, 27:30, 39:41, 50, 51, 60, 69)], 
       col=c("firebrick"), pch=19, cex=1.5)#noGE/noGE - steelblue1 - #283618
points(distgeobisH_dist[c(6:13, 19:26, 31:38, 42:49, 52:59, 61:68, 77, 84, 90, 95, 99, 102, 104, 105)], 
       PFst_microsatH[c(6:13, 19:26, 31:38, 42:49, 52:59, 61:68, 77, 84, 90, 95, 99, 102, 104, 105)], 
       col=c("gold1"), pch=19, cex=1.5)#GE/noGE - mediumorchid
axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(side=1, at=c(0, 25000, 50000, 75000, 100000, 125000))
title(xlab="Geographical distance (m)", 
      ylab="Genetic distances (FST)", cex.lab=1.5)


#LOG_TRANSFORMED
plot(log(distgeobisH_dist), PFst_microsatH, type="p", pch=1, cex=1.5,
     xlab="", ylab="", axes =F)
points(log(distgeobisH_dist[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)]), 
       PFst_microsatH[c(70:76, 78:83, 85:89, 91:94, 96:98, 100, 101, 103)], 
       col=c("orange1"), pch=19, cex=1.5)#GE/GE - olivedrab3 - #606C38
points(log(distgeobisH_dist[c(1:5, 14:18, 27:30, 39:41, 50, 51, 60, 69)]), 
       PFst_microsatH[c(1:5, 14:18, 27:30, 39:41, 50, 51, 60, 69)], 
       col=c("firebrick"), pch=19, cex=1.5)#noGE/noGE - steelblue1 - #283618
points(log(distgeobisH_dist[c(6:13, 19:26, 31:38, 42:49, 52:59, 61:68, 77, 84, 90, 95, 99, 102, 104, 105)]), 
       PFst_microsatH[c(6:13, 19:26, 31:38, 42:49, 52:59, 61:68, 77, 84, 90, 95, 99, 102, 104, 105)], 
       col=c("gold1"), pch=19, cex=1.5)#GE/noGE - mediumorchid
axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(side=1, at=c(5, 6, 7, 8, 9, 10, 11))
title(xlab="Log of Geographical distance (m)", 
      ylab="Genetic distances (FST)", cex.lab=1.5)



## 3.6 Hierarchical analysis of variance ####

#Add Host column
Host <- as.data.frame(c(rep("bats",59),rep("humans", 262)))
colnames(Host) = "Host"
microsat_17pop2 <- cbind(microsat_17pop, Host)
microsat_17pop2 <- microsat_17pop2[,c(1,12, 2:11)]

#Add Individuals column
Ind <- as.data.frame(c(1:321))
colnames(Ind) = "Ind"
microsat_17pop3 <- cbind(microsat_17pop, Ind)
microsat_17pop3 <- microsat_17pop3[,c(12,1, 2:11)]


#Estimation of the variance components
#For each locus
#Overall
#Matrix of hierarchical Fstat
test_varcomp <- varcomp.glob(levels=microsat_17pop2[, c(2,1)], loci=microsat_17pop2[, -c(1,2)])

##Test effect of Host
test_betweenHost <- test.between(microsat_17pop2[, -c(1,2)],rand.unit = microsat_17pop2$Pop,test = microsat_17pop2$Host,nperm = 10000)
#p=0.19

##Test effect of populations within host
test_withinHost_betweenPop <- test.within(microsat_17pop2[, -c(1,2)],test.lev = microsat_17pop2$Pop,within = microsat_17pop2$Host,nperm =1000)
#p=0.001

##Test effect of individuals within populations
test_withinPop_betweenInd <- test.within(microsat_17pop3[, -c(1,2)],test.lev = microsat_17pop3$Ind,within = microsat_17pop3$Pop,nperm =10000)
#p=1

#bootstrap values for variance components
boot.vc(levels=microsat_17pop2[, c(2,1)], loci=microsat_17pop2[, -c(1,2)], nboot=100)


#########################################################
#### 4. Analysis on human dataset only #####
#########################################################
## 4.1 Download table ####
microsatH <- microsat_all[-(1:59),]


## 4.2 Pairwise kinship per pop ####

#Dosage matrix
dos_microsatH <- fstat2dos(microsatH[,2:ncol(microsatH)], diploid=TRUE)
r_matrixH <- beta.dosage(dos_microsatH)
svg("./image_kinship_HUMAN.svg", width = 10, height = 5)
image(r_matrixH)
dev.off()

r_matrixfs_H <- fs.dosage(dos_microsatH, microsatH[,1]) #meme info que relatedness mais average par pop
svg("./image_kinshipaveraged_HUMAN.svg", width = 10, height = 5)
image(r_matrixfs_H$FsM)
dev.off()


## 4.3 PCA all human-associated individuals ####

pcaH <- indpca(microsatH)

# % of genetic variance that explain each axis 
pcaH_eig <- 100*pcaH$ipca$eig/sum(pcaH$ipca$eig) #PC1 = 21.02% & PC2 = 11.70%
barplot(pcaH_eig[1:10])
sum(pcaH_eig[1:2])#32.73
sum(pcaH_eig[1:10])#76.56

plot(pcaH)

#Per pop ellpise
s.class(pcaH$ipca$li, as.factor(microsatH$Pop))
svg("./PCA_human_ellispe.svg", width = 10, height = 5)
s.class(pcaH$ipca$li, as.factor(microsatH$Pop), grid=FALSE)
dev.off()


## 4.4 Pairwise FST per human shelters ####

PFst_microsatH <- genet.dist(microsatH, method="WC84") #faire une matrice pairwise FST sans NA
pheatmap(PFst_microsatH)
