
                                                                
                                                                                            ############################################
                                                                                            ############################################
                                                                                            ############ GENOMIC PIPELINE ##############
                                                                                            ############################################
                                                                                            ############################################



###########################################################
################# 1- GENOME PROCESSING ####################
###########################################################


########################
# 1) INDEXING
########################


# téléchargement du génome : NCBI --> Myotis daubentonii --> FTP --> copie le lien : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/259/705/GCF_963259705.1_mMyoDau2.1/GCF_963259705.1_mMyoDau2.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/259/705/GCF_963259705.1_mMyoDau2.1/GCF_963259705.1_mMyoDau2.1_cds_from_genomic.fna.gz

# atassel/genome/index.sh


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name index
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/11.4.0
module load bwa/0.7.17

gunzip -c /work/FAC/FBM/DEE/pchriste/default/atassel/genome/GCF_963259705.1_mMyoDau2.1_genomic.fna.gz > /work/FAC/FBM/DEE/pchriste/default/atassel/genome/MyoDau2.1_genomic.fna
bwa index /work/FAC/FBM/DEE/pchriste/default/atassel/genome/MyoDau2.1_genomic.fna -p MyoDau2.1_genomic.fna

# Submitted batch job 47252041
# Slurm Job_id=47252041 Name=index Began, Queued time 00:00:17



###########################################################
################ 2- LIBRARY PROCESSING ####################
###########################################################


###############################
########### SPOCK 103 #########
###############################


########################
# 1) TRIMMING
########################


# atassel/v2/trimming/trimming.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name trimming
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

set -x
set -e

module load gcc/11.4.0
module load bbmap/39.01

output_lane1="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1"
output_lane2="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2"
adaptors="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/adapters.fa"

for f1 in /work/FAC/FBM/DEE/pchriste/default/atassel/rawdata/lane1/*.fastq.gz
do
    outname1=$(basename "$f1" .fastq.gz)
    echo " processing " $outname1

    bbduk.sh in=${f1} out=${output_lane1}/${outname1}.bbduk.fastq.gz ref=${adaptors} ktrim=r k=20 mink=10 hdist=1 tbo overwrite=true bhist=${output_lane1}/bhist_${outname1}.txt qhist=${output_lane1}/qhist_${outname1}.txt aqhist=${output_lane1}/aqhist_${outname1}.txt lhist=${output_lane1}/lhist_${outname1}.txt >& ${output_lane1}/${outname1}.bbduk.log
    echo " Done"

done

for f2 in /work/FAC/FBM/DEE/pchriste/default/atassel/rawdata/lane2/*.fastq.gz
do
    outname2=$(basename "$f2" .fastq.gz)
    echo " processing " $outname2 

    bbduk.sh in=${f2} out=${output_lane2}/${outname2}.bbduk.fastq.gz ref=${adaptors} ktrim=r k=20 mink=10 hdist=1 tbo overwrite=true bhist=${output_lane2}/bhist_${outname2}.txt qhist=${output_lane2}/qhist_${outname2}.txt aqhist=${output_lane2}/aqhist_${outname2}.txt lhist=${output_lane2}/lhist_${outname2}.txt >& ${output_lane2}/${outname2}.bbduk.log
    echo " Done"

done



# Submitted batch job 47299515
# Slurm Job_id=47299515 Name=trimming Began, Queued time 00:00:01
# Slurm Job_id=47299515 Name=trimming Ended, Run time 00:17:31, COMPLETED, ExitCode 0




#######
### QUAL STAT
#######

# atassel/trimming/qual_stat.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name QualStat
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE





module load r/4.3.2


COMPO="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/bhist_BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5.txt"
QUAL="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/qhist_BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5.txt"
LENGTH="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5_reads.length.txt" 

  gzip -dc /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5.bbduk.fastq.gz | awk '(FNR%4)==2{l[length($0)]+=1}END{for (i in l) print i, l[i];}' > /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5_reads.length.txt




R --vanilla <<EOF


# GRAPHIQUE QUALITE NUCLEOTIDE

table_compo <- read.table("$COMPO", header=F, row.names=1)
prop <- t(as.matrix(table_compo))
table_quality <- read.table("$QUAL", header=F, row.names=1)
pdf(paste0("/work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/qual_stat/","BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5","_quality_nucleotide.pdf"), width=10, height=15)
layout(matrix(c(1:2), ncol = 1))
plot(table_quality[,1] ~ as.numeric(row.names(table_quality)), xlim=c(0,(nrow(table_quality)+3)), ylim=c(0,(max(table_quality[,1])+10)), col="blue", pch=20, xlab="Nucleotide position along the read", ylab="Quality score", main="BAT_630_568_531_552__G6_L1_R1_001_CL6tjXKbW3H5")
  

# GRAPHIQUE COMPOSITION EN NUCLEOTIDE LE LONG DU READ

barplot(prop, col=c("blue", "white", "green", "pink", "black"), xlab="nucleotide composition along the read")
legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))


# GRAPHIQUE TAILLE DES READS

table_compo_2 <- read.table("$LENGTH", header=F, row.names=1)
tutu <- as.character(row.names(table_compo_2)[order(as.numeric(row.names(table_compo_2)))])
table_compo_2 <- as.data.frame(table_compo_2[order(as.numeric(row.names(table_compo_2))),])
row.names(table_compo_2)<-tutu
par(mar=c(8, 6, 4, 6))
b = barplot(table_compo_2[,1], col = "blue", xlab="Reads length after adapters trimming", names.arg = row.names(table_compo_2), las = 3, cex.names=0.5, axes = F)
axis(2, las=2)
mtext("Nb reads", side = 2, line = 4, cex.lab = 1, las = 3)


dev.off()


EOF



########################
# 2) DEMULTIPLEXING
########################

#######
### LANE 1
#######

# atassel/v2/demultiplexed/demultiplex_lane1.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name demultiplex_lane1
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53



process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*A9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_A9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*B9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_B9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*C9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_C9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*D9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_D9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*E9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_E9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*F9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_F9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*G9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_G9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane1/*H9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/process_radtags.lane1_H9.log

# Submitted batch job 47301644
# Slurm Job_id=47301644 Name=demultiplex_lane1 Began, Queued time 00:00:01
# Slurm Job_id=47301644 Name=demultiplex_lane1 Ended, Run time 00:53:48, COMPLETED, ExitCode 0

# OUTPUT
1373761 total sequences
 132965 barcode not found drops (9.7%)
  45676 low quality read drops (3.3%)
   8921 RAD cutsite not found drops (0.6%)
1186199 retained reads (86.3%)

## nbr of reads : 143 704 128


#######
### LANE 2
#######

# atassel/v2/demultiplexed/demultiplex_lane2.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name demultiplex_lane2
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*A9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_A9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*B9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_B9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*C9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_C9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*D9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_D9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*E9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_E9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*F9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_F9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*G9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_G9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v2/trimming/lane2/*H9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/process_radtags.lane2_H9.log



# Submitted batch job 47304894
# Slurm Job_id=47304894 Name=demultiplex_lane2 Began, Queued time 00:00:01
# Slurm Job_id=47304894 Name=demultiplex_lane2 Ended, Run time 00:50:50, COMPLETED, ExitCode 0

# OUTPUT 
1438594 total sequences
 242342 barcode not found drops (16.8%)
  56562 low quality read drops (3.9%)
   8915 RAD cutsite not found drops (0.6%)
1130775 retained reads (78.6%)



## nbr reads = 140 112 435

########################
# 3) MERGE LANE 1 ET 2
########################


#######
### LANE 1
#######


for seq in /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane1/*.fq.gz
do
    outname1=$(basename ${seq} .fq.gz)
    echo " processing " ${outname1}

    cp ${seq} /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/
    mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/${outname1}.fq.gz /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/${outname1}_L1.fq.gz

done



#######
### LANE 2
#######

for seq in /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/lane2/*.fq.gz
do
    outname1=$(basename ${seq} .fq.gz)
    echo " processing " ${outname1}

    cp ${seq} /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/
    mv /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/${outname1}.fq.gz /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/${outname1}_L2.fq.gz

done




#######
### MERGING
####### 


# atassel/v2/merged/merging.sh


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name merged
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for seq1 in /work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/*_L1.fq.gz
do
    outname=$(basename "${seq1}" _L1.fq.gz)
    
    seq2="/work/FAC/FBM/DEE/pchriste/default/atassel/v2/demultiplexed/all/${outname}_L2.fq.gz"
    echo "processing " ${outname}    
    cat <(zcat "$seq1") <(zcat "$seq2") | gzip > "/work/FAC/FBM/DEE/pchriste/default/atassel/v2/merged/${outname}_merged.fq.gz" 
    echo ${outname} " merging finished"

done


# Submitted batch job 47305424
# Slurm Job_id=47305424 Name=merged Began, Queued time 00:00:01
# Slurm Job_id=47305424 Name=merged Ended, Run time 01:31:00, COMPLETED, ExitCode 0
                                                        



###############################
########### SPOCK 108 #########
###############################


########################
# 1) TRIMMING
########################


# atassel/v3/trimming/trimming.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name trimming
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

set -x
set -e

module load gcc/11.4.0
module load bbmap/39.01

output_lane1="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1"
output_lane2="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2"
adaptors="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/adapters.fa"

for f1 in /work/FAC/FBM/DEE/pchriste/default/atassel/rawdata_SPOCK103/lane1/*.fastq.gz
do
    outname1=$(basename "$f1" .fastq.gz)
    echo " processing " $outname1

    bbduk.sh in=${f1} out=${output_lane1}/${outname1}.bbduk.fastq.gz ref=${adaptors} ktrim=r k=20 mink=10 hdist=1 tbo overwrite=true bhist=${output_lane1}/bhist_${outname1}.txt qhist=${output_lane1}/qhist_${outname1}.txt aqhist=${output_lane1}/aqhist_${outname1}.txt lhist=${output_lane1}/lhist_${outname1}.txt >& ${output_lane1}/${outname1}.bbduk.log
    echo " Done"

done



# LANE 2

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 4:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name trimming
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

set -x
set -e

module load gcc/11.4.0
module load bbmap/39.01


output_lane2="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis"
adaptors="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/adapters.fa"

for f2 in /work/FAC/FBM/DEE/pchriste/default/atassel/rawdata_SPOCK103/lane2/*.fastq.gz
do
    outname2=$(basename "$f2" .fastq.gz)
    echo " processing " $outname2 

    bbduk.sh in=${f2} out=${output_lane2}/${outname2}.bbduk.fastq.gz ref=${adaptors} ktrim=r k=20 mink=10 hdist=1 tbo overwrite=true bhist=${output_lane2}/bhist_${outname2}.txt qhist=${output_lane2}/qhist_${outname2}.txt aqhist=${output_lane2}/aqhist_${outname2}.txt lhist=${output_lane2}/lhist_${outname2}.txt >& ${output_lane2}/${outname2}.bbduk.log
    echo " Done"

done



# Submitted batch job 47339548
# Slurm Job_id=47339548 Name=trimming Began, Queued time 00:00:00
# Slurm Job_id=47339548 Name=trimming Ended, Run time 00:17:07, COMPLETED, ExitCode 0




#######
### QUAL STAT
#######

# atassel/v3/trimming/qual_stat.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name QualStat
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE





module load r/4.3.2


COMPO="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/bhist_BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6.txt"
QUAL="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/qhist_BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6.txt"
LENGTH="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6_reads.length.txt" 

  gzip -dc /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6.bbduk.fastq.gz | awk '(FNR%4)==2{l[length($0)]+=1}END{for (i in l) print i, l[i];}' > /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6_reads.length.txt




R --vanilla <<EOF


# GRAPHIQUE QUALITE NUCLEOTIDE

table_compo <- read.table("$COMPO", header=F, row.names=1)
prop <- t(as.matrix(table_compo))
table_quality <- read.table("$QUAL", header=F, row.names=1)
pdf(paste0("/work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/qual_stat/","BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6","_quality_nucleotide.pdf"), width=10, height=15)
layout(matrix(c(1:2), ncol = 1))
plot(table_quality[,1] ~ as.numeric(row.names(table_quality)), xlim=c(0,(nrow(table_quality)+3)), ylim=c(0,(max(table_quality[,1])+10)), col="blue", pch=20, xlab="Nucleotide position along the read", ylab="Quality score", main="BAT_630_568_531_552__G6_L1_R1_001_1tGtFK6TfKi6")
  

# GRAPHIQUE COMPOSITION EN NUCLEOTIDE LE LONG DU READ

barplot(prop, col=c("blue", "white", "green", "pink", "black"), xlab="nucleotide composition along the read")
legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))


# GRAPHIQUE TAILLE DES READS

table_compo_2 <- read.table("$LENGTH", header=F, row.names=1)
tutu <- as.character(row.names(table_compo_2)[order(as.numeric(row.names(table_compo_2)))])
table_compo_2 <- as.data.frame(table_compo_2[order(as.numeric(row.names(table_compo_2))),])
row.names(table_compo_2)<-tutu
par(mar=c(8, 6, 4, 6))
b = barplot(table_compo_2[,1], col = "blue", xlab="Reads length after adapters trimming", names.arg = row.names(table_compo_2), las = 3, cex.names=0.5, axes = F)
axis(2, las=2)
mtext("Nb reads", side = 2, line = 4, cex.lab = 1, las = 3)


dev.off()


EOF


# Submitted batch job 47331019
# Slurm Job_id=47331019 Name=QualStat Began, Queued time 00:00:01
# Slurm Job_id=47331019 Name=QualStat Ended, Run time 00:00:06, COMPLETED, ExitCode 0




#######
### POUR AVOIR LE NBR DE READ LANE 1
#######

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name BAT_lane1_nbreads
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE


zcat *.fastq.gz | grep @ | wc -l > nbreads_lane1.txt



# Slurm Job_id=47310496 Name=BAT_lane1_nbreads Began, Queued time 00:00:01
# Slurm Job_id=47310496 Name=BAT_lane1_nbreads Ended, Run time 00:13:40, COMPLETED, ExitCode 0

# OUTPUT: 512587474



#######
### POUR AVOIR LE NBR DE READ LANE 2
#######

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name BAT_lane2_nbreads
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE


zcat *.fastq.gz | grep @ | wc -l > nbreads_lane2.txt


# Slurm Job_id=47310520 Name=BAT_lane2_nbreads Began, Queued time 00:00:00
# Slurm Job_id=47310520 Name=BAT_lane2_nbreads Ended, Run time 00:13:55, COMPLETED, ExitCode 0

# OUTPUT : 527675731

### TOTAL : 1 040 263 205




########################
# 2) DEMULTIPLEXING
########################

#######
### LANE 1
#######

# atassel/v3/demultiplexed/lane1/dem_lane1.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name demultiplex_lane1
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53


process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*A9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_A9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*B9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_B9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*C9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_C9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*D9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_D9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*E9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_E9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*F9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_F9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*G9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_G9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane1/*H9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/process_radtags.lane1_H9.log


# Submitted batch job 47327985
# Slurm Job_id=47327985 Name=demultiplex_lane1 Began, Queued time 00:00:01
# Slurm Job_id=47327985 Name=demultiplex_lane1 Ended, Run time 03:31:18, COMPLETED, ExitCode 0

# OUTPUT 
4341471 total sequences
  73200 barcode not found drops (1.7%)
   6825 low quality read drops (0.2%)
  30064 RAD cutsite not found drops (0.7%)
4231382 retained reads (97.5%)



#######
### LANE 2
#######

# atassel/v3/demultiplexed/dem_lane2.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 4:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name demultiplex_lane2
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53


process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*A9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_A9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_A9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*B9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_B9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_B9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*C9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_C9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_C9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*D9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_D9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_D9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*E9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_E9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_E9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*F9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_F9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_F9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*G9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_G9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_G9.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H1_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H1.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H1.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H10_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H10.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H10.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H11_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H11.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H11.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H12_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H12.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H12.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H2_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H2.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H2.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H3_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H3.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H3.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H4_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H4.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H4.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H5_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H5.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H5.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H6_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H6.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H6.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H7_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H7.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H7.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H8_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H8.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H8.log
process_radtags -f /work/FAC/FBM/DEE/pchriste/default/atassel/v3/trimming/lane2_bis/*H9_*.fastq.gz -b /work/FAC/FBM/DEE/pchriste/default/atassel/barcode_txt_files_BAT/barcode_txt_files_BAT_H9.txt -o /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2 -crq -e ecoRI --renz_2 mseI; mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_bis.log  /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/process_radtags.lane2_H9.log

# Submitted batch job 47339894 
# Slurm Job_id=47339894  Name=demultiplex_lane2 Began, Queued time 00:00:00
# Slurm Job_id=47339894  Name=demultiplex_lane1 Ended, Run time 00:16:46, COMPLETED, ExitCode 0


# OUTPUT : 
4352735 total sequences
  73937 barcode not found drops (1.7%)
  75747 low quality read drops (1.7%)
  30562 RAD cutsite not found drops (0.7%)
4172489 retained reads (95.9%)



########################
# 3) MERGE LANE 1 ET 2
########################


#######
### LANE 1
#######


for seq in /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane1/*.fq.gz
do
    outname1=$(basename ${seq} .fq.gz)
    echo " processing " ${outname1}

    cp ${seq} /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/
    mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/${outname1}.fq.gz /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/${outname1}_L1.fq.gz

done



#######
### LANE 2
#######

for seq in /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/lane2/*.fq.gz
do
    outname1=$(basename ${seq} .fq.gz)
    echo " processing " ${outname1}

    cp ${seq} /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/
    mv /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/${outname1}.fq.gz /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/${outname1}_L2.fq.gz

done




#######
### MERGING
####### 


# atassel/v3/merged/merging.sh


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 5:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name merged
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for seq1 in /work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/*_L1.fq.gz
do
    outname=$(basename "${seq1}" _L1.fq.gz)
    
    seq2="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/demultiplexed/all/${outname}_L2.fq.gz"
    echo "processing " ${outname}    
    cat <(zcat "$seq1") <(zcat "$seq2") | gzip > "/work/FAC/FBM/DEE/pchriste/default/atassel/v3/merged/${outname}_merged.fq.gz" 
    echo ${outname} " merging finished"

done


# Submitted batch job 47353753
# Slurm Job_id=47353753 Name=merged Began, Queued time 00:00:03




#######
### MERGING SPOCK 103 et SPOCK 108
####### 


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 30:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name merged
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for seq1 in /work/FAC/FBM/DEE/pchriste/default/atassel/v2/merged/*_merged.fq.gz
do
    outname=$(basename "${seq1}" _merged.fq.gz)
    
    seq2="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/merged/${outname}_merged.fq.gz"
    echo "processing " ${outname}    
    cat <(zcat "$seq1") <(zcat "$seq2") | gzip > "/work/FAC/FBM/DEE/pchriste/default/atassel/v3/Merged_ALL/${outname}_merged.fq.gz" 
    echo ${outname} " merging finished"

done


# Submitted batch job 47365521
# Slurm Job_id=47365521 Name=merged Began, Queued time 00:00:01
# Slurm Job_id=47365521 Name=merged Ended, Run time 06:41:16, COMPLETED, ExitCode 0



########################
# 4) MAPPING
########################


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name mapping
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE


module load gcc/11.4.0
module load bwa/0.7.17
module load samtools/1.17

genome="/work/FAC/FBM/DEE/pchriste/default/atassel/genome/MyoDau2.1_genomic.fna"
path_reads="/work/FAC/FBM/DEE/pchriste/default/atassel/v3/Merged_ALL"
outpath1="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped"
path_map_read="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped"
outpath2="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/output_stat"


for ind in ${path_reads}/*_merged.fq.gz
do 
    outname1=$(basename "$ind" _merged.fq.gz)
    echo " processing for bwa mem" ${outname1}
    bwa mem -t 16 ${genome} ${ind} | samtools view -@ 16 -Sb | samtools sort -@ 16 -O BAM > ${outpath1}/${outname1}.bam
    echo " done bwa mem"

done


for read in ${path_map_read}/*.bam
do 
    outname2=$(basename "$read" .bam)
    echo " processing for flagstat" ${outname2}
    samtools flagstat ${read} > ${outpath2}/${outname2}_flag.txt
    echo " done flagstat"
    

done


# Submitted batch job 47722651
# Slurm Job_id=47722651 Name=mapping Began, Queued time 00:00:23
# Slurm Job_id=47722651 Name=mapping Ended, Run time 23:54:25, COMPLETED, ExitCode 0


##########################################################################################
#################### 3 - GENOTYPING PIPELINE AVEC POPULATION (1) 0.8 #######################
##########################################################################################



########################
# 1) GSTACKS
########################


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 15:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name gstacks
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE


module load gcc/11.4.0
module load stacks/2.53


path_map_read="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/"
outpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/"
pop_map="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/pop_map.txt"


gstacks -I ${path_map_read} -O ${outpath} -M ${pop_map}
echo "processing"


# Submitted batch job 47782387
# Slurm Job_id=47782387 Name=gstacks Began, Queued time 00:00:01
# Slurm Job_id=47782387 Name=gstacks Ended, Run time 02:54:05, COMPLETED, ExitCode 0


#OUTPUT:
Genotyped 1070589 loci:
  effective per-sample coverage: mean=6.2x, stdev=2.1x, min=1.0x, max=12.7x
  mean number of sites per locus: 131.2
  a consistent phasing was found for 9150783 of out 9341566 (98.0%) diploid loci needing phasing



########################
# 2) POPULATIONS
########################


## on refait avec 0.8 pour être consistent

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name populations
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE


module load gcc/11.4.0
module load stacks/2.53

inpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/"
outpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/"
pop_map="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/pop_map.txt"

populations -P ${inpath} --popmap ${pop_map} -O ${outpath} -p 18 -r 0.8 -R 0.8 -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp


# Submitted batch job 449983010
# Slurm Job_id=49983010 Name=populations_0.8 Began, Queued time 00:00:31
# Slurm Job_id=49983010 Name=populations_0.8 Ended, Run time 01:13:56, COMPLETED, ExitCode 0


# OUTPUT
Removed 955000 loci that did not pass sample/population constraints from 1070589 loci.
Kept 115589 loci, composed of 16679051 sites; 130813 of those sites were filtered, 96305 variant sites remained.
    16642923 genomic sites, of which 36128 were covered by multiple loci (0.2%).
Mean genotyped sites per locus: 144.30bp (stderr 0.01).


########################
# 3) VCF TOOLS
########################

#../atassel/SNP_filtering/population_0.8

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf /work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/min_cov_5


# Resultat
After filtering, kept 395 out of 395 Individuals
Outputting VCF file...
After filtering, kept 96305 out of a possible 96305 Sites
Run Time = 61.00 seconds



#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8

Changing_allFIELDS_from_doubleDots_Geno(){

VCFin=$1
VCFout=$2

grep '#' ${VCFin} > ${VCFout}

grep -v '#' ${VCFin} | awk '{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t"
for (i=10; i<=NF; i++){
split($i,a,":")
if(i < NF){
if(a[1] == "./.") {
printf a[1]":.:.:.:.\t"}
else {
printf $i "\t"}}
else {
if(a[1] == "./.") {
print a[1]":.:.:.:.\t"}
else {
print $i"\t"}}}}' >> ${VCFout}

}

#Copy and Paste the ENTIRE FUNCTION
#Then just use it with FIRST RGUMENT = INPUT VCF; second ARGUMENT = OUTPUT VCF --> Changing_allFIELDS_from_doubleDots_Geno input_VCF.vcf output_VCF.vcf


Changing_allFIELDS_from_doubleDots_Geno min_cov_5.recode.vcf min_cov_5_clean.recode.vcf



########################
# 4) SNP sharing
########################


#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8

module load gcc/11.4.0
module load vcftools/0.1.16


#######
### SNPS SHARING - 85%
#######

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.85 --recode --recode-INFO-all --out min_cov_5_clean_85miss


# Resultat
After filtering, kept 395 out of 395 Individuals
Outputting VCF file...
After filtering, kept 52262 out of a possible 96305 Sites
Run Time = 34.00 seconds





########################
# 5) QUALITY CHECK
########################

module load gcc/11.4.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/min_cov_5_clean_85miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/min_cov_5_clean_85miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het

## il y a 29 individus qui montre plus de 30% de missing data



########################
# POPULATIONS (2)
########################

# on refait population mais en enlevant les individus qui ont plus de 30% de missing data 


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --job-name populations
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE


module load gcc/11.4.0
module load stacks/2.53

inpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/"
outpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/population_2eme"
pop_map="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/population_2eme/pop_map.txt"

populations -P ${inpath} --popmap ${pop_map} -O ${outpath} -p 18 -r 0.8 -R 0.8 -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp


# Submitted batch job 51027974
# Slurm Job_id=51027974 Name=populations Began, Queued time 00:00:28
# Slurm Job_id=51027974 Name=populations Ended, Run time 00:59:36, COMPLETED, ExitCode 0



#OUTPUT:
Removed 935384 loci that did not pass sample/population constraints from 1070589 loci.
Kept 135205 loci, composed of 19508347 sites; 114417 of those sites were filtered, 113439 variant sites remained.
    19459458 genomic sites, of which 48889 were covered by multiple loci (0.3%).
Mean genotyped sites per locus: 144.29bp (stderr 0.01).



########################
# VCF TOOLS (2)
########################


module load gcc/13.2.0
module load vcftools/0.1.16

vcftools --vcf /work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/population_2eme/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_2eme/min_cov_5

# Resultat
After filtering, kept 366 out of 366 Individuals
Outputting VCF file...
After filtering, kept 113439 out of a possible 113439 Sites
Run Time = 68.00 seconds


#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_2eme

Changing_allFIELDS_from_doubleDots_Geno(){

VCFin=$1
VCFout=$2

grep '#' ${VCFin} > ${VCFout}

grep -v '#' ${VCFin} | awk '{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t"
for (i=10; i<=NF; i++){
split($i,a,":")
if(i < NF){
if(a[1] == "./.") {
printf a[1]":.:.:.:.\t"}
else {
printf $i "\t"}}
else {
if(a[1] == "./.") {
print a[1]":.:.:.:.\t"}
else {
print $i"\t"}}}}' >> ${VCFout}

}

#Copy and Paste the ENTIRE FUNCTION
#Then just use it with FIRST RGUMENT = INPUT VCF; second ARGUMENT = OUTPUT VCF --> Changing_allFIELDS_from_doubleDots_Geno input_VCF.vcf output_VCF.vcf


Changing_allFIELDS_from_doubleDots_Geno min_cov_5.recode.vcf min_cov_5_clean.recode.vcf




##############################
# SNP sharing (2) MAF = 0.005
##############################


#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_2eme/maf_0.005

module load gcc/11.4.0
module load vcftools/0.1.16


#######
### SNPS SHARING - 90%
#######

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.85 --recode --recode-INFO-all --out min_cov_5_clean_85miss

# Resultat
After filtering, kept 366 out of 366 Individuals
Outputting VCF file...
After filtering, kept 65702 out of a possible 113439 Sites
Run Time = 40.00 seconds


## quality check pour SNP sharing 90

module load gcc/13.2.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_2eme/min_cov_5_clean_85miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_2eme/min_cov_5_clean_85miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het


### Encore 7 individus avec des missings data




########################
# POPULATIONS (3)
########################

# on refait population mais en enlevant les individus qui ont plus de 30% de missing data 


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --job-name populations
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE


module load gcc/13.2.0
module load stacks/2.53

inpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/"
outpath="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/population_3eme"
pop_map="/work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/pop_map_2.txt"

populations -P ${inpath} --popmap ${pop_map} -O ${outpath} -p 18 -r 0.8 -R 0.8 -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp


# Submitted batch job 51031379
# Slurm Job_id=51031379 Name=populations Began, Queued time 00:00:19
# Slurm Job_id=51031379 Name=populations Ended, Run time 00:53:22, COMPLETED, ExitCode 0


# OUTPUT
Removed 939773 loci that did not pass sample/population constraints from 1070589 loci.
Kept 130816 loci, composed of 18875824 sites; 128098 of those sites were filtered, 107274 variant sites remained.
    18830267 genomic sites, of which 45557 were covered by multiple loci (0.2%).
Mean genotyped sites per locus: 144.29bp (stderr 0.01).



########################
# VCF TOOLS (3)
########################


module load gcc/13.2.0
module load vcftools/0.1.16

vcftools --vcf /work/FAC/FBM/DEE/pchriste/default/atassel/mapped/gstacks/population_0.8/final/population_3eme/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5

# Resultat
After filtering, kept 359 out of 359 Individuals
Outputting VCF file...
After filtering, kept 107274 out of a possible 107274 Sites
Run Time = 64.00 seconds


#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_2eme

Changing_allFIELDS_from_doubleDots_Geno(){

VCFin=$1
VCFout=$2

grep '#' ${VCFin} > ${VCFout}

grep -v '#' ${VCFin} | awk '{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t"
for (i=10; i<=NF; i++){
split($i,a,":")
if(i < NF){
if(a[1] == "./.") {
printf a[1]":.:.:.:.\t"}
else {
printf $i "\t"}}
else {
if(a[1] == "./.") {
print a[1]":.:.:.:.\t"}
else {
print $i"\t"}}}}' >> ${VCFout}

}

#Copy and Paste the ENTIRE FUNCTION
#Then just use it with FIRST RGUMENT = INPUT VCF; second ARGUMENT = OUTPUT VCF --> Changing_allFIELDS_from_doubleDots_Geno input_VCF.vcf output_VCF.vcf


Changing_allFIELDS_from_doubleDots_Geno min_cov_5.recode.vcf min_cov_5_clean.recode.vcf




##############################
# SNP sharing (3) 
##############################


#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_3eme/maf_0.005

module load gcc/11.4.0
module load vcftools/0.1.16


#######
### SNPS SHARING - 85%
#######

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.85 --recode --recode-INFO-all --out min_cov_5_clean_85miss

# Resultat
After filtering, kept 359 out of 359 Individuals
Outputting VCF file...
After filtering, kept 66785 out of a possible 107274 Sites
Run Time = 40.00 seconds


## quality check pour SNP sharing 85

module load gcc/13.2.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het


### il y a plus de missing individuals --> verifier si les mêmes individus que filtres avant



#############################
# SNP sharing (3) MAC = 5
##############################

#In /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_3eme/mac_5


#######
### Minor allele frequency filter 
#######


vcftools --vcf min_cov_5_clean_85miss.recode.vcf --mac 5 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_85miss_mac

# Resultat
After filtering, kept 359 out of 359 Individuals
Outputting VCF file...
After filtering, kept 44078 out of a possible 66785 Sites
Run Time = 28.00 seconds


#######
### Filtre gène paralogue
#######


vcftools --vcf min_cov_5_85miss_mac.recode.vcf --min-meanDP 10 --max-meanDP 26 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_wo_hwe


# Resultat
After filtering, kept 359 out of 359 Individuals
Outputting VCF file...
After filtering, kept 38756 out of a possible 44078 Sites
Run Time = 23.00 seconds


#######
### Filtre equilibre de hardy weinberg
#######

vcftools --vcf min_cov_5_clean_85miss_meanDP_wo_hwe.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_hwe

# Resultat
After filtering, kept 359 out of 359 Individuals
Outputting VCF file...
After filtering, kept 35919 out of a possible 38756 Sites
Run Time = 20.00 seconds


####################################
# QUALITY CHECK (3) pour MAC = 5
####################################

## without HWE

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_wo_hwe.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_wo_we

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het


## with HWE

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_hwe.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/final/population_3eme/min_cov_5_clean_85miss_meanDP_hwe


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het


vcftools --gzvcf $VCF --out $OUT --counts --derived


####### TADJIMA

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user alyzee.tassel@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --job-name populations
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

module load gcc/13.2.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_85miss.recode.vcf --TajimaD  1000  --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/missing_85/population_4eme/Tajima_resultat


# Submitted batch job 51163303
# Slurm Job_id=51163303 Name=populations Began, Queued time 00:00:08
# Slurm Job_id=51163303 Name=populations Ended, Run time 00:00:20, COMPLETED, ExitCode 0



module load gcc/13.2.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --TajimaD  1000  --out /work/FAC/FBM/DEE/pchriste/default/atassel/SNP_filtering/population_0.8/population_3eme/Tajima_resultat




