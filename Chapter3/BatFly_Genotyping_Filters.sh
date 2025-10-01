          #################################################
					##### Genomic Pipeline PE - NYCT - ANALYSIS #####
					#################################################
					# De novo Pipeline

## To copy files or folder (-r) from my computer to the cluster
#scp /PATH/TO/MY/FILES ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/FOLDER/FILE
#scp -r /PATH/TO/FOLDER ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/FOLDER/



## To check version
#module spider NAMEOFTHESOFTWARE

## Here there is no reference genome --This is a de novo pipeline


#######
### Download sequences from LIMS ###
#######

## Create a file with all the links from LIMS

wget -i /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/rawdata/download_links_NYCT.txt





#################################################################################################################

			###############################
			##### 1 - Reads filtering #####
			###############################

#################################################################################################################


## The folder ./rawdata is now on the /nas

#######
### Trim the Illumina adapter sequences ###
#######

mkdir /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS	## Create a folder for each step
mkdir /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource	## For the adapters sequences


###
### ./1.FILTERED_READS/Filter_reads.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name NYCT_Adapter_Trimming
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

set -x
set -e

module load gcc/11.4.0
module load bbmap/39.01
module load r/4.3.2

for f in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/rawdata/*_R1_*.fastq.gz
do
        lib=${f/_R1_*.fastq.gz}
        outname=`basename $lib`
        echo "  processing      " $lib " ...."

        R1=`ls /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/rawdata/${outname}_R1_*.fastq.gz`
        R2=`ls /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/rawdata/${outname}_R2_*.fastq.gz`

        echo ${R1}
        echo ${R2}

        mkdir -p /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/

        bbduk.sh in1=${R1} in2=${R2} out1=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R1.bbduk.fastq.gz out2=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R2.bbduk.fastq.gz ref=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/adapters.fa ktrim=r k=20 mink=10 hdist=1 tpe tbo overwrite=true bhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/bhist_${outname}.txt qhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/qhist_${outname}.txt aqhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/aqhist_${outname}.txt lhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/lhist_${outname}.txt >& /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}.bbduk.log

    mv /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R1.bbduk.fastq.gz /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R1.cut.fastq.gz
    mv /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R2.bbduk.fastq.gz /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R2.cut.fastq.gz
    echo "      Done"


echo /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}


#######
###UNZIP_FOR_LENGTH_HIST
#######

gzip -dc /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/${outname}_R1.cut.fastq.gz | awk '(FNR%4)==2{l[length($0)]+=1}END{for (i in l)
print i, l[i];}' > /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/reads.length.txt

LENGTH=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/reads.length.txt


#######
###LENGTH_HISTOGRAMS
#######

R --vanilla <<EOF

  table_compo <- read.table("$LENGTH", header=F, row.names=1)
  tutu <- as.character(row.names(table_compo)[order(as.numeric(row.names(table_compo)))])
  table_compo <- as.data.frame(table_compo[order(as.numeric(row.names(table_compo))),])
  row.names(table_compo)<-tutu

  pdf("/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/Reads_length.pdf", width = 10, height = 7)

  par(mar=c(8, 6, 4, 6))

  b = barplot(table_compo[,1], col = "blue", xlab="Reads length after adapters trimming", names.arg = row.names(table_compo), las = 3, cex.names=0.5, axes = F)
  axis(2, las=2)
  mtext("Nb reads", side = 2, line = 4, cex.lab = 1, las = 3)
  toto <- NULL
  for (i in as.numeric(row.names(table_compo))){
        toto <- c(toto, sum(table_compo[as.numeric(row.names(table_compo))<i,1]))
  }
  tutu <- 1-(toto/sum(table_compo[,1]))

  par(new=TRUE)
  plot(1, type="n", xlim = c(min(b), max(b)), ylim = c(0, 1), axes = F, xlab="", ylab="")
  c = b + 0.5
  points(tutu~(c), pch = 20, type = "p", cex = 0.7, col = "grey70")
  points(tutu~(c), pch = 20, type = "l", cex = 0.7, col = "grey70")
  axis(4, las = 2)
  mtext("Percent of retained reads", side=4, line=3, cex.lab=1, las=3)

  dev.off()

EOF

done



#Submitted batch job 46093078 
#Slurm Job_id=46093078 Name=NYCT_Adapter_Trimming Began, Queued time 00:00:27





#################################################################################################################

			###########################################
			##### 2 - Quality statistics on reads #####
			###########################################

#################################################################################################################



###
### ./1.FILTERED_READS/Quality_nucl.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name NYCT_QualStat
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load r/4.3.2

for forward in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*/*_R1.cut.fastq.gz
   do
	lib=${forward/_R1.cut.fastq.gz/}	## Find all the files ending w/ _R1.cut.fastq.gz
	outname=`basename ${lib}`
	COMPO=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/bhist_${outname}.txt
 	QUAL=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/qhist_${outname}.txt
 

R --vanilla <<EOF


  table_compo <- read.table("$COMPO", header=F, row.names=1)
  prop <- t(as.matrix(table_compo[c(1:(nrow(table_compo)/2)),]))
  table_quality <- read.table("$QUAL", header=F, row.names=1)
  
  pdf(paste0("/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/${outname}/","${outname}","_quality_nucleotide.pdf"), width=10, height=15)
  layout(matrix(c(1:2), ncol = 1))
  plot(table_quality[,1] ~ as.numeric(row.names(table_quality)), xlim=c(0,(nrow(table_quality)+3)), ylim=c(0,(max(table_quality[,1])+10)), col="blue", pch=20, xlab="Nucleotide position along the read", ylab="Quality score", main="${outname}")
  ## Nucleotid proportion at each reads position
  barplot(prop, col=c("blue", "white", "green", "pink", "black"), xlab="nucleotide composition along the read")
  legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))
  dev.off()


EOF

done

#Submitted batch job 46111050
#Slurm Job_id=46111050 Name=NYCT_QualStat Ended, Run time 00:00:20, COMPLETED, ExitCode 0





#################################################################################################################

			######################################
			##### 3 - Library Demultiplexing #####
			######################################

#################################################################################################################



#######
### Preparing demultiplexing ###
#######

## Creating new directories

mkdir PATH/2.DEMULTIPLEXED


#######
### Process RADtags ###
#######

## After Process_radtags_Nyct.R since it was sequenced w/ UDI + EcoRI barcodes

###
### ./2.DEMULTIPLEXED/Process_radtags_NYCT.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name NYCT_demultiplexing
#SBATCH -o %j.stdoutq
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A1*/*A1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A1*/*A1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A10*/*A10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A10*/*A10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A11*/*A11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A11*/*A11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A11/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A12*/*A12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A12*/*A12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A12/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A2*/*A2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A2*/*A2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A3*/*A3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A3*/*A3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A3/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A4*/*A4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A4*/*A4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A4/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A5*/*A5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A5*/*A5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A5/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A6*/*A6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A6*/*A6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A6/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A7*/*A7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A7*/*A7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A7/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A8*/*A8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A8*/*A8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A8/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A9*/*A9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*A9*/*A9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_A9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_A9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B1*/*B1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B1*/*B1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B10*/*B10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B10*/*B10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B11*/*B11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B11*/*B11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B11/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B12*/*B12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B12*/*B12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B12/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B2*/*B2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B2*/*B2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B3*/*B3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B3*/*B3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B3/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B4*/*B4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B4*/*B4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B4/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B5*/*B5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B5*/*B5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B5/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B6*/*B6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B6*/*B6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B6/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B7*/*B7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B7*/*B7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B7/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B8*/*B8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B8*/*B8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B8/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B9*/*B9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*B9*/*B9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_B9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_B9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*C10*/*C10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*C10*/*C10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_C10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_C10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*C9*/*C9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*C9*/*C9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_C9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_C9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*D9*/*D9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*D9*/*D9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_D9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_D9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E1*/*E1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E1*/*E1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E10*/*E10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E10*/*E10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E11*/*E11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E11*/*E11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E11/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E12*/*E12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E12*/*E12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E12/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E2*/*E2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E2*/*E2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E3*/*E3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E3*/*E3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E3/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E4*/*E4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E4*/*E4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E4/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E5*/*E5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E5*/*E5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E5/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E6*/*E6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E6*/*E6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E6/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E7*/*E7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E7*/*E7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E7/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E8*/*E8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E8*/*E8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E8/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E9*/*E9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*E9*/*E9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_E9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_E9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F1*/*F1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F1*/*F1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F10*/*F10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F10*/*F10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F11*/*F11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F11*/*F11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F11/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F12*/*F12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F12*/*F12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F12/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F2*/*F2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F2*/*F2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F3*/*F3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F3*/*F3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F3/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F4*/*F4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F4*/*F4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F4/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F5*/*F5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F5*/*F5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F5/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F6*/*F6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F6*/*F6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F6/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F7*/*F7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F7*/*F7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F7/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F8*/*F8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F8*/*F8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F8/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F9*/*F9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*F9*/*F9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_F9/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_F9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*G10*/*G10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*G10*/*G10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_G10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_G10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*G12*/*G12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*G12*/*G12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_G12/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_G12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*H10*/*H10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/1.FILTERED_READS/*H10*/*H10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_H10/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/ressource/barcode_NYCT/barcodes_NYCT_H10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90


#Submitted batch job 47274446
#Slurm Job_id=47274446 Name=NYCT_demultiplexing Began, Queued time 00:00:28
#Slurm Job_id=47274446 Name=NYCT_demultiplexing Ended, Run time 04:49:32, COMPLETED, ExitCode 0



#######
### Nb of reads/samples ###
#######

## Copy the files in another folder in order to copy them on my computer

mkdir To_copy

for i in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/NYCT_*/process_radtags*.log
do
    outname=$(basename $(dirname $i))
    cp $i /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/To_copy/process_radtags.${outname}.log
done


## Copy in computer
scp -r ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED/To_copy C:/Users/ccastex/Documents/5.Thesis/Project/7.RADseq/Project_RAD2024/Parasites/Nyct/V2_Results/Process_radtags





# ## In R do :
# TAB=read.csv("Process_radtags_retained_reads.csv", header=T, sep=";")
# TAB_Order = TAB[order(TAB[,6]),]
# TAB_Reduced = rbind(TAB_Order[,6] , TAB_Order[,3]-TAB_Order[,6])
# pdf("./Reads_filtration.pdf", width = 30, height = 15)
# layout(matrix(c(1,1,1,2),ncol = 2), widths=c(2,1))
# x = barplot(height = TAB_Reduced, col = c("blue", "red"), main = paste0("Number of reads for each individual"), las=2, cex.names = 0.4, space=0.2, border=NA)
# text(cex=0.8, x=x, y=-max(TAB_Order[,4])*0.02, TAB_Order[,2], xpd=TRUE, srt=90, adj=1)
# legend(x = "topright", legend = c("Total", "Retained"), fill = c("red","blue"), bty = "n")
# plot.new()
# pie(rowSums(TAB_Reduced), col = c("blue","red"), radius = 1, init.angle = 90, labels = round(c(rowSums(TAB_Reduced)/sum(rowSums(TAB_Reduced)))*100, digits = 2), main = "Mean prop.")
# dev.off()



### Move rem files

mkdir rem_files

base_dir="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED"

for subdir in "$base_dir"/NYCT_*/; do
    for file in "${subdir}"*.rem.*.fq.gz; do
        if [[ -e "$file" ]]; then
            mv "$file" "$base_dir/rem_files/"
        fi
    done
done


#######
### Merge paired-end ###
#######

###
### ./3.MERGED_PE/concatenate_fastq.py
###

#!/usr/bin/env python3

import argparse
import gzip

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence using basic Python."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def concatenate_fastq(forward_fastq, reverse_fastq, output_fastq):
    with gzip.open(forward_fastq, 'rt') as forward, gzip.open(reverse_fastq, 'rt') as reverse, gzip.open(output_fastq, 'wt') as output:
        while True:
            # Read four lines from the forward and reverse FASTQ files
            forward_id = forward.readline().strip()
            forward_seq = forward.readline().strip()
            forward_plus = forward.readline().strip()
            forward_quality = forward.readline().strip()

            reverse_id = reverse.readline().strip()
            reverse_seq = reverse.readline().strip()
            reverse_plus = reverse.readline().strip()
            reverse_quality = reverse.readline().strip()

            # Check if we've reached the end of the file
            if not forward_id or not reverse_id:
                break

            # Reverse complement the reverse read sequence
            reverse_seq_rc = reverse_complement(reverse_seq)

            # Concatenate sequences and quality scores
            concatenated_seq = forward_seq + reverse_seq_rc
            concatenated_quality = forward_quality + reverse_quality

            # Write the concatenated reads in FASTQ format
            output.write(f"{forward_id}\n{concatenated_seq}\n{forward_plus}\n{concatenated_quality}\n")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Concatenate paired-end FASTQ files with reverse complement for reverse reads.")
    parser.add_argument("-F", "--forward", required=True, help="Forward FASTQ file (R1)")
    parser.add_argument("-R", "--reverse", required=True, help="Reverse FASTQ file (R2)")
    parser.add_argument("-O", "--output", required=True, help="Output FASTQ file")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to concatenate the FASTQ files
    concatenate_fastq(args.forward, args.reverse, args.output)

if __name__ == "__main__":
    main()


## Then to lauch this file :

###
### ./3.MERGED_PE/python_concatenate.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name NYCT_concatenate_fastq
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load python/3.11.6

base_dir="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/2.DEMULTIPLEXED"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/3.MERGED_PE"

for sub_dir in ${base_dir}/NYCT_*
do
  if [[ -d "${sub_dir}" ]]; then
    for forward in ${sub_dir}/*.1.fq.gz; do
      sample=${forward/.1.fq.gz}
      sample_id=$(basename "$sample")
      echo "Processing ..." ${sample_id}
      if [[ -f "${forward}" && -s "${forward}" ]]; then
        reverse="${forward/.1.fq.gz/.2.fq.gz}"
        if [[ -f "${reverse}" && -s "${reverse}" ]]; then
          echo "Processing..." ${forward} " and " ${reverse}
          python concatenate_fastq.py -F ${forward} -R ${reverse} -O ${OUTPATH}/${sample_id}.merged.fq.gz
          echo "Merged sequences of sample ${sample_id}"
        else
          echo "Missing of empty reverse files for sample ${sample_id}, skipping..."
        fi
      else
        echo "Missing or empty forward file for sample ${sample_id}, skipping..."
      fi
    done
  else
    echo "Missing or empty ${sub_dir}"
  fi
done

#Submitted batch job 47329113
#Slurm Job_id=47329113 Name=NYCT_concatenate_fastq Began, Queued time 00:00:01
#Slurm Job_id=47329113 Name=NYCT_concatenate_fastq Ended, Run time 06:39:39, COMPLETED, ExitCode 0



#################################################################################################################

      #########################################
      ##### 4 - SNP Assembly For Analysis #####
      #########################################

#################################################################################################################



##############################################################################################################
##############################################################################################################

#####################
##### ON CONCAT #####
#####################

mkdir 4.RUN_STACKS/ANALYSIS
# The analysis for the NYCT will happen in this directory


#######
### USTACKS - stacks assembly ### ## ON CONCAT + USING DEFAULT PARAMETERS
#######

mkdir M3m4

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/ustacks_analysis_concat_M3m4.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name NYCT_ustacks_CONCAT_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4"
base_dir="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/3.MERGED_PE"

rm -f ${OUTPATH}/Stat_ustacks_M4m5.log
rm -f ${OUTPATH}/ustacks_err_M4m5
rm -f ${OUTPATH}/cov_M4m5
rm -f ${OUTPATH}/nbstacks_M4m5

echo Ind nb_unique_stacks nb_merged_stacks mean_cov > ${OUTPATH}/Stat_ustacks_M3m4.log

i="0"

for merged_sample in ${base_dir}/*.merged.fq.gz; do
  echo `basename ${merged_sample}`
  ustacks -p 8 -f ${merged_sample} -i ${i} -o ${OUTPATH} -t gzfastq -M 3 -m 4 2> ${OUTPATH}/ustacks_err_M3m4
  i=$[${i}+1]
  awk '/Assembled/' ${OUTPATH}/ustacks_err_M3m4 | awk '{print $2,$5}' | tr -d "\n" | awk '{print $1,$3}' > ${OUTPATH}/nbstacks_M3m4
  awk '/Final coverage/ {print $3}' ${OUTPATH}/ustacks_err_M3m4 | sed -e 's/mean=//g' | sed -e 's/;//g'  > ${OUTPATH}/cov_M3m4
  chem=${merged_sample/.R1.fq.gz/}
  name=${chem/${base_dir}\//}
  echo ${name} `cat ${OUTPATH}/nbstacks_M3m4` `cat ${OUTPATH}/cov_M3m4` >> ${OUTPATH}/Stat_ustacks_M3m4.log
done


#Submitted batch job 48402885
#Slurm Job_id=48402885 Name=NYCT_ustacks_CONCAT_M3m4 Began, Queued time 00:00:10
#Slurm Job_id=48402885 Name=NYCT_ustacks_CONCAT_M3m4 Ended, Run time 08:41:45, COMPLETED, ExitCode 0


## Mean cov = 21.99
## Mean nb unique stacks = 47774.625
## Mean nb merged stacks = 30856.307


# ## In R do :

# TABLE = read.csv("Stat_ustacks_CONCAT.csv", header = T, sep=";")
# str(TABLE)
# TABLE = TABLE[order(TABLE[,2]),]
# names=unlist(lapply(strsplit(as.character(TABLE[,1]), "./", fixed=TRUE), function(x) x[2]))
# pdf("Indiv_FragCover.pdf", width=(5+round(nrow(TABLE)/5)), height=10)
# layout(matrix(c(1,1,1,2),ncol = 2), widths=c(2,1))
# b = barplot(TABLE[,2], col = "darkblue", names.arg = names, las = 3, cex.names = 0.7, axes = F)
# axis(2,las = 2)
# barplot(TABLE[,3], col = "lightblue", add = T, axes = F)
# mtext("Nb Stacks", side=2, line=3, cex.lab=1,las=3, col="lightblue")
# legend(x = "topright", legend = c("Nb unique stacks", "Nb merged stacks"), fill = c("darkblue","lightblue"), bty = "n")
# par(new=TRUE)
# plot(1,type = "n", xlim = c(min(b), max(b)+1), ylim = c(0, max(TABLE[,4])), axes = F, xlab ="", ylab="")
# c = b+0.5
# points(TABLE[,4]~(c), pch = 20, type = "b")
# axis(4,las = 2)
# mtext("Median cov", side=4, line=3, cex.lab=1, las=3)
# dev.off()




#######
### CSTACKS - build catalog ### ## ON CONCAT
#######


## Here : -n 4 (nb of differences within an individual + 1)

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/cstacks_analysis_concat_n4.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 25
#SBATCH --mem 100G
#SBATCH --job-name NYCT_cstacks_AnalysisCONCAT_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4"
TEMP_DIR="${OUTPATH}/TEMP"

mkdir -p ${TEMP_DIR}


LIST_FILE="${TEMP_DIR}/list"

if [ -e  "${LIST_FILE}" ]
then 
  rm ${LIST_FILE}
fi

for forward in ${OUTPATH}/*.snps.tsv.gz
    do
        base_name=$(basename ${forward} .snps.tsv.gz)
        echo "-s ${OUTPATH}/${base_name}" >> ${LIST_FILE}
    done

cstacks -o ${OUTPATH} -p 25 -n 4 $(cat ${LIST_FILE})


#Submitted batch job 48452864
#Slurm Job_id=48452864 Name=NYCT_cstacks_AnalysisCONCAT_M3m4 Began, Queued time 00:00:21
#Slurm Job_id=48452864 Name=NYCT_cstacks_AnalysisCONCAT_M3m4 Ended, Run time 1-18:41:24, COMPLETED, ExitCode 0

Final catalog contains 506072 loci.


#######
### SSTACKS - match catalog ### ## ON CONCAT
#######

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/sstacks_concat_Analysis_M3m4.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 25
#SBATCH --mem 100G
#SBATCH --job-name NYCT_sstacksAnalysis_CONCAT_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4"

for forward in ${OUTPATH}/*.snps.tsv.gz
do
    if [[ ! $forward =~ catalog ]]
    then
        fwd=${forward/.snps.tsv.gz/}
        sstacks -p 25 -c ${OUTPATH} -s ${fwd} -o ${OUTPATH}/
    fi
done

#Submitted batch job 48573599
#Slurm Job_id=48573599 Name=NYCT_sstacksAnalysis_CONCAT_M3m4 Began, Queued time 00:00:28
#Slurm Job_id=48573599 Name=NYCT_sstacksAnalysis_CONCAT_M3m4 Ended, Run time 1-21:33:20, COMPLETED, ExitCode 0



#######
### TSV2BAM - convert files ### ## ON CONCAT
#######

mkdir ./M3m4/CATALOG

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/tsv2bam_concat_Analysis_M3m4.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 72G
#SBATCH --job-name NYCT_tsv2bam_CONCATAnalysis_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/CATALOG/NYCT_popmapAnalysis_concat.txt"

tsv2bam -P ${INPATH} -M ${POP_INFOS} -t 16
mv "${INPATH}"/*.bam "${OUTPATH}"/

#Submitted batch job 48616713
#Slurm Job_id=48616713 Name=NYCT_tsv2bam_CONCATAnalysis_M3m4 Began, Queued time 00:00:12
#Slurm Job_id=48616713 Name=NYCT_tsv2bam_CONCATAnalysis_M3m4 Ended, Run time 00:06:22, COMPLETED, ExitCode 0



#######
### GSTACKS - genotyping ###  ## ON CONCAT
#######

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/gstacks_concat_Analysis_M3m4.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 72G
#SBATCH --job-name NYCT_gstacksAnalysis_CONCAT_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53


OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/CATALOG/NYCT_popmapAnalysis_concat.txt"

gstacks -P ${OUTPATH} -M ${POP_INFOS} -t 16

#Submitted batch job 48618858
#Slurm Job_id=48618858 Name=NYCT_gstacksAnalysis_CONCAT_M3m4 Began, Queued time 00:00:29
#Slurm Job_id=48618858 Name=NYCT_gstacksAnalysis_CONCAT_M3m4 Ended, Run time 00:24:20, COMPLETED, ExitCode 0


## RESULTS:

Genotyped 437712 loci:
  effective per-sample coverage: mean=27.3x, stdev=13.5x, min=7.0x, max=107.0x
  mean number of sites per locus: 180.5
  a consistent phasing was found for 1507139 of out 1797082 (83.9%) diploid loci needing phasing



#######
### POPULATIONS - SNPs calling ###  ## ON CONCAT
#######

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/populations_concat_Analysis_M3m4_p14.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 72G
#SBATCH --job-name NYCT_populations_CONCATAnalysis_M3m4
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53


INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/CATALOG/NYCT_popmapAnalysis_concat.txt"

populations -P ${INPATH} --popmap ${POP_INFOS} -O ${OUTPATH} -p 14 -r 0.6 -R 0.6 -f p_value -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp



#### Test with -p 14

#Submitted batch job 48621795
#Slurm Job_id=48621795 Name=NYCT_populations_CONCATAnalysis_M3m4 Began, Queued time 03:24:24
#Slurm Job_id=48621795 Name=NYCT_populations_CONCATAnalysis_M3m4 Ended, Run time 00:16:30, COMPLETED, ExitCode 0


# RESULTS:
Removed 426927 loci that did not pass sample/population constraints from 437712 loci.
Kept 10785 loci, composed of 1955448 sites; 6072 of those sites were filtered, 10768 variant sites remained.
Mean genotyped sites per locus: 181.31bp (stderr 0.04).


#### Test with -p 10
#### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/p10/populations_concat_Analysis_M3m4_p10.sh
#### #SBATCH --job-name NYCT_populations_CONCATAnalysis_M3m4_p10
#### OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10/"


#Submitted batch job 48621790
#Slurm Job_id=48621790 Name=NYCT_populations_CONCATAnalysis_M3m4_p10 Began, Queued time 03:09:01
#Slurm Job_id=48621790 Name=NYCT_populations_CONCATAnalysis_M3m4_p10 Ended, Run time 00:17:16, COMPLETED, ExitCode 0


# RESULTS:
Removed 424887 loci that did not pass sample/population constraints from 437712 loci.
Kept 12825 loci, composed of 2325117 sites; 6576 of those sites were filtered, 12795 variant sites remained.
Mean genotyped sites per locus: 181.30bp (stderr 0.03).


## KEEP -p 10



#################################################################################################################

      ################################################
      ##### 5 - SNP Filtering For Analysis - p10 #####
      ################################################

#################################################################################################################

mkdir ./5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4


#######
### MIN COVERAGE
#######

module load gcc/11.4.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out ${OUTPATH}/min_cov_5

# OUTPATH:
After filtering, kept 12795 out of a possible 12795 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

Changing_allFIELDS_from_doubleDots_Geno(){

#First argument passed to the script is input VCF
VCFin=$1
#Second argument passed to the script is what you want as output
VCFout=$2

#First step is to extract header
grep '#' ${VCFin} > ${VCFout}

#Get the rest and take the first columns (not INDV/GENO, unvariable (such as CHR, POS, etc.)
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

#Then just use it with FIRST RGUMENT = INPUT VCF; second ARGUMENT = OUTPUT VCF -->

#Changing_allFIELDS_from_doubleDots_Geno input_VCF.vcf output_VCF.vcf

### Code ###

Changing_allFIELDS_from_doubleDots_Geno min_cov_5.recode.vcf min_cov_5_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out min_cov_5_clean_70miss

# OUTPUT:
After filtering, kept 11248 out of a possible 12795 Sites


#######
### SNPS SHARING - 80%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out min_cov_5_clean_80miss

# OUTPUT:
After filtering, kept 8958 out of a possible 12795 Sites


#######
### SNPS SHARING - 90%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out min_cov_5_clean_90miss

# OUTPUT:
After filtering, kept 5142 out of a possible 12795 Sites



###############
###############
#### ALL THE NEXT ANALYSES ARE DONE ON THE 80% SHARED SNPS BUT I WILL ALSO TRY 70% 
###############
###############


#######
### Test maf 0.01 sur 80% shared snps
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

mkdir Filter_maf

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss.recode.vcf --maf 0.01 --recode --recode-INFO-all --out ./Filter_maf/min_cov_5_clean_80miss_maf

# OUTPUT:
After filtering, kept 2018 out of a possible 8958 Sites


#######
### Test mac 5 sur 80% shared snps
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4

mkdir Filter_mac

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss.recode.vcf --mac 5 --recode --recode-INFO-all --out ./Filter_mac/min_cov_5_clean_80miss_mac

# OUTPUT:
After filtering, kept 2719 out of a possible 8958 Sites



#######
### QUALITY CHECK WITHOUT MAF OR MAC
#######

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/min_cov_5_clean_80miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/min_cov_5_clean_80miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het



# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/min_cov_5_clean_80miss* ./M3m4/nomac/

### Mean depth site = 32.09X
### Mean depth individuals = 31.17X
### Nb individuals >25% missing data = 30
### ID = 401N2 ; 419N2 ; 428N1 ; 428N2 ; 431N2 ; 436N2 ; 442N1 ; 443N2 ; 449N2 ; 451N1 ;
###      469N1 ; 479N1 ; 484N1 ; 487N2 ; 616N2 ; 635N1 ; 646N1 ; 646N2 ; 651N1 ; 674N1 ;
###      676N2 ; 679N1 ; 690N2 ; 692N1 ; 697N1 ; 708N1 ; 719N1 ; 719N2 ; 770N2 ; 792N2



##############################################################################################################

## Redo from populations on the individuals with <25% missing data

# New population map in : /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10/Ind25miss



#######
### POPULATIONS - SNPs calling ###  ## ON CONCAT
#######

###
### ./4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/p10/Ind25miss/populations_concat_Analysis_M3m4_p10_25miss.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 72G
#SBATCH --job-name NYCT_populations_CONCATAnalysis_M3m4_2
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/11.4.0
module load stacks/2.53


INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10/Ind25miss"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10/Ind25miss/NYCT_popmapAnalysis_concat25miss.txt"

populations -P ${INPATH} --popmap ${POP_INFOS} -O ${OUTPATH} -p 10 -r 0.6 -R 0.6 -f p_value -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp

#Submitted batch job 48647061
#Slurm Job_id=48647061 Name=NYCT_populations_CONCATAnalysis_M3m4_2 Began, Queued time 00:00:18
#Slurm Job_id=48647061 Name=NYCT_populations_CONCATAnalysis_M3m4_2 Ended, Run time 00:15:10, COMPLETED, ExitCode 0

Removed 424069 loci that did not pass sample/population constraints from 437712 loci.
Kept 13643 loci, composed of 2473307 sites; 6726 of those sites were filtered, 13596 variant sites remained.
Mean genotyped sites per locus: 181.29bp (stderr 0.03).




mkdir ./5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

#######
### MIN COVERAGE
#######

module load gcc/11.4.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/4.RUN_STACKS/ANALYSIS/CONCAT/M3m4/CATALOG/POPULATIONS/p10/Ind25miss"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out ${OUTPATH}/min_cov_5

# OUTPUT:
After filtering, kept 13596 out of a possible 13596 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

Changing_allFIELDS_from_doubleDots_Geno(){

#First argument passed to the script is input VCF
VCFin=$1
#Second argument passed to the script is what you want as output
VCFout=$2

#First step is to extract header
grep '#' ${VCFin} > ${VCFout}

#Get the rest and take the first columns (not INDV/GENO, unvariable (such as CHR, POS, etc.)
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

#Then just use it with FIRST RGUMENT = INPUT VCF; second ARGUMENT = OUTPUT VCF -->

#Changing_allFIELDS_from_doubleDots_Geno input_VCF.vcf output_VCF.vcf

### Code ###

Changing_allFIELDS_from_doubleDots_Geno min_cov_5.recode.vcf min_cov_5_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out min_cov_5_clean_70miss

# OUTPUT:
After filtering, kept 12088 out of a possible 13596 Sites


#######
### SNPS SHARING - 80%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out min_cov_5_clean_80miss

# OUTPUT:
After filtering, kept 10016 out of a possible 13596 Sites


#######
### SNPS SHARING - 90%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out min_cov_5_clean_90miss

# OUTPUT:
After filtering, kept 6623 out of a possible 13596 Sites



#######
### MINIMUM COV 10 / MAXIMUM COV 65
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss.recode.vcf --min-meanDP 10 --max-meanDP 65 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65

# OUTPUT :
After filtering, kept 9900 out of a possible 10016 Sites


#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65_hweq

# OUTPUTS :
After filtering, kept 7597 out of a possible 9900 Sites



#######
### QUALITY CHECK WITHOUT MAF OR MAC BUT WITH HWE
#######

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq* ./


#######
### MAF - 0.01
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65_hweq.recode.vcf --maf 0.01 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65_hweq_maf

# OUTPUTS :
After filtering, kept 1122 out of a possible 7597 Sites


#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65_hweq_mac

# OUTPUTS :
After filtering, kept 1544 out of a possible 7597 Sites


#######
### MAC - 5 - HETEROZYGOSITY & DEPTH
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/Ind28miss

module load gcc/13.2.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_hweq_mac* ./


#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65_hweq_mac.recode.vcf --TajimaD 1000 --out mac_tajima_results


## To do in R 
# Check significativity


#######
### QUALITY CHECK WITHOUT HWE
#######

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het
vcftools --gzvcf $VCF --out $OUT --depth

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65* ./


## To do in R 
# HWE graph 
# Het/cov graph
# Maf graph 



##############################################################################################################

# Do the analysis on the SNPs with no HWE


#######
### QUALITY CHECK WITHOUT HWE
#######

#cf upper for the quality check


#######
### MAF - 0.01
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65.recode.vcf --maf 0.01 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65_maf

# OUTPUTS :

After filtering, kept 2303 out of a possible 9900 Sites


#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65.recode.vcf --mac 5 --recode --recode-INFO-all --out min_cov_5_clean_80miss_cov10_65_mac

# OUTPUTS :

After filtering, kept 3061 out of a possible 9900 Sites


#######
### NO HWE - MAC 5 - HETEROZYGOSITY & DEPTH
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load samtools/1.17
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/min_cov_5_clean_80miss_cov10_65_mac* ./


#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/11.4.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65_mac.recode.vcf --TajimaD 1000 --out mac_nohwe_tajima_results

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/mac_nohwe_tajima_results.Tajima.D ./



########################
##
## Tajima's D has to be done without mac filter!!!
##

#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss

module load gcc/13.2.0
module load vcftools/0.1.16

vcftools --vcf min_cov_5_clean_80miss_cov10_65_hweq.recode.vcf --TajimaD 1000 --out hwe_nomac_tajima_results

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/NYCT/5.SNPS_FILTERING/ANALYSIS/CONCAT/M3m4/Ind25miss/hwe_nomac_tajima_results.Tajima.D ./
