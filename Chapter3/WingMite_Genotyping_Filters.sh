					##################################################################
					##### Genomic Pipeline PE - SPX - CONCAT ONLY - RESEQUENCING #####
					##################################################################
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

wget -i /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/rawdata/download_links_resequencing2025_SPX.txt





#################################################################################################################

			###############################
			##### 1 - Reads filtering #####
			###############################

#################################################################################################################



#######
### Trim the Illumina adapter sequences ###
#######

mkdir /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025	## Create a folder for each step
mkdir /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource	## For the adapters sequences


###
### ./1.FILTERED_READS/filter_reads2025.sh
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
#SBATCH --job-name SPX_Adapter_Trimming2025
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

set -x
set -e

module load gcc/12.3.0
module load bbmap/39.01
module load r-light/4.4.1

for f in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/rawdata/*_R1_*.fastq.gz
do
        lib=${f/_R1_*.fastq.gz}
        outname=`basename $lib`
        echo "  processing      " $lib " ...."

        R1=`ls /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/rawdata/${outname}_R1_*.fastq.gz`
        R2=`ls /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/rawdata/${outname}_R2_*.fastq.gz`

        echo ${R1}
        echo ${R2}

        mkdir -p /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/

        bbduk.sh in1=${R1} in2=${R2} out1=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R1.bbduk.fastq.gz out2=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R2.bbduk.fastq.gz ref=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/adapters.fa ktrim=r k=20 mink=10 hdist=1 tpe tbo overwrite=true bhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/bhist_${outname}.txt qhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/qhist_${outname}.txt aqhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/aqhist_${outname}.txt lhist=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/lhist_${outname}.txt >& /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}.bbduk.log

    mv /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R1.bbduk.fastq.gz /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R1.cut.fastq.gz
    mv /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R2.bbduk.fastq.gz /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R2.cut.fastq.gz
    echo "      Done"


echo /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}


#######
###UNZIP_FOR_LENGTH_HIST
#######

gzip -dc /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/${outname}_R1.cut.fastq.gz | awk '(FNR%4)==2{l[length($0)]+=1}END{for (i in l)
print i, l[i];}' > /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/reads.length.txt

LENGTH=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/reads.length.txt


#######
###LENGTH_HISTOGRAMS
#######

R --vanilla <<EOF

  table_compo <- read.table("$LENGTH", header=F, row.names=1)
  tutu <- as.character(row.names(table_compo)[order(as.numeric(row.names(table_compo)))])
  table_compo <- as.data.frame(table_compo[order(as.numeric(row.names(table_compo))),])
  row.names(table_compo)<-tutu

  pdf("/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/Reads_length.pdf", width = 10, height = 7)

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

#Submitted batch job 50480219
#Slurm Job_id=50480219 Name=SPX_Adapter_Trimming2025 Began, Queued time 00:00:27
#Slurm Job_id=50480219 Name=SPX_Adapter_Trimming2025 Ended, Run time 01:28:59, COMPLETED, ExitCode 0





#################################################################################################################

			###########################################
			##### 2 - Quality statistics on reads #####
			###########################################

#################################################################################################################



###
### ./RESEQUENCING/1.FILTERED_READS_2025/Quality_nucl_2025.sh
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
#SBATCH --job-name SPX_QualStat_2025
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load r-light/4.4.1

for forward in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*/*_R1.cut.fastq.gz
   do
	lib=${forward/_R1.cut.fastq.gz/}	## Find all the files ending w/ _R1.cut.fastq.gz
	outname=`basename ${lib}`
	COMPO=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/bhist_${outname}.txt
 	QUAL=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/qhist_${outname}.txt
 

R --vanilla <<EOF


  table_compo <- read.table("$COMPO", header=F, row.names=1)
  prop <- t(as.matrix(table_compo[c(1:(nrow(table_compo)/2)),]))
  table_quality <- read.table("$QUAL", header=F, row.names=1)
  
  pdf(paste0("/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_${outname}/","${outname}","_quality_nucleotide.pdf"), width=10, height=15)
  layout(matrix(c(1:2), ncol = 1))
  plot(table_quality[,1] ~ as.numeric(row.names(table_quality)), xlim=c(0,(nrow(table_quality)+3)), ylim=c(0,(max(table_quality[,1])+10)), col="blue", pch=20, xlab="Nucleotide position along the read", ylab="Quality score", main="${outname}")
  ## Nucleotid proportion at each reads position
  barplot(prop, col=c("blue", "white", "green", "pink", "black"), xlab="nucleotide composition along the read")
  legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))
  dev.off()


EOF

done

#Submitted batch job 50490343
#Slurm Job_id=50490343 Name=SPX_QualStat_2025 Began, Queued time 00:00:01
#Slurm Job_id=50490343 Name=SPX_QualStat_2025 Ended, Run time 00:01:07, COMPLETED, ExitCode 0



###
### ./RESEQUENCING/1.FILTERED_READS_2025/Mean_quality_nucl_2025.sh
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
#SBATCH --job-name SPX_MeanQualStat_2025
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load r-light/4.4.1

# Define paths
ROOT_DIR="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025"
OUTPUT_PDF="${ROOT_DIR}/Mean_Quality_Composition.pdf"

# Temporary files for data (without headers)
ALL_COMPO="${ROOT_DIR}/all_compo.txt"
ALL_QUAL="${ROOT_DIR}/all_qual.txt"

# Clear files if they exist
> "$ALL_COMPO"
> "$ALL_QUAL"

first_compo=1
first_qual=1

# Loop through all folders and collect data
for forward in ${ROOT_DIR}/*/*_R1.cut.fastq.gz
do
    lib=${forward/_R1.cut.fastq.gz/}
    outname=$(basename "${lib}")

    COMPO="${ROOT_DIR}/2025_${outname}/bhist_${outname}.txt"
    QUAL="${ROOT_DIR}/2025_${outname}/qhist_${outname}.txt"

    # Process COMPO file (skip header after the first file)
    if [[ -f "$COMPO" ]]; then
        if [[ $first_compo -eq 1 ]]; then
            cat "$COMPO" > "$ALL_COMPO"
            first_compo=0
        else
            tail -n +2 "$COMPO" >> "$ALL_COMPO"
        fi
    fi

    # Process QUAL file (skip header after the first file)
    if [[ -f "$QUAL" ]]; then
        if [[ $first_qual -eq 1 ]]; then
            cat "$QUAL" > "$ALL_QUAL"
            first_qual=0
        else
            tail -n +2 "$QUAL" >> "$ALL_QUAL"
        fi
    fi
done

# Run R to compute means and create a single PDF
R --vanilla <<EOF
# Read and process composition table
compo_data <- read.table("$ALL_COMPO", header=TRUE, sep="\t")
compo_matrix <- t(as.matrix(compo_data))
mean_compo <- colMeans(compo_matrix, na.rm=TRUE)

# Read and process quality table
qual_data <- read.table("$ALL_QUAL", header=TRUE, sep="\t")
mean_qual <- rowMeans(qual_data, na.rm=TRUE)

# Generate single PDF
pdf("$OUTPUT_PDF", width=10, height=15)
layout(matrix(c(1:2), ncol = 1))

# Quality plot
plot(mean_qual ~ as.numeric(names(mean_qual)), 
     xlab="Nucleotide position along the read", 
     ylab="Mean Quality Score", 
     main="Mean Quality Score across all folders", 
     col="blue", pch=20)

# Composition plot
barplot(mean_compo, col=c("blue", "white", "green", "pink", "black"),
        xlab="Nucleotide composition along the read")
legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))

dev.off()
EOF

echo "PDF created at: $OUTPUT_PDF"


#Submitted batch job 50498907
#Slurm Job_id=50498907 Name=SPX_MeanQualStat_2025 Began, Queued time 00:00:10
#Slurm Job_id=50498907 Name=SPX_MeanQualStat_2025 Ended, Run time 00:00:03, COMPLETED, ExitCode 0


## I did the same for the 1st batch of sequencing
#Submitted batch job 50500385
#Slurm Job_id=50500385 Name=SPX_MeanQualStat Began, Queued time 00:00:23
#Slurm Job_id=50500385 Name=SPX_MeanQualStat Ended, Run time 00:00:04, COMPLETED, ExitCode 0



### Copy all the pdf into my computer
mkdir To_copy

for i in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_SPX*/*_quality_nucleotide.pdf
do
    outname=$(basename $(dirname $i))
    cp $i /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/To_copy/${outname}_quality_nucleotide.pdf
done

for i in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_SPX*/Reads_length.pdf
do
    outname=$(basename $(dirname $i))
    cp $i /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/To_copy/${outname}_Reads_length.pdf
done


## Copy in computer
scp -r ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/To_copy ./
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/all_compo.txt ./
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/all_qual.txt ./
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/1.FILTERED_READS/all_compo.txt ./
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/1.FILTERED_READS/all_qual.txt ./







#################################################################################################################

			######################################
			##### 3 - Library Demultiplexing #####
			######################################

#################################################################################################################



#######
### Preparing demultiplexing ###
#######

## Creating new directories

mkdir PATH/RESEQUENCING/2.DEMULTIPLEXED_2025

## Copy barcodes files into the cluster



#######
### Process RADtags ###
#######

## After Process_radtags_Spx.R since it was sequenced w/ UDI + EcoRI barcodes

## Create new directories

for dir in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/2025_SPX_*_Pool*; do
    if [ -d "$dir" ]; then
        new_name=$(basename "$dir" | sed 's/_L$//')
        mkdir /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/"$new_name"
    fi
done

###
### ./RESEQUENCING/2.DEMULTIPLEXED_2025/Process_radtags_SPX_2025.sh
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
#SBATCH --job-name SPX_demultiplexing_2025
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53

process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*A9_Pool_L1/*A9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*A9_Pool_L1/*A9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_A9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_A9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B11_Pool_L1/*B11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B11_Pool_L1/*B11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_B11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_B11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B9_Pool_L1/*B9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B9_Pool_L1/*B9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_B9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_B9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C1_Pool_L1/*C1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C1_Pool_L1/*C1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C1_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C10_Pool_L1/*C10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C10_Pool_L1/*C10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C11_Pool_L1/*C11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C11_Pool_L1/*C11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C12_Pool_L1/*C12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C12_Pool_L1/*C12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C12_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C2_Pool_L1/*C2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C2_Pool_L1/*C2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C2_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C3_Pool_L1/*C3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C3_Pool_L1/*C3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C3_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C4_Pool_L1/*C4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C4_Pool_L1/*C4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C4_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C5_Pool_L1/*C5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C5_Pool_L1/*C5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C5_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C6_Pool_L1/*C6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C6_Pool_L1/*C6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C6_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C7_Pool_L1/*C7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C7_Pool_L1/*C7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C7_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C8_Pool_L1/*C8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C8_Pool_L1/*C8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C8_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C9_Pool_L1/*C9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C9_Pool_L1/*C9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D1_Pool_L1/*D1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D1_Pool_L1/*D1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D1_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D10_Pool_L1/*D10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D10_Pool_L1/*D10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D11_Pool_L1/*D11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D11_Pool_L1/*D11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D12_Pool_L1/*D12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D12_Pool_L1/*D12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D12_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D2_Pool_L1/*D2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D2_Pool_L1/*D2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D2_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D3_Pool_L1/*D3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D3_Pool_L1/*D3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D3_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D4_Pool_L1/*D4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D4_Pool_L1/*D4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D4_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D5_Pool_L1/*D5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D5_Pool_L1/*D5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D5_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D6_Pool_L1/*D6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D6_Pool_L1/*D6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D6_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D7_Pool_L1/*D7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D7_Pool_L1/*D7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D7_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D8_Pool_L1/*D8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D8_Pool_L1/*D8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D8_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D9_Pool_L1/*D9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D9_Pool_L1/*D9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E10_Pool_L1/*E10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E10_Pool_L1/*E10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_E10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_E10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E12_Pool_L1/*E12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E12_Pool_L1/*E12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_E12_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_E12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F10_Pool_L1/*F10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F10_Pool_L1/*F10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F11_Pool_L1/*F11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F11_Pool_L1/*F11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F12_Pool_L1/*F12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F12_Pool_L1/*F12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F12_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F9_Pool_L1/*F9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F9_Pool_L1/*F9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G1_Pool_L1/*G1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G1_Pool_L1/*G1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G1_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G10_Pool_L1/*G10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G10_Pool_L1/*G10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G11_Pool_L1/*G11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G11_Pool_L1/*G11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G2_Pool_L1/*G2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G2_Pool_L1/*G2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G2_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G3_Pool_L1/*G3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G3_Pool_L1/*G3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G3_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G4_Pool_L1/*G4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G4_Pool_L1/*G4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G4_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G5_Pool_L1/*G5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G5_Pool_L1/*G5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G5_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G6_Pool_L1/*G6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G6_Pool_L1/*G6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G6_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G7_Pool_L1/*G7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G7_Pool_L1/*G7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G7_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G8_Pool_L1/*G8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G8_Pool_L1/*G8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G8_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G9_Pool_L1/*G9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G9_Pool_L1/*G9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H1_Pool_L1/*H1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H1_Pool_L1/*H1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H1_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H10_Pool_L1/*H10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H10_Pool_L1/*H10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H10_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H11_Pool_L1/*H11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H11_Pool_L1/*H11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H11_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H12_Pool_L1/*H12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H12_Pool_L1/*H12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H12_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H2_Pool_L1/*H2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H2_Pool_L1/*H2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H2_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H3_Pool_L1/*H3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H3_Pool_L1/*H3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H3_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H4_Pool_L1/*H4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H4_Pool_L1/*H4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H4_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H5_Pool_L1/*H5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H5_Pool_L1/*H5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H5_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H6_Pool_L1/*H6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H6_Pool_L1/*H6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H6_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H7_Pool_L1/*H7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H7_Pool_L1/*H7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H7_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H8_Pool_L1/*H8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H8_Pool_L1/*H8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H8_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H9_Pool_L1/*H9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H9_Pool_L1/*H9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H9_Pool_L1/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*A9_Pool_L2/*A9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*A9_Pool_L2/*A9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_A9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_A9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B11_Pool_L2/*B11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B11_Pool_L2/*B11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_B11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_B11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B9_Pool_L2/*B9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*B9_Pool_L2/*B9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_B9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_B9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C1_Pool_L2/*C1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C1_Pool_L2/*C1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C1_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C10_Pool_L2/*C10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C10_Pool_L2/*C10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C11_Pool_L2/*C11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C11_Pool_L2/*C11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C12_Pool_L2/*C12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C12_Pool_L2/*C12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C12_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C2_Pool_L2/*C2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C2_Pool_L2/*C2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C2_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C3_Pool_L2/*C3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C3_Pool_L2/*C3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C3_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C4_Pool_L2/*C4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C4_Pool_L2/*C4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C4_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C5_Pool_L2/*C5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C5_Pool_L2/*C5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C5_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C6_Pool_L2/*C6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C6_Pool_L2/*C6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C6_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C7_Pool_L2/*C7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C7_Pool_L2/*C7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C7_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C8_Pool_L2/*C8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C8_Pool_L2/*C8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C8_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C9_Pool_L2/*C9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*C9_Pool_L2/*C9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_C9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_C9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D1_Pool_L2/*D1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D1_Pool_L2/*D1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D1_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D10_Pool_L2/*D10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D10_Pool_L2/*D10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D11_Pool_L2/*D11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D11_Pool_L2/*D11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D12_Pool_L2/*D12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D12_Pool_L2/*D12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D12_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D2_Pool_L2/*D2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D2_Pool_L2/*D2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D2_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D3_Pool_L2/*D3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D3_Pool_L2/*D3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D3_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D4_Pool_L2/*D4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D4_Pool_L2/*D4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D4_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D5_Pool_L2/*D5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D5_Pool_L2/*D5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D5_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D6_Pool_L2/*D6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D6_Pool_L2/*D6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D6_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D7_Pool_L2/*D7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D7_Pool_L2/*D7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D7_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D8_Pool_L2/*D8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D8_Pool_L2/*D8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D8_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D9_Pool_L2/*D9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*D9_Pool_L2/*D9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_D9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_D9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E10_Pool_L2/*E10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E10_Pool_L2/*E10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_E10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_E10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E12_Pool_L2/*E12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*E12_Pool_L2/*E12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_E12_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_E12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F10_Pool_L2/*F10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F10_Pool_L2/*F10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F11_Pool_L2/*F11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F11_Pool_L2/*F11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F12_Pool_L2/*F12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F12_Pool_L2/*F12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F12_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F9_Pool_L2/*F9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*F9_Pool_L2/*F9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_F9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_F9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G1_Pool_L2/*G1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G1_Pool_L2/*G1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G1_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G10_Pool_L2/*G10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G10_Pool_L2/*G10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G11_Pool_L2/*G11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G11_Pool_L2/*G11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G2_Pool_L2/*G2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G2_Pool_L2/*G2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G2_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G3_Pool_L2/*G3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G3_Pool_L2/*G3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G3_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G4_Pool_L2/*G4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G4_Pool_L2/*G4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G4_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G5_Pool_L2/*G5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G5_Pool_L2/*G5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G5_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G6_Pool_L2/*G6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G6_Pool_L2/*G6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G6_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G7_Pool_L2/*G7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G7_Pool_L2/*G7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G7_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G8_Pool_L2/*G8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G8_Pool_L2/*G8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G8_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G9_Pool_L2/*G9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*G9_Pool_L2/*G9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_G9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_G9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H1_Pool_L2/*H1_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H1_Pool_L2/*H1_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H1_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H1.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H10_Pool_L2/*H10_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H10_Pool_L2/*H10_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H10_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H10.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H11_Pool_L2/*H11_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H11_Pool_L2/*H11_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H11_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H11.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H12_Pool_L2/*H12_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H12_Pool_L2/*H12_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H12_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H12.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H2_Pool_L2/*H2_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H2_Pool_L2/*H2_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H2_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H2.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H3_Pool_L2/*H3_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H3_Pool_L2/*H3_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H3_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H3.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H4_Pool_L2/*H4_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H4_Pool_L2/*H4_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H4_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H4.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H5_Pool_L2/*H5_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H5_Pool_L2/*H5_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H5_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H5.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H6_Pool_L2/*H6_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H6_Pool_L2/*H6_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H6_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H6.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H7_Pool_L2/*H7_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H7_Pool_L2/*H7_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H7_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H7.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H8_Pool_L2/*H8_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H8_Pool_L2/*H8_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H8_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H8.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90
process_radtags -1 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H9_Pool_L2/*H9_*_R1.cut.fastq.gz -2 /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/1.FILTERED_READS_2025/*H9_Pool_L2/*H9_*_R2.cut.fastq.gz -D -o /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_H9_Pool_L2/ -b /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/ressource/barcode_SPX/barcode_SPX_H9.txt --renz_1 ecoRI --renz_2 mseI -crq -t 90

#Submitted batch job 50521763
#Slurm Job_id=50521763 Name=SPX_demultiplexing_2025 Began, Queued time 00:00:11
#Slurm Job_id=50521763 Name=SPX_demultiplexing_2025 Ended, Run time 07:57:20, COMPLETED, ExitCode 0



### Move rem files

mkdir rem_files

base_dir="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025"

for subdir in "$base_dir"/2025_SPX_*/; do
    for file in "${subdir}"*.rem.*.fq.gz; do
        if [[ -e "$file" ]]; then
            mv "$file" "$base_dir/rem_files/"
        fi
    done
done



#######
### Nb Retained Reads After Process_radtags ###
#######

#in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/

rm process_radtag.summary.txt
for FILE in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025/2025_SPX_*/process*.log;
do
    WELL=`echo $FILE | awk -F "_" '{print $3}'`
    LANE=`echo $FILE | awk -F "_" '{print $5}' | sed 's&/process&&g'`
    awk 'NF==6 {print}' $FILE | awk 'NR>3 {print}'> $$
    awk -v well="$WELL" -v lane="$LANE" '{print well, lane, $0}' $$ >> process_radtag.summary.txt
done

#Copy ./process_radtags.summary.txt to computer



#######
### Merge paired-end ###
#######

mkdir ./RESEQUENCING/3.MERGED_PE_2025/

###
### ./RESEQUENCING/3.MERGED_PE_2025/concatenate_fastq_spx2025.py
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
### ./3.MERGED_PE_2025/python_concatenate_spx2025_lane1.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name SPX_concatenate_fastq_Lane1
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load python/3.12.1

base_dir="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/2.DEMULTIPLEXED_2025"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/3.MERGED_PE_2025"

for sub_dir in ${base_dir}/2025_SPX_*_L1
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
          python concatenate_fastq_spx2025.py -F ${forward} -R ${reverse} -O ${OUTPATH}/${sample_id}.merged_L1.fq.gz
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

#Submitted batch job 50903625
#Slurm Job_id=50903625 Name=SPX_concatenate_fastq_Lane1 Began, Queued time 00:00:11
#Slurm Job_id=50903625 Name=SPX_concatenate_fastq_Lane1 Ended, Run time 06:03:06, COMPLETED, ExitCode 0


### Idem but with _L2

#Submitted batch job 50903633
#Slurm Job_id=50903633 Name=SPX_concatenate_fastq_Lane1 Began, Queued time 00:00:23
#Slurm Job_id=50903633 Name=SPX_concatenate_fastq_Lane1 Ended, Run time 05:41:52, COMPLETED, ExitCode 0




#######
### Merge sequences from the 2 lanes ###
#######

mkdir ./RESEQUENCING/4.MERGED_LANES


###
### ./4.MERGED_LANES/merging_lanes.sh
###


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name merging_lanes
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for seq1 in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/3.MERGED_PE_2025/*_L1.fq.gz
do
    outname=$(basename "${seq1}" _L1.fq.gz)
    
    seq2="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/3.MERGED_PE_2025/${outname}_L2.fq.gz"
    echo "processing " ${outname}    
    cat <(zcat "$seq1") <(zcat "$seq2") | gzip > "/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/4.MERGED_LANES/${outname}_merged2025.fq.gz" 
    echo ${outname} " merging finished"

done

#Submitted batch job 50932316
#Slurm Job_id=50932316 Name=merging_lanes Began, Queued time 00:00:06
#Slurm Job_id=50932316 Name=merging_lanes Ended, Run time 06:39:27, COMPLETED, ExitCode 0



#######
### Merge 2024 & 2025 ###
#######

mkdir 5.MERGED_ALL_SEQUENCINGS


###
### ./5.MERGED_ALL_SEQUENCINGS/merging_SPX_resequencing.sh
###


#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name merging_2024_2025
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for seq1 in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/4.MERGED_LANES/*.merged_merged2025.fq.gz
do
    outname=$(basename "${seq1}" .merged_merged2025.fq.gz)
    
    seq2="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/3.MERGED_PE/${outname}.merged.fq.gz"
    echo "processing " ${outname}    
    cat <(zcat "$seq1") <(zcat "$seq2") | gzip > "/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS/${outname}_mergedALL.fq.gz" 
    echo ${outname} " merging finished"

done

#Submitted batch job 50940940
#Slurm Job_id=50940940 Name=merging_2024_2025 Began, Queued time 00:00:04
#Slurm Job_id=50940940 Name=merging_2024_2025 Ended, Run time 11:05:33, COMPLETED, ExitCode 0


#######
### Nb of reads/samples ###
#######

# ./RESEQUENCING/nbreads_merged.sh

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu 
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name nbreads
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr 
#SBATCH --export NONE

for file_ALL in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS/*.fq.gz; do
    echo -n "$file_ALL: " >> /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS/nbreads_ALLSEQ.txt
    zcat "$file_ALL" | wc -l | awk '{print $1/4}' >> /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS/nbreads_ALLSEQ.txt
done

for file_2025 in /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/4.MERGED_LANES/*.fq.gz; do
    echo -n "$file_2025: " >> /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/4.MERGED_LANES/nbreads_2025lanes.txt
    zcat "$file_2025" | wc -l | awk '{print $1/4}' >> nbreads_2025lanes.txt
done

#Submitted batch job 50996205
#Slurm Job_id=50996205 Name=nbreads Began, Queued time 00:00:07
#Slurm Job_id=50996205 Name=nbreads Ended, Run time 00:57:02, COMPLETED, ExitCode 0

## Copy in computer

##Check nbreads/indiv... etc



#################################################################################################################

      #########################################
      ##### 4 - SNP Assembly For Analysis #####
      #########################################

#################################################################################################################


mkdir ./RESEQUENCING/6.RUN_STACKS_2025
# The analysis for the concatenated sequences of SPX will happen in this directory


#######
### USTACKS - stacks assembly ###
#######

## Here we kept parameters with best #SNPs (2024 analysis): -M 5 -m 3

###
### ./RESEQUENCING/6.RUN_STACKS_2025/ustacks_analysis_m3M5_2025.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --job-name SPX_ustacks_m3M5
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53
module load r-light/4.4.1

OUTPATH1="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/"
base_dir1="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS"

rm -f ${OUTPATH1}/Stat_ustacks_m3M5_2025.log
rm -f ${OUTPATH1}/ustacks_err_m3M5_2025
rm -f ${OUTPATH1}/cov_m3M5_2025
rm -f ${OUTPATH1}/nbstacks_m3M5_2025

echo Ind nb_unique_stacks nb_merged_stacks mean_cov > ${OUTPATH1}/Stat_ustacks_m3M5_2025.log

i="0"

for sample in ${base_dir1}/*_mergedALL.fq.gz; do
    echo `basename ${sample}`
    ustacks -p 8 -f ${sample} -i ${i} -o ${OUTPATH1} -t gzfastq -M 5 -m 3 2> ${OUTPATH1}/ustacks_err_m3M5_2025
    i=$[${i}+1]
    awk '/Assembled/' ${OUTPATH1}/ustacks_err_m3M5_2025 | awk '{print $2,$5}' | tr -d "\n" | awk '{print $1,$3}' > ${OUTPATH1}/nbstacks_m3M5_2025
    awk '/Final coverage/ {print $3}' ${OUTPATH1}/ustacks_err_m3M5_2025 | sed -e 's/mean=//g' | sed -e 's/;//g'  > ${OUTPATH1}/cov_m3M5_2025
    echo ${sample} `cat ${OUTPATH1}/nbstacks_m3M5_2025` `cat ${OUTPATH1}/cov_m3M5_2025` >> ${OUTPATH1}/Stat_ustacks_m3M5_2025.log
done


#Submitted batch job 51043150
#Slurm Job_id=51043150 Name=SPX_ustacks_m3M5 Began, Queued time 00:00:21
#Slurm Job_id=51043150 Name=SPX_ustacks_m3M5 Ended, Run time 1-02:53:56, COMPLETED, ExitCode 0ssh c

## Mean cov = 16.67
## Mean nb unique stacks = 232137
## Mean nb merged stacks = 164556


## Copy in computer
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/Stat_ustacks_m3M5_2025.log ./


# ## In R do :

# TABLE = read.csv("Stat_ustacks_SPX_CONCAT.csv", header = T, sep=";")
# str(TABLE)
# TABLE = TABLE[order(TABLE[,2]),]
# names=unlist(lapply(strsplit(as.character(TABLE[,1]), "./", fixed=TRUE), function(x) x[2]))
# pdf("Indiv_FragCover_SPX_CONCAT.pdf", width=(5+round(nrow(TABLE)/5)), height=10)
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

# mean(TABLE[,4])


### EXPLICATION DU GRAPHE

# Ustacks: créer dabord des stacks uniques (XXXAXXX/XXXAXXX/... et XXXTXXX/XXXTXXX) puis merge ces stacks selon -M 2

# Ici on veut 3 reads minimum pour créer un stacks (-m 3). La conséquence est que avec la couverture moyenne de 7, 
# il y a des locus qui ne peuvent pas créer de stacks uniques et donc ne pas les merger. A la fin on obtient
# des individus avec peu de stacks mais qui sont que uniques et donc très surement homozygotes et donc ca va créer
# un pb pour htz~cov.

# Pour savoir lequel changer entre -M et -m, il faut faire varier l'un en fixant l'autre.



#######
### CSTACKS - build catalog ###
#######

####
#### Use Checkpoints for > 3 days job
####

## Here : -n 6 (nb of differences within an individual + 1)

###
### ./RESEQUENCING/6.RUN_STACKS_2025/cstacks_analysis_n6.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --job-name SPX_cstacks_defaultCONCAT
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr

module load gcc/12.3.0
module load stacks/2.53

OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/"
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

cstacks -o ${OUTPATH} -p 4 -n 6 $(cat ${LIST_FILE})


###
### ./RESEQUENCING/6.RUN_STACKS_2025/cstacks_checkpoints_n6.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 12:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --job-name SPX_checkpoints_cstacks
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr

source /dcsrsoft/spack/external/ckptslurmjob/scripts/ckpt_methods.sh
setup_ckpt

launch_app ./"cstacks_analysis_n6.sh"

#Don't forget to make the job executable
chmod +x ./cstacks_analysis_n6.sh

## Do this BEFORE launching a job
export SBATCH_OPEN_MODE="append"
export SBATCH_SIGNAL=B:USR1@60
sbatch cstacks_checkpoints_n6.sh


#Submitted batch job 51208382
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:52:01
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:58:48
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 04:11:20
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:59:04
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:09:58
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:59:07
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:09:57
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:58:58
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:02:07
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:59:15
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:10:04
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:58:54
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:10:05
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:59:05
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:09:57
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:59:16
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:09:58
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Requeued, Run time 11:58:53
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Began, Queued time 00:10:02
#Slurm Job_id=51208382 Name=SPX_checkpoints_cstacks Ended, Run time 08:09:30, COMPLETED, ExitCode 0

Final catalog contains 1693395 loci.


#######
### SSTACKS - match catalog ### ## ON CONCAT
#######

###
### ./RESEQUENCING/6.RUN_STACKS_2025/sstacks_analysis.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 72:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --job-name SPX_sstacks_resequencing
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr

module load gcc/12.3.0
module load stacks/2.53

OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025"

for forward in ${OUTPATH}/*.snps.tsv.gz
do
    if [[ ! $forward =~ catalog ]]
    then
        fwd=${forward/.snps.tsv.gz/}
        sstacks -p 4 -c ${OUTPATH} -s ${fwd} -o ${OUTPATH}/
    fi
done


###
### ./RESEQUENCING/6.RUN_STACKS_2025/sstacks_resequencing_checkpoints.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 6:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --job-name SPX_sstacks_resequencing_checkpoints
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr

source /dcsrsoft/spack/external/ckptslurmjob/scripts/ckpt_methods.sh

CKPT_LOG=chkpt_logs/checkpoint-$REAL_JOB_ID.log
CKPT_DIR=chkpt_files/checkpoint-$REAL_JOB_ID-files

setup_ckpt

launch_app ./"sstacks_analysis.sh"


#Don't forget to make the job executable
chmod +x ./sstacks_analysis.sh

## Do this BEFORE launching a job
export SBATCH_OPEN_MODE="append"
export SBATCH_SIGNAL=B:USR1@60
sbatch sstacks_resequencing_checkpoints.sh


#TRY THE SCRIPT OF THE HELPDESK: sstacks_resequencing_checkpoints-fix.sh


## Do this BEFORE launching a job
export SBATCH_OPEN_MODE="append"
export SBATCH_SIGNAL=B:USR1@60
sbatch sstacks_resequencing_checkpoints-fix.sh

#Submitted batch job 51378207
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:00:01
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:03
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:20
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:21
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:19
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:07
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:23
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:18
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:02
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:00
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:14
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:09
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:07
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:04
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:05
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:03
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:18
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:20
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:00
#Slurm Job_id=51378207 Name=SPX_sstacks_resequencing_checkpoints Failed, Run time 06:00:22, TIMEOUT, ExitCode 0
#Slurm Job_id=51573189 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:00:03
#Slurm Job_id=51573189 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:10:01
#Slurm Job_id=51573189 Name=SPX_sstacks_resequencing_checkpoints Began, Queued time 00:09:58
#Slurm Job_id=51573189 Name=SPX_sstacks_resequencing_checkpoints Ended, Run time 04:40:27, COMPLETED, ExitCode 0



#######
### TSV2BAM - convert files ###
#######

mkdir ./CATALOG

###
### ./RESEQUENCING/6.RUN_STACKS_2025/CATALOG/tsv2bam_SPX_resequencing.sh
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
#SBATCH --job-name SPX_tsv2bam_m3M5
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/SPX_popmap_resequencing.txt"

tsv2bam -P ${INPATH} -M ${POP_INFOS} -t 16
mv "${INPATH}"/*.bam "${OUTPATH}"/

#Submitted batch job 51583182
#Slurm Job_id=51583182 Name=SPX_tsv2bam_m3M5 Began, Queued time 00:00:01
#Slurm Job_id=51583182 Name=SPX_tsv2bam_m3M5 Ended, Run time 00:16:11, COMPLETED, ExitCode 0



#######
### GSTACKS - genotyping ###
#######

###
### ./RESEQUENCING/6.RUN_STACKS_2025/CATALOG/gstacks_SPX_resequencing.sh
###

#!/bin/bash

#SBATCH --account pchriste_default
#SBATCH --mail-user clara.castex@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --job-name SPX_gstacks_resequencing
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53


OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/SPX_popmap_resequencing.txt"

gstacks -P ${OUTPATH} -M ${POP_INFOS} -t 4

#Submitted batch job 51586644
#Slurm Job_id=51586644 Name=SPX_gstacks_resequencing Began, Queued time 00:00:20
#Slurm Job_id=51586644 Name=SPX_gstacks_resequencing Ended, Run time 01:53:18, COMPLETED, ExitCode 0


## RESULTS:
Genotyped 1620767 loci:
  effective per-sample coverage: mean=14.5x, stdev=8.1x, min=4.0x, max=59.7x
  mean number of sites per locus: 180.3
  a consistent phasing was found for 9219590 of out 9945388 (92.7%) diploid loci needing phasing



#######
### POPULATIONS - SNPs calling ###
#######

###
### ./RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/populations_SPX_resequencing_p14.sh
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
#SBATCH --job-name SPX_populations14_reseq
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53


INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/SPX_popmap_resequencing.txt"

populations -P ${INPATH} --popmap ${POP_INFOS} -O ${OUTPATH} -p 14 -r 0.6 -R 0.6 -f p_value -t 8 --vcf --fstats --max-obs-het 0.5 --write-single-snp



#### Test with -p 14

#Submitted batch job 51597892
#Slurm Job_id=51597892 Name=SPX_populations14_reseq Began, Queued time 00:00:12
#Slurm Job_id=51597892 Name=SPX_populations14_reseq Ended, Run time 00:33:13, COMPLETED, ExitCode 0


# RESULTS:
Removed 1602876 loci that did not pass sample/population constraints from 1620767 loci.
Kept 17891 loci, composed of 3238251 sites; 41029 of those sites were filtered, 17887 variant sites remained.
Mean genotyped sites per locus: 181.00bp (stderr 0.03).


#### Test with -p 10
#### ./RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/populations_SPX_resequencing_p10.sh
#### #SBATCH --job-name SPX_populations10_reseq
#### OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/"


#Submitted batch job 51598150
#Slurm Job_id=51598150 Name=SPX_populations10_reseq Began, Queued time 00:00:03
#Slurm Job_id=51598150 Name=SPX_populations10_reseq Ended, Run time 00:46:41, COMPLETED, ExitCode 0


# RESULTS:
Removed 1578451 loci that did not pass sample/population constraints from 1620767 loci.
Kept 42316 loci, composed of 7659162 sites; 54589 of those sites were filtered, 42305 variant sites remained.
Mean genotyped sites per locus: 181.00bp (stderr 0.02).


## KEEP -p 10



#################################################################################################################

      ################################################
      ##### 5 - SNP Filtering For Analysis - p10 #####
      ################################################

#################################################################################################################


mkdir ./RESEQUENCING/7.SNPS_FILTERING_2025
# The analysis for the concatenated sequences of SPX will happen in this directory


#######
### MIN COVERAGE
#######

module load gcc/12.3.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out ${OUTPATH}/SPX_mincov5

# OUTPATH:
After filtering, kept 42305 out of a possible 42305 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025

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

Changing_allFIELDS_from_doubleDots_Geno SPX_mincov5.recode.vcf SPX_mincov5_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPX_mincov5_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out SPX_mincov5_clean_70miss

# OUTPUT:
After filtering, kept 14647 out of a possible 42305 Sites


#######
### SNPS SHARING - 80%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPX_mincov5_clean.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out SPX_mincov5_clean_80miss

# OUTPUT:
After filtering, kept 303 out of a possible 42305 Sites


#######
### SNPS SHARING - 90%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPX_mincov5_clean.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out SPX_mincov5_clean_90miss

# OUTPUT:
After filtering, kept 12 out of a possible 42305 Sites



###############
###############
#### ALL THE NEXT ANALYSES ARE DONE ON THE 70% SHARED SNPS 
###############
###############



#######
### QUALITY CHECK WITHOUT MAF OR MAC
#######

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/SPX_mincov5_clean_70miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/SPX_mincov5_clean_70miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/SPX_mincov5_clean_70miss* ./Spx/New_sequencing/4.SNP_Filtering

### Mean depth site = 21.61
### Mean depth individuals = 18.84
### Nb individuals >25% missing data = 169
### ID = Ind_missingdata_SPX_25miss.txt

### Nb individuals >30% missing data = 149
### ID = Ind_missingdata_SPX_30miss.txt

### Nb individuals >50% missing data = 91
### ID = Ind_missingdata_SPX_50miss.txt


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/5.MERGED_ALL_SEQUENCINGS/nbreads_ALLSEQ.txt ./Spx/New_sequencing/4.SNP_Filtering



##############################################################################################################

## Redo from populations on the individuals with <25% missing data

## Now we have 287 individuals

# New population map in : /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss/SPX_popmap25miss.txt



#######
### POPULATIONS - SNPs calling ###
#######

###
### ./RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss/populations_SPX2025_25miss.sh
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
#SBATCH --job-name SPX_2025_pop25miss
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr
#SBATCH --export NONE

module load gcc/12.3.0
module load stacks/2.53


INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss"
POP_INFOS="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss/SPX_popmap25miss.txt"

populations -P ${INPATH} --popmap ${POP_INFOS} -O ${OUTPATH} -p 10 -r 0.6 -R 0.6 -f p_value -t 16 --vcf --fstats --max-obs-het 0.5 --write-single-snp

#Submitted batch job 51630389
#Slurm Job_id=51630389 Name=SPX_2025_pop25miss Began, Queued time 00:05:05
#Slurm Job_id=51630389 Name=SPX_2025_pop25miss Ended, Run time 00:52:24, COMPLETED, ExitCode 0


Removed 1527926 loci that did not pass sample/population constraints from 1620767 loci.
Kept 92841 loci, composed of 16801451 sites; 37860 of those sites were filtered, 92807 variant sites remained.
Mean genotyped sites per locus: 180.97bp (stderr 0.01).




mkdir ./RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss
# The analysis for the concatenated sequences of SPX will happen in this directory


#######
### MIN COVERAGE - 25miss
#######

module load gcc/12.3.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 5 --recode --recode-INFO-all --out ${OUTPATH}/SPXmiss_mincov5

# OUTPATH:
After filtering, kept 92807 out of a possible 92807 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss

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

Changing_allFIELDS_from_doubleDots_Geno SPXmiss_mincov5.recode.vcf SPXmiss_mincov5_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_70miss

# OUTPUT:
After filtering, kept 76859 out of a possible 92807 Sites


#######
### SNPS SHARING - 80%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss

# OUTPUT:
After filtering, kept 56707 out of a possible 92807 Sites


#######
### SNPS SHARING - 90%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_90miss

# OUTPUT:
After filtering, kept 20528 out of a possible 92807 Sites



###############
###############
#### ALL THE NEXT ANALYSES ARE DONE ON THE 80% SHARED SNPS 
###############
###############



#######
### QUALITY CHECK WITHOUT MAF OR MAC
#######

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss* ./Spx/New_sequencing/4.SNP_Filtering/Ind25miss

# Mean depth site = 21.78X
# Mean deapth individuals = 21.25X

## Still 32 individuals having more than 25% missing data
## First I will do the filtering with these individuals and then without and see the impact.



#######
### MINIMUM COV 10 / MAXIMUM COV 40
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss.recode.vcf --min-meanDP 10 --max-meanDP 40 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_cov10_65

# OUTPUT :
After filtering, kept 56465 out of a possible 56707 Sites


#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_cov10_65.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_cov10_65_hweq

# OUTPUTS :
After filtering, kept 37276 out of a possible 56465 Sites


#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_cov10_65_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results


#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_cov10_65_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_cov10_65_hweq_mac

# OUTPUTS :
After filtering, kept 10479 out of a possible 37276 Sites


#######
### MAC 5 - HETEROZYGOSITY & DEPTH
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss_cov10_65_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss_cov10_65_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss_cov10_65_hweq_mac* ./Spx/New_sequencing/5.Analysis
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPX_tajima_results* ./Spx/New_sequencing/5.Analysis




##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

### Redo analysis with minDP 10 because of the correlation heterozygosity/coverage.
### I will remove SNPs with low coverage so the heterozygosity will change.



mkdir ./RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10
# The analysis for the concatenated sequences of SPX will happen in this directory


#######
### MIN COVERAGE - 25miss
#######

module load gcc/12.3.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 10 --recode --recode-INFO-all --out ${OUTPATH}/SPXmiss_mincov10

# OUTPATH:
After filtering, kept 92807 out of a possible 92807 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10

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

Changing_allFIELDS_from_doubleDots_Geno SPXmiss_mincov10.recode.vcf SPXmiss_mincov10_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss

# OUTPUT:
After filtering, kept 37199 out of a possible 92807 Sites


#######
### SNPS SHARING - 80%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_80miss

# OUTPUT:
After filtering, kept 9165 out of a possible 92807 Sites


#######
### SNPS SHARING - 90%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_90miss

# OUTPUT:
After filtering, kept 172 out of a possible 92807 Sites



###############
###############
#### ALL THE NEXT ANALYSES ARE DONE ON THE 70% SHARED SNPS 
###############
###############



#######
### QUALITY CHECK WITHOUT MAF OR MAC
#######

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss* ./Spx/New_sequencing/4.SNP_Filtering/Ind25miss/minDP10

# Mean depth site = 25.69X
# Mean deapth individuals = 24.39X

## Still 116 individuals having more than 25% missing data
## First I will do the filtering with these individuals and then without and see the impact.



#######
### MINIMUM COV 10 / MAXIMUM COV 65
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss.recode.vcf --min-meanDP 10 --max-meanDP 50 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov10_50

# OUTPUT :
After filtering, kept 37071 out of a possible 37199 Sites



#######
### MINIMUM COV 15 / MAXIMUM COV 65
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss.recode.vcf --min-meanDP 15 --max-meanDP 50 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov15_50

# OUTPUT :
After filtering, kept 36881 out of a possible 37199 Sites



#######
### MINIMUM COV 20 / MAXIMUM COV 65
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss.recode.vcf --min-meanDP 20 --max-meanDP 50 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov20_50

# OUTPUT :
After filtering, kept 13521 out of a possible 37199 Sites



#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov10_50.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov10_50_hweq

# OUTPUTS :
After filtering, kept 30608 out of a possible 37071 Sites



#######
### HWE - min-meanDP 15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov15_50.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov15_50_hweq

# OUTPUTS :
After filtering, kept 30442 out of a possible 36881 Sites



#######
### HWE - min-meanDP 20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov20_50.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov20_50_hweq

# OUTPUTS :
After filtering, kept 11549 out of a possible 13521 Sites



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov10_50_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES - min-meanDP 15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov15_50_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_min15



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES - min-meanDP 20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov20_50_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_min20



#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov10_50_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov10_50_hweq_mac

# OUTPUTS :
After filtering, kept 9913 out of a possible 30608 Sites



#######
### MAC - 5 - min-meanDP 15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov15_50_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov15_50_hweq_mac

# OUTPUTS :
After filtering, kept 9857 out of a possible 30442 Sites



#######
### MAC - 5 - min-meanDP 20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_cov20_50_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_cov20_50_hweq_mac

# OUTPUTS :
After filtering, kept 3608 out of a possible 11549 Sites


#######
### MAC 5 - HETEROZYGOSITY & DEPTH
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov10_50_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov10_50_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov10_50_hweq_mac* ./Spx/New_sequencing/5.Analysis/minDP10
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPX_tajima_results* ./Spx/New_sequencing/5.Analysis/minDP10


#######
### MAC 5 - HETEROZYGOSITY & DEPTH - min-meanDP 15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov15_50_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov15_50_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov15_50_hweq_mac* ./Spx/New_sequencing/5.Analysis/minDP10
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPX_tajima_results* ./Spx/New_sequencing/5.Analysis/minDP10


#######
### MAC 5 - HETEROZYGOSITY & DEPTH - min-meanDP 20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov20_50_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov20_50_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 


# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss_cov20_50_hweq_mac* ./Spx/New_sequencing/5.Analysis/minDP10
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPX_tajima_results* ./Spx/New_sequencing/5.Analysis/minDP10




##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

### Redo analysis with minDP 15 because of the correlation heterozygosity/coverage.
### I will remove SNPs with low coverage so the heterozygosity will change.



mkdir ./RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15
# The analysis for the concatenated sequences of SPX will happen in this directory


#######
### MIN COVERAGE - 25miss
#######

module load gcc/12.3.0
module load vcftools/0.1.16

INPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/6.RUN_STACKS_2025/CATALOG/POPULATIONS/p10/Ind25miss"
OUTPATH="/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15"

vcftools --vcf ${INPATH}/populations.snps.vcf --minDP 15 --recode --recode-INFO-all --out ${OUTPATH}/SPXmiss_mincov15

# OUTPATH:
After filtering, kept 92807 out of a possible 92807 Sites


### After this step I should use the script of Eleonore to remove all the DP files from hte missing genotypes. ### Otherwise when I use mean-minDP it stil consider those sites and create artifacts

#######
### Change the DP field to 0 after min/maxDP filter
#######

#### ELU script to remove all the DP files from hte missing genotypes
#In the frontend
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15

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

Changing_allFIELDS_from_doubleDots_Geno SPXmiss_mincov15.recode.vcf SPXmiss_mincov15_clean.recode.vcf


#######
### SNPS SHARING - 70%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_70miss

# OUTPUT:
After filtering, kept 1913 out of a possible 92807 Sites


#######
### SNPS SHARING - 60%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss

# OUTPUT:
After filtering, kept 15711 out of a possible 92807 Sites


#######
### SNPS SHARING - 50%
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_50miss

# OUTPUT:
After filtering, kept 41691 out of a possible 92807 Sites



###############
###############
#### ALL THE NEXT ANALYSES ARE DONE ON THE 60% SHARED SNPS 
###############
###############



#######
### QUALITY CHECK WITHOUT MAF OR MAC
#######

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16


VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss


vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 
vcftools --gzvcf $VCF --out $OUT --depth 
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 
vcftools --gzvcf $VCF --out $OUT --site-quality 
vcftools --gzvcf $VCF --out $OUT --missing-indv 
vcftools --gzvcf $VCF --out $OUT --missing-site 
vcftools --gzvcf $VCF --out $OUT --het

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss* ./Spx/New_sequencing/4.SNP_Filtering/Ind25miss/minDP15

# Mean depth site = 31.11X
# Mean deapth individuals = 29.08X

## Still 173 individuals having more than 25% missing data
## First I will do the filtering with these individuals and then without and see the impact.



#######
### MINIMUM COV 10 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss.recode.vcf --min-meanDP 10 --max-meanDP 60 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov10_60

# OUTPUT :
After filtering, kept 15638 out of a possible 15711 Sites



#######
### MINIMUM COV 25 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss.recode.vcf --min-meanDP 25 --max-meanDP 60 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov25_60

# OUTPUT :
After filtering, kept 828 out of a possible 15711 Sites



#######
### MINIMUM COV 15 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss.recode.vcf --min-meanDP 15 --max-meanDP 60 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov15_60

# OUTPUT :
After filtering, kept 15638 out of a possible 15711 Sites



#######
### MINIMUM COV 20 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss.recode.vcf --min-meanDP 20 --max-meanDP 60 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov20_60

# OUTPUT :
After filtering, kept 5342 out of a possible 15711 Sites


#######
### HWE - min-meanDP15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov15_60.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov15_60_hweq

# OUTPUTS :
After filtering, kept 14379 out of a possible 15638 Sites


#######
### HWE - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov20_60.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov20_60_hweq

# OUTPUTS :
After filtering, kept 4940 out of a possible 5342 Sites



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES - min-meanDP15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov15_60_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_meanDP15



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov20_60_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_meanDP20



#######
### MAC - 5 - min-meanDP15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov15_60_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov15_60_hweq_mac

# OUTPUTS :
After filtering, kept 4510 out of a possible 14379 Sites



#######
### MAC - 5 - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_cov20_60_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_cov20_60_hweq_mac

# OUTPUTS :
After filtering, kept 1404 out of a possible 4940 Sites



#######
### MAC 5 - HETEROZYGOSITY & DEPTH - min-meanDP15
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov15_60_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov15_60_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth 



#######
### MAC 5 - HETEROZYGOSITY & DEPTH - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov20_60_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov20_60_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth
vcftools --gzvcf $VCF --out $OUT --missing-indv 

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov15_60_hweq_mac* ./Spx/New_sequencing/5.Analysis/minDP15
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPX_tajima_results* ./Spx/New_sequencing/5.Analysis/minDP15
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss_cov20_60_hweq_mac* ./Spx/New_sequencing/5.Analysis/minDP15



#################################################################################################################

      #######################################################################
      ##### 6 - SNP Filtering For Analysis - p10 - without missing data #####
      #######################################################################

#################################################################################################################


mkdir ./RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data


#################################
### minDP 5 - max-missing 0.8 ###
#################################

mkdir ./minDP5


#######
### KEEP IND WITHOUT MISSING DATA
#######

#From R get the list of individuals to keep and import it in the cluster

module load gcc/12.3.0
module load vcftools/0.1.16

SAMPLES=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/No_miss25_minDP5_TOKEEP.txt
VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/SPXmiss_mincov5_clean_80miss.recode.vcf
OUTPUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata

vcftools --vcf $VCF --keep $SAMPLES --recode --recode-INFO-all --out $OUTPUT

## MEAN DEPTH
module load gcc/12.3.0
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata

vcftools --vcf $VCF --out $OUT --depth

#Mean depth = 22.35174

scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/



#######
### MINIMUM COV 10 / MAXIMUM COV 45
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_nomissdata.recode.vcf --min-meanDP 10 --max-meanDP 45 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45

# OUTPUT :
After filtering, kept 56507 out of a possible 56707 Sites


#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq

# OUTPUTS :
After filtering, kept 39368 out of a possible 56507 Sites



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_meanDP5_nomissdata



#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq_mac

# OUTPUTS :
After filtering, kept 11104 out of a possible 39368 Sites



#######
### HETEROZYGOSITY & DEPTH - min-meanDP10
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth
vcftools --gzvcf $VCF --out $OUT --missing-indv 

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPXmiss_mincov5_clean_80miss_nomissdata_cov10_45_hweq_mac* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP5/SPX_tajima_results_meanDP5_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss





##################################
### minDP 10 - max-missing 0.7 ###
##################################

mkdir ./minDP10


#######
### KEEP IND WITHOUT MISSING DATA
#######

#From R get the list of individuals to keep and import it in the cluster

module load gcc/12.3.0
module load vcftools/0.1.16

SAMPLES=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/No_miss25_minDP10_TOKEEP.txt
VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP10/SPXmiss_mincov10_clean_70miss.recode.vcf
OUTPUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata

vcftools --vcf $VCF --keep $SAMPLES --recode --recode-INFO-all --out $OUTPUT

## MEAN DEPTH
module load gcc/12.3.0
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata

vcftools --vcf $VCF --out $OUT --depth

#Mean depth = 28.88543

scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP10



#######
### MINIMUM COV 20 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_nomissdata.recode.vcf --min-meanDP 20 --max-meanDP 60 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60

# OUTPUT :
After filtering, kept 36658 out of a possible 37199 Sites


#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq

# OUTPUTS :
After filtering, kept 32811 out of a possible 36658 Sites



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_meanDP10_nomissdata



#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac

# OUTPUTS :
After filtering, kept 9525 out of a possible 32811 Sites



#######
### HETEROZYGOSITY & DEPTH - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth
vcftools --gzvcf $VCF --out $OUT --missing-indv 

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPXmiss_mincov10_clean_70miss_nomissdata_cov20_60_hweq_mac* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP10
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP10/SPX_tajima_results_meanDP10_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP10


##################################
### minDP 15 - max-missing 0.6 ###
##################################

mkdir ./minDP15


#######
### KEEP IND WITHOUT MISSING DATA
#######

#From R get the list of individuals to keep and import it in the cluster

module load gcc/12.3.0
module load vcftools/0.1.16

SAMPLES=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/No_miss25_minDP15_TOKEEP.txt
VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/minDP15/SPXmiss_mincov15_clean_60miss.recode.vcf
OUTPUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata

vcftools --vcf $VCF --keep $SAMPLES --recode --recode-INFO-all --out $OUTPUT

## MEAN DEPTH
module load gcc/12.3.0
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata

vcftools --vcf $VCF --out $OUT --depth

#Mean depth = 36.84822

scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP15



#######
### MINIMUM COV 20 / MAXIMUM COV 60
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_nomissdata.recode.vcf --min-meanDP 20 --max-meanDP 80 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80

# OUTPUT :
After filtering, kept 15623 out of a possible 15711 Sites


#######
### HWE
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80.recode.vcf --hwe 0.05 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq

# OUTPUTS :
After filtering, kept 14996 out of a possible 15623 Sites



#######
### TAJIMA - TEST FOR EXCESS OF RARE ALLELES
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq.recode.vcf --TajimaD 1000 --out SPX_tajima_results_meanDP15_nomissdata



#######
### MAC - 5
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15

module load gcc/12.3.0
module load vcftools/0.1.16

vcftools --vcf SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq.recode.vcf --mac 5 --recode --recode-INFO-all --out SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq_mac

# OUTPUTS :
After filtering, kept 3559 out of a possible 14996 Sites



#######
### HETEROZYGOSITY & DEPTH - min-meanDP20
#######
#In /work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15

module load gcc/12.3.0
module load samtools/1.19.2
module load vcftools/0.1.16

VCF=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq_mac.recode.vcf
OUT=/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq_mac

vcftools --vcf $VCF --out $OUT --het
vcftools --vcf $VCF --out $OUT --depth
vcftools --gzvcf $VCF --out $OUT --site-mean-depth
vcftools --gzvcf $VCF --out $OUT --missing-indv 

# Copy the file on your desktop
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPXmiss_mincov15_clean_60miss_nomissdata_cov20_80_hweq_mac* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP15
scp ccastex@curnagl.dcsr.unil.ch:/work/FAC/FBM/DEE/pchriste/default/ccastex/RADseq_24/SPX/RESEQUENCING/7.SNPS_FILTERING_2025/Ind25miss/Nomiss_data/minDP15/SPX_tajima_results_meanDP15_nomissdata* ./Spx/New_sequencing/4.SNPS_Filtering/Ind25miss/minDP15







