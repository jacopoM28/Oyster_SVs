#!/bin/bash -l
#SBATCH -J pbsv.SVs
#SBATCH -o pbsv.SVs.output
#SBATCH -e pbsv.SVs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 70:00:00
#SBATCH -A naiss2024-22-5      
#SBATCH -p core
#SBATCH -n 16

module load conda bioinfo-tools samtools/1.14
conda activate Python_env

GENOME=$1
READS=$2
SPECIE=$3
TANDEM_BEDFILE=$4	#Tandem repeats annotation in bed format. 

pbmm2 align --preset CCS -j 16 "$GENOME" "$READS" "$GENOME".bam --log-level DEBUG --log-level INFO --sort --median-filter --sample sample1 #PacBio reads, add --CCS for HiFi reads
#Calculate median coverage, just to know it
samtools index -@ 6 "$GENOME".bam
mosdepth -t 16 -n --fast-mode --by 500 "$GENOME" "$GENOME".bam
#Get summary statistics of mapping
samtools flagstat -@ 4 "$GENOME".bam > "$GENOME".flagstat
#Identify SVs and call them
pbsv discover --tandem-repeats "$TANDEM_BEDFILE" "$GENOME".bam "$SPECIE".svsig.gz 
pbsv call --ccs --log-level TRACE -j 16 "$GENOME" "$SPECIE".svsig.gz "$SPECIE".var.vcf #--ccs for HiFi reads
#Extract heterozigous INS and DEL (ALT homozigous INDELS are likely artifacts)
grep DEL "$SPECIE".var.vcf | grep PASS | grep -v '1/1' | awk '{print $1"\t"$2-1"\t"$2+length($4)-1"\t"$3}' | awk '{print $0"\t"$3-$2"\t1"}' > "$SPECIE".var.DEL.6field.bed
grep INS "$SPECIE".var.vcf | grep PASS | grep -v -P '1/1' | awk '{print $1"\t"$2-1"\t"$2+length($5)-1"\t"$3}' | awk '{print $0"\t"$3-$2"\t1"}' > "$SPECIE".var.INS.6field.bed
