#!/bin/bash -l
#SBATCH -J TEMP2
#SBATCH -o TEMP2.output
#SBATCH -e TEMP2.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH -A snic2022-22-1207     
#SBATCH -p core
#SBATCH -n 18

READS1=$1
READS2=$2
GENOME=$3
TE_LIB=$4
TE_BED=$5
BAM=$6

module load python/3.11.4 bioinfo-tools BEDTools/2.29.2 samtools/1.19 bwa/0.7.17

/crex/proj/sllstore2017073/private/Jacopo/Software/TEMP2/TEMP2 insertion2 -d -l "$READS1" -r "$READS2" -i "$BAM" -g "$GENOME" -R "$TE_LIB" -t "$TE_BED" -c 18
