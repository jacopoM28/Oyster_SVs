#!/bin/bash -l
#SBATCH -J SNV_Calling
#SBATCH -o SNV_Calling.output
#SBATCH -e SNV_Calling.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 240:00:00
#SBATCH -A snic2022-22-1207	
#SBATCH -p core
#SBATCH -n 20

module load conda bioinfo-tools samtools/1.17

conda activate Platypus

BAM_LIST=$1
OUT=$2
GENOME=$3

platypus callVariants --bamFiles="$1" --nCPU=11 --refFile="$GENOME" \
	--logFileName=platypus102.log --minMapQual=20 \
       	--minBaseQual=20 --maxVariants=10 --filterReadsWithUnmappedMates=0 --filterReadsWithDistantMates=0 \
       	--filterReadPairsWithSmallInserts=0 --output="$OUT" \
	--maxReads 50000000 --maxSize 50
