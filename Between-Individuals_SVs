#!/bin/bash -l
#SBATCH -J Reseq.SVs
#SBATCH -o Reseq.SVs.output
#SBATCH -e Reseq.SVs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 150:00:00
#SBATCH -A snic2022-22-1207
#SBATCH -p core
#SBATCH -n 20

GENOME=$1
SPECIE=$2

module load bioinfo-tools manta/1.6.0

 
mkdir MANTA

for i in *bam; do
	configManta.py --bam "$i" --referenceFasta "$GENOME" --runDir MANTA/"${i/.sorted.bam/.MANTA}"
	MANTA/"${i/.sorted.bam/.MANTA}"/runWorkflow.py -j 20
done
