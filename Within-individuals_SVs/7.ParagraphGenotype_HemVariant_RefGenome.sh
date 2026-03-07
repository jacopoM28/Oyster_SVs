#!/bin/bash -l
#SBATCH -J Ref.Hem.SVs_Geno
#SBATCH -o Ref.Hem.SVs_Geno.output
#SBATCH -e Ref.Hem.SVs_Geno.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A naiss2024-22-5
#SBATCH -p core
#SBATCH -n 16

####Re-genotype hemyzygous variants in the reference genome using short reads data coming from the same sample

module load conda bioinfo-tools SURVIVOR

GENOME=$1 #Reference genome
BAM=$2	#sorted bam files with reads mapped to the reference genome
OUT=$3	#Output prefix (i.e. Sample name)
SPECIE=$4	#Specie abbreviation
VCF=$5 #VCF files with variants closest than read length from end of scaffold removed
READ_LENGTH=$6
MEAN_GENCOV=$7

#Create Manifest file for Paragraph
BAM_PATH=$(realpath "$BAM") 
echo -e "id\tpath\tdepth\tread length" > "$OUT"_Manifest
echo -e "$OUT\t$BAM_PATH\t$MEAN_GENCOV\t$READ_LENGTH" >> "$OUT"_Manifest
#Run Paragraph
conda activate Paragraph
MAX_COV=$( echo "$MEAN_GENCOV*20" | bc )
multigrmpy.py -i "$VCF" -M "$MAX_COV" -t 16 -m "$OUT"_Manifest -r "$GENOME" -o Paragraph_OUT/"$OUT"
