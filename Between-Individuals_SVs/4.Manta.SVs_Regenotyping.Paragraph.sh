#!/bin/bash -l
#SBATCH -J Manta.Paragraph.SVs
#SBATCH -o Manta.Paragraph.SVs.output
#SBATCH -e Manta.Paragraph.SVs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 240:00:00
#SBATCH -A snic2022-22-1207	
#SBATCH -p core
#SBATCH -n 16

module load conda

DIR=$1	#Directory with bam files
SPECIE=$2 #Specie ID
VCF=$3
CHR=$4
READ_LENGTH=$5
GENOME=$6

conda activate Paragraph

cd "$DIR"
mkdir "$CHR".MANTA.SVs_Paragraph.OUT
mkdir "$CHR".MANTA.SVs_Manifest

#Genotype the merged set of MANTA derived SVs separately for each sample
for i in *.bam; do
	varSAMPLE=$( echo "$i" | cut -d"." -f1)
	echo -e "id\tpath\tdepth\tread length" > "$varSAMPLE"_Manifest
	varCOV=$( tail -n 1 /crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_Cov/"$varSAMPLE".mosdepth.summary.txt | awk -F"\t" '{print$4}')
	varPATH=$( realpath "$i" )
echo -e "$varSAMPLE\t$varPATH\t$varCOV\t$READ_LENGTH" >> "$varSAMPLE"_Manifest
	MAX_COV=$( echo "$varCOV*20" | bc | sed 's/\..*$//')
	multigrmpy.py -i "$VCF" -m "$varSAMPLE"_Manifest -M "$MAX_COV" -r "$GENOME" -o "$CHR"."$varSAMPLE".MANTA.SVs_Paragraph.OUT  -t 20
	mv "$CHR"."$varSAMPLE".MANTA.SVs_Paragraph.OUT "$CHR".MANTA.SVs_Paragraph.OUT
	mv "$CHR"."$varSAMPLE"_Manifest "$CHR".MANTA.SVs_Manifest
done;

