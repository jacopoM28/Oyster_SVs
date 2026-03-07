#!/bin/bash -l
#SBATCH -J SNIFFLEs
#SBATCH -o SNIFFLEs.output
#SBATCH -e SNIFFLEs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 100:00:00
#SBATCH -A naiss2024-22-5       
#SBATCH -p core
#SBATCH -n 16

#Protocol from https://www.nature.com/articles/s41477-019-0507-8

GENOME=$1
READS=$2
SPECIE=$3
TANDEM=$4

module load conda bioinfo-tools minimap2/2.24-r1122 samtools/1.14
conda activate ngmlr

#Align reads with minimap2
#PacBio Reads RSII
minimap2 -t 16 --cs --MD -ax map-pb "$GENOME" "$READS" | samtools view -Sb - > "$SPECIE".minimap.bam
#PacBio Reads HiFi 
minimap2 -t 16 --cs --MD -ax map-hifi "$GENOME" "$READS" | samtools view -Sb - > "$SPECIE".minimap.bam
samtools sort -@6 "$SPECIE".minimap.bam > "$SPECIE".minimap.sorted.bam
#rm "$SPECIE".minimap.bam
samtools index "$SPECIE".minimap.sorted.bam 
#Call SVs with long reads aligned with Minimap2
sniffles --tandem-repeats "$TANDEM" -i "$SPECIE".minimap.sorted.bam -v "$SPECIE".SNIFFLES.minimap.vcf -t 16
