#!/bin/bash -l
#SBATCH -J NGLMR
#SBATCH -o NGLMR.output
#SBATCH -e NGLMR.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 200:00:00
#SBATCH -A naiss2024-22-5        
#SBATCH -p core
#SBATCH -n 20

GENOME=$1
READS=$2
SPECIE=$3
TANDEM=$4

#Note: Ngmlr has problem with samtools version > 1.9
module load conda bioinfo-tools minimap2/2.16 samtools/1.9
conda activate ngmlr

#Map reads with nglmr and reformat output
ngmlr -t 20 -x pacbio -r "$GENOME" -q "$READS" --bam-fix | samtools view -bh > "$GENOME".ngmlr.bam #PacBio reads
samtools sort -@ 6 "$GENOME".ngmlr.bam > "$GENOME".ngmlr.sort.bam
rm "$GENOME".ngmlr.bam
samtools index "$GENOME".ngmlr.sort.bam
#Call SVs with sniffles
sniffles --tandem-repeats "$TANDEM" -i "$GENOME".ngmlr.sort.bam -v "$SPECIE".SNIFFLES.nglmr.vcf -t 16
#Call SVs with pbsv
conda deactivate
conda activate Python_env
pbsv discover --sample sample1 --tandem-repeats "$TANDEM" "$GENOME".ngmlr.sort.bam "$SPECIE".ngmlr.svsig.gz
#pbsv call --log-level TRACE -j 16 "$GENOME" "$SPECIE".ngmlr.svsig.gz "$SPECIE".ngmlr.var.vcf
pbsv call --log-level TRACE -j 16 "$GENOME" --ccs "$SPECIE".ngmlr.svsig.gz "$SPECIE".ngmlr.var.vcf #HiFi reads

