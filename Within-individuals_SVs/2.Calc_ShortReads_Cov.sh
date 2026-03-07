#!/bin/bash -l
#SBATCH -J sr-cov
#SBATCH -o sr-cov.output
#SBATCH -e sr-cov.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH -A snic2022-22-1207       
#SBATCH -p core
#SBATCH -n 16

GENOME=$1
READS1=$2
READS2=$3
SPECIE=$4
DEL_BAM=$5
ADAPTERS=$6

module load conda bioinfo-tools picard/2.27.5 bcftools/1.14 BEDTools/2.29.2 BEDOPS/2.4.39 jellyfish/2.3.0 bbmap/38.08 samtools/1.14 bwa/0.7.17

#Prepare output dir
mkdir -p /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/Hemizigosity/"$SPECIE"/Sr_Cov
#Clean reads
#bbduk.sh minlen=35 in1="$READS1" in2="$READS2" out1="${READS1/.fastq/.trimmed.fastq}" out2="${READS2/.fastq/.trimmed.fastq}" ref="$ADAPTERS" ktrim=r k=25 mink=11 hdist=1 qtrim=r trimq=30 tpe tbo threads=16
#Map reads, format conversion and extraction of mapped reads
bwa index "$GENOME" "${GENOME/.fa/}"
bwa mem -t 16 "$GENOME" "$READS1" "$READS2" | samtools sort -@16 -o "$SPECIE"-sr.sorted.bam -
#"${READS1/.fastq/.trimmed.fastq}" "${READS2/.fastq/.trimmed.fastq}" | samtools sort -@16 -o "$SPECIE"-sr.sorted.bam -
samtools index "$SPECIE"-sr.sorted.bam
samtools view -@ 8 -F 4 -h "$SPECIE"-sr.sorted.bam > "$SPECIE"-sr.sorted.sam
reformat.sh in="$SPECIE"-sr.sorted.sam out="$SPECIE"-sr.mapped.fastq
#Coun kmers and produce histograms of mapped reads
jellyfish count -t 8 -C -m 24 -s 16G "$SPECIE"-sr.mapped.fastq -o "$SPECIE"-sr.mapped.jf
jellyfish histo -o "$SPECIE"-sr.mapped.histo "$SPECIE"-sr.mapped.jf
#Extract reads that match to deletions
samtools view -@ 8 -F 4 -h -b -L $DEL_BAM "$SPECIE"-sr.sorted.bam > "$SPECIE"-sr.DEL.sorted.bam
bedmap --echo --fraction-map 1 <(bam2bed <"$SPECIE"-sr.DEL.sorted.bam) "$SPECIE"-sr.sorted.bam > "$SPECIE"-sr.DEL.sorted.bed
cut -f1-6 "$SPECIE"-sr.DEL.sorted.bed > "$SPECIE"-sr.DEL.sorted.6field.bed
fastaFromBed -fi $GENOME -bed "$SPECIE"-sr.DEL.sorted.6field.bed -fo "$SPECIE"-sr.DEL.fa
#Coun kmers and produce histograms of reads mapped to deletions
jellyfish count -t 8 -C -m 24 -s 16G "$SPECIE"-sr.DEL.fa  -o "$SPECIE"-sr.DEL.jf
jellyfish histo -o "$SPECIE"-sr.DEL.histo "$SPECIE"-sr.DEL.jf
#BAM files of reads mapped to deletions and infer short reads coverage
cut -f4 "$SPECIE"-sr.DEL.sorted.bed | sort -u > "$SPECIE"-sr.DEL.sorted.names
samtools view "$SPECIE"-sr.sorted.bam | fgrep -w -f "$SPECIE"-sr.DEL.sorted.names > "$SPECIE"-sr.DEL.sorted.sam
samtools view -H "$SPECIE"-sr.sorted.bam > "$SPECIE"-sr.DEL.header.sam
cat "$SPECIE"-sr.DEL.sorted.sam >> "$SPECIE"-sr.DEL.header.sam
samtools view -@8 -S -b "$SPECIE"-sr.DEL.header.sam >"$SPECIE"-sr.DEL.header.bam
samtools index -@8 "$SPECIE"-sr.DEL.header.bam
conda activate Python_env
mosdepth -t 8 -m -b "$DEL_BAM" "$GENOME".DEL "$SPECIE"-sr.DEL.header.bam
#Whole genome and DEL specific reads count
mosdepth -t 8 -m -b 1000 "$GENOME"-sr "$SPECIE"-sr.sorted.bam
zcat "$GENOME"-sr.regions.bed.gz | awk 'BEGIN{OFS=FS="\t"}{$5=sprintf("%.0f",$4) }1' | cut -f5 | sort -n | uniq -c | sed -e 's/^ *//' -e 's/\ /\t/' | awk '{print $2"\t"$1}' | head -200 > "$SPECIE"-sr.coverage.txt
zcat "$GENOME".DEL.regions.bed.gz | awk 'BEGIN{OFS=FS="\t"}{$6=sprintf("%.0f",$5) }1' | cut -f6 | sort -n | uniq -c | sed -e 's/^ *//' -e 's/\ /\t/' | awk '{print $2"\t"$1}' | head -200 > "$SPECIE"-sr_DEL_coverage.txt
#Calculate Nucleotide level heterozigosity
varTMP=$( tail -1 "$GENOME"-sr.mosdepth.summary.txt | cut -f4) 
varCOV=$(echo "$varTMP*2" | bc)
bcftools mpileup --annotate INFO/AD,FORMAT/AD --threads 8 -Ou -f "$GENOME" "$SPECIE"-sr.sorted.bam | bcftools call --annotate GQ -Ou -mv | bcftools view --max-alleles 2 |  bcftools filter -s LowQual -e "%QUAL<20 || DP>$varCOV || MQBZ < -4 || DP<5" > "$SPECIE"-sr.SNP.vcf
