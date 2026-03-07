#!/bin/bash -l
#SBATCH -J SNPs
#SBATCH -o SNPs.output
#SBATCH -e SNPs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
#SBATCH -A snic2022-5-613	
#SBATCH -p core
#SBATCH -n 1

VCF=$1
SPECIE=$1

module load bioinfo-tools plink2 vcftools/0.1.16 bcftools/1.17

#Gzip and index platypus vcf
bgzip -c "$VCF" > "$VCF".gz
tabix -p vcf "$VCF".gz
#Keep only biallelic SNPs
bcftools view -m2 -M2 -v snps -Oz -o "${VCF/.vcf.gz/.biallelic.vcf.gz}" "$VCF"
tabix -p vcf "${VCF/.vcf.gz/.biallelic.vcf.gz}"
#Keep only SNPs called on assembled chromosomes and remove sites with FILTER != PASS or .
bcftools view --apply-filters .,PASS --regions CM035811.1,CM035812.1,CM035813.1,CM035814.1,CM035815.1,CM035816.1,CM035817.1,CM035818.1,CM035819.1,CM035820.1 -Oz -o "${VCF/.vcf.gz/.biallelic.ChrPASSED.vcf.gz}" "${VCF/.vcf.gz/.biallelic.vcf.gz}"
tabix -p vcf "${VCF/.vcf.gz/.biallelic.ChrPASSED.vcf.gz}"
#Calculate numbero of heterozigous variants
plink2 --out "${VCF/.vcf.gz/.biallelic.ChrPASSED}" --gzvcf "${VCF/.vcf.gz/.biallelic.ChrPASSED.vcf.gz}" --sample-counts cols=homalt,het --allow-extra-chr
#Calculate MAF on variants genotyped in at least 30% of the samples
vcftools --max-missing 0.3 --gzvcf "${VCF/.vcf.gz/.biallelic.ChrPASSED.vcf.gz}" --freq2 --out "${VCF/.vcf.gz/.biallelic.ChrPASSED}" --max-alleles 2
#Final filtering of individual variants based on MAF (min MAF=0.05) and number of succesfully genotyped samples (at least 30%) for pop gen analyses
vcftools --gzvcf "${VCF/.vcf.gz/.biallelic.ChrPASSED.vcf.gz}" --maf 0.05 --max-missing 0.3 --recode --stdout > "${VCF/.vcf.gz/.biallelic.ChrPASSED_POPGen.vcf}"

