#!/bin/bash

module load bioinfo-tools plink2 vcftools/0.1.16 bcftools/1.17

VCF=$1
TANDEM_BED=$2
SPECIE=$3

#First I want to set genotype of non-PASSED variants to missing (.), secondo I want to remove allv ariants called within tandem repeats
vcftools --remove-filtered-geno-all --exclude-bed "$TANDEM_BED" --vcf "${VCF/.vcf/.filtered.vcf}" --stdout --recode > "$SPECIE"_Recoded.NoTR.vcf
#Keep only INDELS
grep "#" "$SPECIE"_Recoded.NoTR.vcf > "$SPECIE"_Recoded.NoTR.INDEL.vcf.tmp
grep "DEL" "$SPECIE"_Recoded.NoTR.vcf | grep -v "#" >> "$SPECIE"_Recoded.NoTR.INDEL.vcf.tmp
grep "INS" "$SPECIE"_Recoded.NoTR.vcf  | grep -v "#" >> "$SPECIE"_Recoded.NoTR.INDEL.vcf.tmp
bcftools sort -Ou "$SPECIE"_Recoded.NoTR.INDEL.vcf.tmp -o "$SPECIE"_Recoded.NoTR.INDEL.vcf
rm "$SPECIE"_Recoded.NoTR.INDEL.vcf.tmp
#Now I want to remove all variants called as homozigous for the reference allele in all samples (likely false variants)
vcftools --non-ref-ac-any 1 --vcf "$SPECIE"_Recoded.NoTR.INDEL.vcf --stdout --recode > "$SPECIE"_Recoded.NoTR.NoHom.vcf
#calculating fraction of variants genotyped as 0/0 in all samples, number of insertions and deletions
varFILT=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom.vcf | wc -l)
varALL=$( grep -v "#" "$SPECIE"_Recoded.NoTR.INDEL.vcf | wc -l)
varINS=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom.vcf | grep -c "INS")
varDEL=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom.vcf | grep -c "DEL")
#Calculate summary statistics (missing genotypes for sample, MAF, Heterozigosoty rates)
vcftools --gzvcf "$SPECIE"_Recoded.NoTR.NoHom.vcf --missing-indv --out "$SPECIE"_Recoded.NoTR.NoHom
vcftools --gzvcf "$SPECIE"_Recoded.NoTR.NoHom.vcf --missing-site --out "$SPECIE"_Recoded.NoTR.NoHom
#Calculate MAF on variants genotyped in at least 30% of the samples
vcftools --max-missing 0.3 --vcf "$SPECIE"_Recoded.NoTR.NoHom.vcf --freq2 --out "$SPECIE"_Recoded.NoTR.NoHom --max-alleles 2
#Final filtering of individual variants based on MAF (min MAF=0.05) and number of succesfully genotyped samples (at least 30%) for pop gen analyses
vcftools --vcf "$SPECIE"_Recoded.NoTR.NoHom.vcf --maf 0.05 --max-missing 0.3 --recode --stdout > "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf
grep "#" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf > "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet_DEL.vcf
grep "DEL" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf | grep -v "#" >> "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet_DEL.vcf
grep "#" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf > "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet_INS.vcf
grep "INS" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf | grep -v "#" >> "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet_INS.vcf
#Calculate final statistics and print out a summary table
varPOPGEN=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf | wc -l)
varPOPGEN_DEL=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf | grep -c "DEL")
varPOPGEN_INS=$( grep -v "#" "$SPECIE"_Recoded.NoTR.NoHom_PoPGenSet.vcf | grep -c "INS")
echo -e "ToT.Variants\tFiltered.Variants\tFiltered.INS\tFiltered.DEL\tPoPGen.Variants\tPoPGen.DEL\tPoPGen.INS" > "$SPECIE"_Stats.tsv
echo -e "$varALL\t$varFILT\t$varINS\t$varDEL\t$varPOPGEN\t$varPOPGEN_DEL\t$varPOPGEN_INS" >> "$SPECIE"_Stats.tsv

#Putting order into chaos
mkdir Stats
mv "$SPECIE"_Stats.tsv Stats
mv "$SPECIE"_Recoded.NoTR.NoHom.frq Stats
mv "$SPECIE"_Recoded.NoTR.NoHom.imiss Stats
mv "$SPECIE"_Recoded.NoTR.NoHom.lmiss Stats
mv "$SPECIE"_Recoded.NoTR.NoHom.vcf.scount Stats
