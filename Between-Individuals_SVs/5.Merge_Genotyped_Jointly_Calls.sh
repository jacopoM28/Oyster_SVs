#!bin/bash

DIR=$1 #Directory in which to search for single-chromosomes and single individual VCFs. In my case I have four directories which results will be finally merged all together to create the final VCF file
SPECIE=$2

mkdir /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph/Merge_Genotyped_JointlyCalls

cd "$DIR"
#sort and index individual VCF files
for i in $(ls | grep "Paragraph.OUT"); do 
	chr=$( echo "$i" | cut -d"." -f1,2); 
	for j in $(ls "$i"); do 
		bcftools sort -o "$i"/"$j"/genotypes.sorted.vcf.gz -Oz "$i"/"$j"/genotypes.vcf.gz; 
		tabix -p vcf "$i"/"$j"/genotypes.sorted.vcf.gz; 
	done;
done

#Merge all single-individual VCF files in one multi-sample chromosome-specific VCF file
for i in $(ls | grep "Paragraph.OUT"); do 
	chr=$( echo "$i" | cut -d"." -f1,2); 
	for j in $(ls "$i"); do 
		realpath "$i"/"$j"/genotypes.sorted.vcf.gz >> /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph/Merge_Genotyped_JointlyCalls/"$DIR"_"$chr".vcf.list;
	done; 
done;
cd /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph/Merge_Genotyped_JointlyCalls
for i in "$DIR"*list; do
	bcftools merge -l "$i" -Oz -o "${i/.vcf.list/.merged.vcf.gz}";
done;

#Concat all chromosome-specific VCF files
bcftools concat -Oz -o "$DIR".merged.vcf.gz "$DIR"*.merged.vcf.gz

####NB: Now it is only necessary to merge all batch-specific VCF files to create the final, merged, multi-sample VCF
