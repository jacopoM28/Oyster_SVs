module load bioinfo-tools SURVIVOR/1.0.3 bcftools/1.17

INDEL1=$1	#PBSV + MINIMAP
INDEL2=$2	#PBSV + NGLMR
INDEL3=$3	#SNIFFLES + MINIMAP
INDEL4=$4	#SNIFFLES + NGLMR
SPECIE=$5

#Filter and reformat original VCF
##Remove variants called as homozigous for the alternate allele
grep -v -P "1/1" "$INDEL1" > "${INDEL1/.vcf/.filtered.vcf}"
grep -v -P "1/1" "$INDEL2" > "${INDEL2/.vcf/.filtered.vcf}"
grep -v -P "1/1" "$INDEL3" > "${INDEL3/.vcf/.filtered.vcf}"
grep -v -P "1/1" "$INDEL4" > "${INDEL4/.vcf/.filtered.vcf}"
##Change sample field in each VCF based on the SVs caller and aligner
ls "${INDEL1/.vcf/.filtered.vcf}" > vcf_ToMerge.txt
ls "${INDEL2/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt
ls "${INDEL4/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt
ls "${INDEL3/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt

#Merge VCF files keeping only variants called by 3 callers
SURVIVOR merge  vcf_ToMerge.txt 1000 3 1 0 1 0 "$SPECIE"_Merged.vcf
#Filter VCF files keeping only variants longer then 50bp and supported by at least for reads
SURVIVOR filter "$SPECIE"_Merged.vcf NA 50 -1 0 4 "$SPECIE"_Merged.FILTERED.vcf
