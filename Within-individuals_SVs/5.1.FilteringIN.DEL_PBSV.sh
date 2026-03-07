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
sed -i 's/SAMPLE/SNIFFLES_MNP/'  "${INDEL3/.vcf/.filtered.vcf}"
sed -i 's/SAMPLE/SNIFFLES_NGLMR/'  "${INDEL4/.vcf/.filtered.vcf}"
sed -i 's/sample1/PBSV_MNP/'  "${INDEL1/.vcf/.filtered.vcf}"
sed -i 's/sample1/PBSV_NGLMR/'  "${INDEL2/.vcf/.filtered.vcf}"
##Add a unique identified in the ID filed of PBSV-generated VCF
sed -i 's/\tpbsv/\tMNP_pbsv/'  "${INDEL1/.vcf/.filtered.vcf}"
sed -i 's/\tpbsv/\tNGLMR_pbsv/' "${INDEL2/.vcf/.filtered.vcf}"

ls "${INDEL1/.vcf/.filtered.vcf}" > vcf_ToMerge.txt
ls "${INDEL2/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt
ls "${INDEL4/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt
ls "${INDEL3/.vcf/.filtered.vcf}" >> vcf_ToMerge.txt

#Merge VCF files keeping only variants called by 3 callers
SURVIVOR merge  vcf_ToMerge.txt 1000 3 1 0 1 0 "$SPECIE"_Merged.FINAL.vcf
#Filter VCF files keeping only variants longer then 50bp and supported by at least four reads
SURVIVOR filter "$SPECIE"_Merged.FINAL.vcf NA 50 -1 0 4 "$SPECIE"_Merged.FINAL.FILTERED.vcf

#Recover PBSV IDs from merged SVs set and create a new VCF files with PBSV variants supported by three VCF
python /crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Scripts/Extract_PBSV_IDs.py --vcf "$SPECIE"_Merged.FINAL.FILTERED.vcf --col1 10 --col2 11 > PBSV_Confirmed.IDs.txt
grep "#" "$INDEL1" > "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf
grep "MNP_pbsv"  PBSV_Confirmed.IDs.txt | grep -w -Ff - "${INDEL1/.vcf/.filtered.vcf}" >> "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf
grep "#" "$INDEL2" > "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf
grep "NGLMR_pbsv"  PBSV_Confirmed.IDs.txt | grep -w -Ff - "${INDEL2/.vcf/.filtered.vcf}" >> "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf
bgzip -c "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf > "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf.gz
bcftools index "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf.gz
bgzip -c "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf > "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf.gz
bcftools index "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf.gz
bcftools concat --allow-overlaps "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp2.vcf.gz "$SPECIE"_Merged.FINAL.FILTERED.PBSV.tmp1.vcf.gz > "$SPECIE"_Merged.FINAL.FILTERED.PBSV.vcf