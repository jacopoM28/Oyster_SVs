#!bin/bash
#Apply further filtering of hemizigous SVs and preper them for sort reads re-genotyping using Paragraph

module load bioinfo-tools SURVIVOR/1.0.3 bcftools/1.17

SPECIE=$1
VCF=$2
GENOME=$3
SAMPLE=$4	#Sample name to keep (i.e. I used always the first one)

eval "$(conda shell.bash hook)"

#Remove SVs closest that 1Kb from end of scaffolds
conda activate Python_env
/crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Scripts/Ending_Coords.py --genome "$GENOME" --out "$SPECIE"_ToExclude.bed --dist 1000
SURVIVOR filter "$VCF" "$SPECIE"_ToExclude.bed -1 -1 0 -1 "$SPECIE"_Min1000Ends.vcf
conda deactivate
#Convert SNIFFLES allele rapresentation from simbolic to sequence-specific ALT alleles
conda activate Paragraph
bayesTyperTools convertAllele -v "$SPECIE"_Min1000Ends.vcf -g "$GENOME" -o Converted_"$SPECIE"_Min1000Ends
#Remove unecessary fields that can crerate bugs with paragraph and keep only one "sample".
bcftools annotate -x FORMAT,^INFO/END,INFO/SEQ,INFO/SVTYPE Converted_"$SPECIE"_Min1000Ends.vcf > Converted_LessFIELDS_"$SPECIE"_Min1000Ends.vcf
#Remove BND,INV and DUP SVs. Keeping only INS and DEL
bcftools view -s "$SAMPLE" Converted_LessFIELDS_"$SPECIE"_Min1000Ends.vcf | grep -v "TYPE=DUP" | grep -v "TYPE=BND" | grep -v "TYPE=INV" > Converted_LessFIELDS_OneSample_"$SPECIE"_Min1000Ends.vcf
conda deactivate
#Fix padded bases of insertions called by SNIFFLES using custom python script or Paragraph won't work
conda activate Python_env
python /crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Scripts/Add_Padded_Bases_vcf.py --vcf Converted_LessFIELDS_OneSample_"$SPECIE"_Min1000Ends.vcf  --genome "$GENOME" > Converted_LessFIELDS_OneSample_"$SPECIE"_Min1000Ends_FINAL.vcf
conda deactivate
#Normalize the variants
bcftools norm -c s --rm-dup exact -f "$GENOME" Converted_LessFIELDS_OneSample_"$SPECIE"_Min1000Ends_FINAL.vcf > Converted_LessFIELDS_OneSample_"$SPECIE"_Min1000Ends_TRUE.vcf
