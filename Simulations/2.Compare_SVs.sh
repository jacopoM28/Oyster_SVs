#!bin/bash

SPECIE=$1	#SPECIE ABBREVIATION
TRUE_SET=$2	#SIMULATED SVs
PBSV_MINIMAP=$3	#IDENTIFIED SVs WITH PBSV + MINIMAP
PBSV_NGLMR=$4	#IDENTIFIED SVs WITH PBSV + NGLMR
SNIFFLES_MINIMAP=$5 #IDENTIFIED SVs WITH SNIFFLES + MINIMAP
SNIFFLES_NGLMR=$6 ##IDENTIFIED SVs WITH SNIFFLES + NGLMR
MERGED_VCF=$7	#MERGED AND FILTERED VCF

module load bioinfo-tools BEDTools/2.29.2 SURVIVOR/1.0.3
eval "$(conda shell.bash hook)"

#PBSV + Minimap2(bpmm2)
##Reformat pbsv output in bed file
grep "DEL" "$PBSV_MINIMAP" | grep "0/1" | awk -v OFS="\t" -F"\t" '{print$1,$2,$8}' | sed 's/IMPRECISE;//' | cut -d";" -f1,3 | sed 's/SVTYPE=//' | sed 's/SVLEN=-//' | sed 's/;/\t/' | awk 'BEGIN{ FS=OFS="\t" } { print $0, ($2+$4) }' | awk -v OFS="\t" -F"\t" '{print$1,$2,$5,$4}' > "$SPECIE"_Sim.var.bed
##Intersection of called SVs with true set (required 80% overlap)
bedtools intersect -a "$SPECIE"_Sim.var.bed -b "$TRUE_SET" -wa -wb  -f 0.80 -r > "$SPECIE"_Sim.var_VS_True.set_80.bed
##Intersection of called SVs with true set (required 90% overlap)
bedtools intersect -a "$SPECIE"_Sim.var.bed -b "$TRUE_SET" -wa -wb  -f 0.90 -r > "$SPECIE"_Sim.var_VS_True.set_90.bed
##Intersection of called SVs with true set (required 100% overlap)
bedtools intersect -a "$SPECIE"_Sim.var.bed -b "$TRUE_SET" -wa -wb  -f 0.99 -r > "$SPECIE"_Sim.var_VS_True.set_99.bed
##Calculating summary statistics on 80% of overlap between called and true set (Remember that we simulated 1000 deletions)
TOT_SVs=$( wc -l "$SPECIE"_Sim.var.bed | cut -d" " -f1 )
TRUE_POS=$( wc -l "$SPECIE"_Sim.var_VS_True.set_80.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_80=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )	#PRECISION (80%)
PBSV_RE_80=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )	#RECALL (80%)
##Calculating summary statistics on 90% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.var_VS_True.set_90.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_90=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' ) #PRECISION (80%)
PBSV_RE_90=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' ) #RECALL (90%)
##Calculating summary statistics on 100% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.var_VS_True.set_99.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_99=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' ) #PRECISION (100%)
PBSV_RE_99=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )	#RECALL (100%)

echo -e "SPECIE\tOVERLAP\tALIGNER\tCALLER\tPRECISION\tRECALL" > "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t80\tMinimap\tPBSV\t$PBSV_PR_80\t$PBSV_RE_80" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t90\tMinimap\tPBSV\t$PBSV_PR_90\t$PBSV_RE_90" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t99\tMinimap\tPBSV\t$PBSV_PR_99\t$PBSV_RE_99" >> "$SPECIE"_SimSVs_Stats.tsv

##PBSV + NGLMR
##Reformat pbsv output in bed file
grep "DEL" "$PBSV_NGLMR" | grep "0/1" | awk -v OFS="\t" -F"\t" '{print$1,$2,$8}' | sed 's/IMPRECISE;//' | cut -d";" -f1,3 | sed 's/SVTYPE=//' | sed 's/SVLEN=-//' | sed 's/;/\t/' | awk 'BEGIN{ FS=OFS="\t" } { print $0, ($2+$4) }' | awk -v OFS="\t" -F"\t" '{print$1,$2,$5,$4}' > "$SPECIE"_Sim.ngmlr.var.bed
##Intersection of called SVs with true set (required 80% overlap)
bedtools intersect -a "$SPECIE"_Sim.ngmlr.var.bed -b "$TRUE_SET" -wa -wb  -f 0.80 -r > "$SPECIE"_Sim.ngmlr.var_VS_True.set_80.bed
##Intersection of called SVs with true set (required 90% overlap)
bedtools intersect -a "$SPECIE"_Sim.ngmlr.var.bed -b "$TRUE_SET" -wa -wb  -f 0.90 -r > "$SPECIE"_Sim.ngmlr.var_VS_True.set_90.bed
##Intersection of called SVs with true set (required 100% overlap)
bedtools intersect -a "$SPECIE"_Sim.ngmlr.var.bed -b "$TRUE_SET" -wa -wb  -f 0.99 -r > "$SPECIE"_Sim.ngmlr.var_VS_True.set_99.bed
##Calculating summary statistics on 80% of overlap between called and true set (Remember that we simulated 500 insetions and 500 deletions)
TOT_SVs=$( wc -l "$SPECIE"_Sim.ngmlr.var.bed | cut -d" " -f1 )
TRUE_POS=$( wc -l "$SPECIE"_Sim.ngmlr.var_VS_True.set_80.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_80=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )  #PRECISION (80%)
PBSV_RE_80=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )  #RECALL (80%)
##Calculating summary statistics on 90% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.ngmlr.var_VS_True.set_90.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_90=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' ) #PRECISION (80%)
PBSV_RE_90=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' ) #RECALL (90%)
##Calculating summary statistics on 100% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.ngmlr.var_VS_True.set_99.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
PBSV_PR_99=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' ) #PRECISION (100%)
PBSV_RE_99=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' ) #RECALL (100%)

echo -e "$SPECIE\t80\tNGLMR\tPBSV\t$PBSV_PR_80\t$PBSV_RE_80" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t90\tNGLMR\tPBSV\t$PBSV_PR_90\t$PBSV_RE_90" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t99\tNGLMR\tPBSV\t$PBSV_PR_99\t$PBSV_RE_99" >> "$SPECIE"_SimSVs_Stats.tsv


#SNIFFLES + Minimap2
##Reformat sniffles output in bed file
grep "DEL" "$SNIFFLES_MINIMAP" | grep "0/1" | awk -v OFS="\t" -F"\t" '{print$1,$2,$8}' | sed 's/PRECISE;//' | cut -d";" -f1,2 | sed 's/SVTYPE=//' | sed 's/SVLEN=-//' | sed 's/;/\t/' | awk 'BEGIN{ FS=OFS="\t" } { print $0, ($2+$4) }' | awk -v OFS="\t" -F"\t" '{print$1,$2,$5,$3}' > "$SPECIE"_Sim.SNIFFLES.minimap.bed
##Intersection of called SVs with true set (required 80% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.minimap.bed -b "$TRUE_SET" -wa -wb  -f 0.80 -r > "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_80.bed
##Intersection of called SVs with true set (required 90% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.minimap.bed -b "$TRUE_SET" -wa -wb  -f 0.90 -r > "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_90.bed
##Intersection of called SVs with true set (required 100% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.minimap.bed -b "$TRUE_SET" -wa -wb  -f 0.99 -r > "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_99.bed
##Calculating summary statistics on 80% of overlap between called and true set
TOT_SVs=$( wc -l "$SPECIE"_Sim.SNIFFLES.minimap.bed | cut -d" " -f1 )
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_80.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_80=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )	#PRECISION
SNIFFLES_RE_80=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )	#RECALL
##Calculating summary statistics on 90% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_90.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_90=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
SNIFFLES_RE_90=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )
##Calculating summary statistics on 100% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.minimap_VS_True.set_99.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_99=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
SNIFFLES_RE_99=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )

echo -e "$SPECIE\t80\tMinimap\tSNIFFLES\t$SNIFFLES_PR_80\t$SNIFFLES_RE_80" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t90\tMinimap\tSNIFFLES\t$SNIFFLES_PR_90\t$SNIFFLES_RE_90" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t99\tMinimap\tSNIFFLES\t$SNIFFLES_PR_99\t$SNIFFLES_RE_99" >> "$SPECIE"_SimSVs_Stats.tsv

#SNIFFLES + NGLMR
grep "DEL" "$SNIFFLES_NGLMR" | grep "0/1" | awk -v OFS="\t" -F"\t" '{print$1,$2,$8}' | sed 's/PRECISE;//' | cut -d";" -f1,2 | sed 's/SVTYPE=//' | sed 's/SVLEN=-//' | sed 's/;/\t/' | awk 'BEGIN{ FS=OFS="\t" } { print $0, ($2+$4) }' | awk -v OFS="\t" -F"\t" '{print$1,$2,$5,$3}' > "$SPECIE"_Sim.SNIFFLES.nglmr.bed
##Intersection of called SVs with true set (required 80% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.nglmr.bed -b "$TRUE_SET" -wa -wb  -f 0.80 -r > "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_80.bed
##Intersection of called SVs with true set (required 90% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.nglmr.bed -b "$TRUE_SET" -wa -wb  -f 0.90 -r > "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_90.bed
##Intersection of called SVs with true set (required 100% overlap)
bedtools intersect -a "$SPECIE"_Sim.SNIFFLES.nglmr.bed -b "$TRUE_SET" -wa -wb  -f 0.99 -r > "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_99.bed
##Calculating summary statistics on 80% of overlap between called and true set
TOT_SVs=$( wc -l "$SPECIE"_Sim.SNIFFLES.nglmr.bed | cut -d" " -f1 )
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_80.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_80=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )      #PRECISION
SNIFFLES_RE_80=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )      #RECALL
##Calculating summary statistics on 90% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_90.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_90=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
SNIFFLES_RE_90=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )
##Calculating summary statistics on 100% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.SNIFFLES.nglmr_VS_True.set_99.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
SNIFFLES_PR_99=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
SNIFFLES_RE_99=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )

echo -e "$SPECIE\t80\tNGLMR\tSNIFFLES\t$SNIFFLES_PR_80\t$SNIFFLES_RE_80" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t90\tNGLMR\tSNIFFLES\t$SNIFFLES_PR_90\t$SNIFFLES_RE_90" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t99\tNGLMR\tSNIFFLES\t$SNIFFLES_PR_99\t$SNIFFLES_RE_99" >> "$SPECIE"_SimSVs_Stats.tsv

#FINAL SET
##Reformat VCF file 
grep "DEL" "$MERGED_VCF" | awk -v OFS="\t" -F"\t" '{print$1,$2,$8}' | cut -d";" -f1,7 | sed 's/SUPP=.;//' | sed 's/END=//' | sed 's/$/\tMERGED/' > "${MERGED_VCF/.vcf/.bed}"
##Intersection of called SVs with true set (required 80% overlap)
bedtools intersect -a "${MERGED_VCF/.vcf/.bed}" -b "$TRUE_SET" -wa -wb  -f 0.80 -r > "$SPECIE"_Sim.Merged_VS_True.set_80.bed
##Intersection of called SVs with true set (required 90% overlap)
bedtools intersect -a "${MERGED_VCF/.vcf/.bed}" -b "$TRUE_SET" -wa -wb  -f 0.90 -r > "$SPECIE"_Sim.Merged_VS_True.set_90.bed
##Intersection of called SVs with true set (required 100% overlap)
bedtools intersect -a "${MERGED_VCF/.vcf/.bed}" -b "$TRUE_SET" -wa -wb  -f 0.99 -r > "$SPECIE"_Sim.Merged_VS_True.set_99.bed
##Calculating summary statistics on 80% of overlap between called and true set
TOT_SVs=$( wc -l "${MERGED_VCF/.vcf/.bed}" | cut -d" " -f1 )
TRUE_POS=$( wc -l "$SPECIE"_Sim.Merged_VS_True.set_80.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
MERGED_PR_80=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )      #PRECISION
MERGED_RE_80=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )      #RECALL
##Calculating summary statistics on 90% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.Merged_VS_True.set_90.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
MERGED_PR_90=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
MERGED_RE_90=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )
##Calculating summary statistics on 100% of overlap between called and true set
TRUE_POS=$( wc -l "$SPECIE"_Sim.Merged_VS_True.set_99.bed | cut -d" " -f1 )
FALSE_POS=$(( TOT_SVs - TRUE_POS ))
DEN=$(( TRUE_POS + FALSE_POS ))
MERGED_PR_99=$( echo "$TRUE_POS / $DEN" | bc -l | sed 's/^/0/' )
MERGED_RE_99=$( echo "$TRUE_POS / 1000" | bc -l | sed 's/^/0/' )

echo -e "$SPECIE\t80\tMERGED\tMERGED\t$MERGED_PR_80\t$MERGED_RE_80" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t90\tMERGED\tMERGED\t$MERGED_PR_90\t$MERGED_RE_90" >> "$SPECIE"_SimSVs_Stats.tsv
echo -e "$SPECIE\t99\tMERGED\tMERGED\t$MERGED_PR_99\t$MERGED_RE_99" >> "$SPECIE"_SimSVs_Stats.tsv
