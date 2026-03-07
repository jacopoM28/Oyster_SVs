#!/bin/bash

module load bioinfo-tools BEDTools/2.29.2

RESULTS_DIR=/crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/Cari
mkdir "$RESULTS_DIR"/Polimorphic_TEs

for i in $(ls | grep "Batch"); do 
	cd "$i";
	for sample in $(ls | grep "P591"); do
		varCOV=$(tail -n1 ../../../../Analyses/PopGen_SVs/Cari/Reseq_Cov/"$sample".mosdepth.summary.txt | cut -f4)
		MinSupportReads=$( echo "$varCOV/10" | bc | sed 's/\..*$//')
		awk '$7 == "1p1"' "$sample"/"$sample".insertion.bed | awk '$5 >= 0.2' | awk -v env_var="$MinSupportReads" '$8 > env_var' > "$RESULTS_DIR"/Polimorphic_TEs/"$sample".TEMP2.filtered.bed
	done;
	cd ../
done
	
cd "$RESULTS_DIR"/Polimorphic_TEs

mkdir Stats
mkdir Filtered_Bed
mkdir Reformatted_Bed

for i in *bed; do 
	sample=$( echo "$i" | cut -d"." -f1);
	Nins=$(wc -l "$i" | awk '{print$1}');
	varCOV=$(tail -n1 ../Reseq_Cov/"$sample".mosdepth.summary.txt | cut -f4)
	echo -e "$sample\t$Nins\t$varCOV" >> Stats/Nins_Summary.txt;
	awk -v var="$sample" '{print$1"\t"$2"\t"$3"\t"var"."$4"\t"0"\t"$6}' "$i" | sed 's/:.*\t0/\t0/' > "${i/.bed/.reformatted.bed}"
	mv "$i" Filtered_Bed
	mv "${i/.bed/.reformatted.bed}" Reformatted_Bed	
done;

cat Reformatted_Bed/*reformatted.bed | bedtools sort -i - > ALL_TEins.bed
bedtools cluster -s -d 50 -i ALL_TEins.bed > ALL_TEins.clustered.bed

while read line; do
	fam=$( echo "$line" | cut -f4 | cut -d"." -f2);
	varClass=$( grep -w "$fam" ../../../../Data/TE_Libs/ALL_Libs.nr.fa | cut -d"#" -f2);
	echo "$varClass";
done< ALL_TEins.clustered.bed
