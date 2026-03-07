#!bin/bash

DIR=$1	#Directory with all MANTAs derived SVs
GENOME=$2	#Absolute path to Reference genome
SPECIE=$3	#Species abbreviation

mkdir /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph
WORKDIR=$( pwd )

for i in $( ls "$DIR" ); do
	#Create sample IDs from "${SAMPLE}".MANTA directories
	varSAMPLE=$( echo "$i" | cut -d"." -f1)
	cp "$DIR"/"$i"/results/variants/diploidSV.vcf.gz /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph/"$varSAMPLE".diploidSV.vcf.gz
	cp "$DIR"/"$i"/results/variants/diploidSV.vcf.gz.tbi /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph/"$varSAMPLE".diploidSV.vcf.gz.tbi
	cd /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph
	#Reformat from MANTA to Paragraph Compatible (https://gist.github.com/KamilSJaron/cc8ef01947c03fd80abe6bac7451469f)
	python /crex/proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Scripts/convertManta2Paragraph_compatible_vcf.py --quality_filtering -r 1000 "$GENOME" "$varSAMPLE".diploidSV.vcf.gz > "$varSAMPLE".diploidSV.Paragraph.vcf
	#Create unique IDs for each variant based on sample ID
	sed -i 's/\tManta/\t'"$varSAMPLE"'./' "$varSAMPLE".diploidSV.Paragraph.vcf
	cd "$WORKDIR"
done

cd /proj/sllstore2017073/private/Jacopo/Ostreida_SVs/Analyses/PopGen_SVs/"$SPECIE"/Reseq_SVs.Paragraph
#Sort and index Paragraph-ready vcfs files
for i in *Paragraph.vcf; do
	var=$( echo "$i" | cut -d"." -f1);
	bcftools sort -O z --output-file "$var".diploidSV.Paragraph.sorted.gz "$i";
	tabix -p vcf "$var".diploidSV.Paragraph.sorted.gz;
done;
#Split individual vcf files by chromosomes, to improve merging speed I will separatly merge all chromosome specific vcf files in parallel
cat *Paragraph.vcf | grep -v "#" | cut -f1 | sort -u > chr.txt
while read C; do 
	mkdir -p Parallel_Merging/"$C"; 
	for i in *Paragraph.vcf; do 
		var=$( echo "$i" | cut -d"." -f1); 
		bcftools view -O z -o Parallel_Merging/"$C"/"$var".${C}.vcf.gz "$var".diploidSV.Paragraph.sorted.gz "${C}"; 
	done; 
done < chr.txt
