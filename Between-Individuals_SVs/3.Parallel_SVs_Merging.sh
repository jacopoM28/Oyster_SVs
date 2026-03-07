#!/bin/bash -l
#SBATCH -J SVs.Merging
#SBATCH -o SVs.Merging.output
#SBATCH -e SVs.Merging.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH -A snic2022-22-1207
#SBATCH -p core
#SBATCH -n 10

module load conda bioinfo-tools bcftools/1.17

input_dir=$1

conda activate Jasmine

gunzip "$input_dir"/*
realpath "$input_dir"/* > "$input_dir"/vcf_to_merge.tsv
cd "$input_dir"
jasmine threads=10 max_dist=20 --pre_normalize --keep_var_ids --clique_merging --normalize_type file_list=vcf_to_merge.tsv  out_file="$input_dir".Merged.SVs.vcf.tmp
bcftools annotate -x FORMAT,^INFO/END,INFO/SEQ,INFO/SVTYPE "$input_dir".Merged.SVs.vcf.tmp | bcftools view -s ^SAMPLE1 > "$input_dir".Merged.SVs.vcf
