#!/bin/bash -l
#SBATCH -J Simulate.SVs
#SBATCH -o Simulate.SVs.output
#SBATCH -e Simulate.SVs.error
#SBATCH --mail-user jacopo.martelossi2@unibo.it
#SBATCH --mail-type=ALL
#SBATCH -t 70:00:00
#SBATCH -A snic2022-22-1207      
#SBATCH -p core
#SBATCH -n 12

CONFIG=$1	#Config file
OUT=$2	#Output prefix

module load bioinfo-tools BioPerl/1.6.924_Perl5.18.4

perl /crex/proj/sllstore2017073/private/Jacopo/Software/Sim-it/Sim-it1.3.4.pl -c "$CONFIG" -o "$OUT"
