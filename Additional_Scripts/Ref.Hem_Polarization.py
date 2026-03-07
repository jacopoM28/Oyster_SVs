#!/usr/bin/env python
# coding: utf-8

# In[73]:


import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Merge SNIFFLES and SVJedy-graph genotypes of reference hemizygous variants  using two outgroup species and create a consensus set')
parser.add_argument('--vcf_svjedy_1', required=True, help='input SVJedy-graph vcf for sepcie 1')
parser.add_argument('--vcf_sniffles_1', required=True, help='input sniffles vcf for sepcie 1')
parser.add_argument('--vcf_svjedy_2', required=True, help='input SVJedy-graph vcf for sepcie 2')
parser.add_argument('--vcf_sniffles_2', required=True, help='input sniffles vcf for sepcie 2')
parser.add_argument('--specie_1', required=True, help='abbreviation for outgroup specie 1')
parser.add_argument('--specie_2', required=True, help='abbreviation for outgroup specie 2')
parser.add_argument('--ref_specie', required=True, help='abbreviation for reference specie')

args = parser.parse_args()
VCF_SVJEDY_1 = args.vcf_svjedy_1
VCF_SNIFFLES_1 = args.vcf_sniffles_1
OUT_SPECIE_1 = args.specie_1
VCF_SVJEDY_2 = args.vcf_svjedy_2
VCF_SNIFFLES_2 = args.vcf_sniffles_2
OUT_SPECIE_2 = args.specie_2
SPECIE = args.ref_specie

id_list = []
length_list = []
type_SVs_list = []
SVJedy_geno_1_list = []
Sniffles_geno_1_list = []
SVJedy_geno_2_list = []
Sniffles_geno_2_list = []


with open(VCF_SVJEDY_1) as f :
    for item in f:
        if item[0] != '#' :
            type_SV =  item.split("\t")[7].split(";")[1].split("=")[1]
            if type_SV == "INS" : 
                length = len(item.split("\t")[4])
            else :
                length = len(item.split("\t")[3])
            id_variants = item.split("\t")[2]
            SVJedy_geno_1 = item.split("\t")[9].split(":")[0]
            id_list.append(id_variants)
            SVJedy_geno_1_list.append(SVJedy_geno_1)
            length_list.append(length)
            type_SVs_list.append(type_SV)
            
with open(VCF_SNIFFLES_1) as f :
    for item in f:
        if item[0] != '#' :
            Sniffles_geno_1 = item.split("\t")[9].split(":")[0]
            Sniffles_geno_1_list.append(Sniffles_geno_1)
            
with open(VCF_SVJEDY_2) as f :
    for item in f:
        if item[0] != '#' :
            SVJedy_geno_2 = item.split("\t")[9].split(":")[0]
            SVJedy_geno_2_list.append(SVJedy_geno_2)
            
with open(VCF_SNIFFLES_2) as f :
    for item in f:
        if item[0] != '#' :
            Sniffles_geno_2 = item.split("\t")[9].split(":")[0]
            Sniffles_geno_2_list.append(Sniffles_geno_2)
            

dict_df = {'id': id_list, 'len' : length_list, 'Type' : type_SVs_list,
           '%s_SVJEDY' %OUT_SPECIE_1 : SVJedy_geno_1_list, 
           '%s_Sniffles' %OUT_SPECIE_1 : Sniffles_geno_1_list,
          '%s_SVJEDY' %OUT_SPECIE_2 : SVJedy_geno_2_list,
          '%s_Sniffles' %OUT_SPECIE_2 : Sniffles_geno_2_list}

df = pd.DataFrame(dict_df)

df_filtered = df[~df.eq("./.").any(1)]
for ind in df_filtered.index:
    if df['%s_SVJEDY' %OUT_SPECIE_1][ind] == df['%s_Sniffles' %OUT_SPECIE_1][ind] :
        if df['%s_SVJEDY' %OUT_SPECIE_2][ind] == df['%s_Sniffles' %OUT_SPECIE_2][ind] :
            if df['%s_SVJEDY' %OUT_SPECIE_1][ind] == df['%s_Sniffles' %OUT_SPECIE_2][ind] :
                Filter = "PASS"
                Genotype = df['%s_Sniffles' %OUT_SPECIE_2][ind]
            else :
                Filter = "DISCORDANT"
                Genotype = "NULL"
        else :
            Filter = "DISCORDANT"
            Genotype = "NULL"
    else :
        Filter = "DISCORDANT"
        Genotype = "NULL"
    df_filtered.loc[ind, 'FILTER'] = Filter
    df_filtered.loc[ind, 'Genotype'] = Genotype
    if df_filtered['Type'][ind] == "DEL" and df_filtered['Genotype'][ind] == "1/1" :
        Polarized = "INSERTION"
    elif df_filtered['Type'][ind] == "DEL" and df_filtered['Genotype'][ind] == "0/0" :
        Polarized = "DELETION"
    elif df_filtered['Type'][ind] == "INS" and df_filtered['Genotype'][ind] == "0/0" :
        Polarized = "INSERTION"
    elif df_filtered['Type'][ind] == "INS" and df_filtered['Genotype'][ind] == "1/1" :
        Polarized = "DELETION"
    elif df_filtered['Type'][ind] == "INS" and df_filtered['Genotype'][ind] == "0/1" :
        Polarized = "Shared"
    elif df_filtered['Type'][ind] == "DEL" and df_filtered['Genotype'][ind] == "0/1" :
        Polarized = "Shared"
    else :
        Polarized = "NULL"
    df_filtered.loc[ind, 'Polarized_SV'] = Polarized

df_filtered.to_csv('Genotypes.tsv', index=False,sep='\t', )
