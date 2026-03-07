#!/usr/bin/env python
###Extract PBSV IDs from VCF file created with SURVIVOR merge because of better accuracy of PBSV breakpoints estimation

import argparse

parser = argparse.ArgumentParser(description='Extract IDs of PBSV from VCF file after intra-sample merging.')
parser.add_argument('--vcf', required=True, help='input vcf')
parser.add_argument('--col1', required=True, help='PBSV + MINIMAP')
parser.add_argument('--col2', required=True, help='PBSV + NGLMR')

args = parser.parse_args()
VCF = args.vcf
COL1 = int(args.col1) - 1
COL2 = int(args.col2) - 1

with open(VCF) as f:
    for item in f:
        if item[0] != '#':
            var_1 = item.split("\t")[COL1].split(":")[7]
            if var_1 == "NaN" :
                var_2 = item.split("\t")[COL2].split(":")[7]
                print(var_2)
            else :
                print(var_1)
