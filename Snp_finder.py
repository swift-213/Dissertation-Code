import itertools
from pysam import VariantFile
import argparse 
from argparse import ArgumentParser
import numpy as np



parser = ArgumentParser()
parser.add_argument("-i", "--input", help="Add in the VCF file under this flag", required=True)
parser.add_argument("-i2", "--input2", help="Add in the VCF file under this flag", required=False)
parser.add_argument("-o", "--output", action='store', help="Assembly information, SNP, INDELs, matches, number of lines", required=False)

args = parser.parse_args()


#Inputs 
vcf = VariantFile(args.input)
chrom_file=args.input2

#Outputs 
output = args.output


#Reading in the chrom_file and turning it into a table 
file1 = open(chrom_file, 'r')
Lines = file1.readlines()
chrom=[]
for line in Lines:
    chrom.append("{}".format(line.strip()))

# fetch the VCF columns as variants
all_variants = vcf.fetch()
variants = itertools.islice(all_variants, None)

#Set up parameters and dictionary for loop
Assembly_info = dict.fromkeys(['number_of_lines', 'matches', 'mismatches'], 0)
previous = -1



for variant in variants:
    if variant.chrom in chrom:
        #Get the number of lines
        if variant.pos != previous:
            Assembly_info['number_of_lines'] += 1
        #Get the number of matches
        if variant.alts == None and variant.info['DP'] != 0:
            Assembly_info['matches'] += 1
        if variant.alts != None and variant.info ['DP'] != 0:
            Assembly_info['mismatches'] += 1
    previous = variant.pos


print(Assembly_info)


table=[]
# print each data item.
for key, value in Assembly_info.items():
    table.append(value)
print(table)

#Assembly info output
with open(output, 'w') as output_file:
    for line in str(table):
        output_file.write(line)
        #output_file.write("\n")
    output_file.close()

