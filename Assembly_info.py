import itertools
from pysam import VariantFile
import argparse 
from argparse import ArgumentParser
import numpy as np


#Definition in order to get us the library names!
def get_text_range_in_hundreds(x):
    lower = int(np.floor(x/100))
    upper = int(np.ceil(x/100))
    return str(lower*100 + 1) + "-" + str(upper*100)

parser = ArgumentParser()
parser.add_argument("-i", "--input", help="Add in the VCF file under this flag", required=True)
parser.add_argument("-i2", "--input2", help="Add in the VCF file under this flag", required=False)
parser.add_argument("-o", "--output", action='store', help="Assembly information, SNP, INDELs, matches, number of lines", required=False)
parser.add_argument("-o2", "--output2", action='store', help="Indel positions", required=False)
parser.add_argument("-o3", "--output3", action='store', help="Indel 100 dictionary", required=False)
parser.add_argument("-o4", "--output4", action='store', help="Indel 10,000 dictionary", required=False)
parser.add_argument("-o5", "--output5", action='store', help="Indel 10,000 dictionary", required=False)

args = parser.parse_args()


#Inputs 
vcf = VariantFile(args.input)
chrom_file=args.input2

#Outputs 
output = args.output
output2 = args.output2
output3 = args.output3
output4 = args.output4
output5 = args.output5


#Reading in the chrom_file and turning it into a table 
file1 = open(chrom_file, 'r')
Lines = file1.readlines()
#print(Lines)
chrom=[]
for line in Lines:
    chrom.append("{}".format(line.strip()))

#print(chrom)



# fetch the VCF columns as variants
all_variants = vcf.fetch()
variants = itertools.islice(all_variants, None)

Assembly_info = dict.fromkeys(['number_of_lines', 'matches', 'mismatches', 'indels'], 0)
Indel_length_dict_100 = dict.fromkeys([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100], 0)
Indel_length_dict_10000 = dict.fromkeys([], 0)
Every_indel_length_dict = dict.fromkeys([], 0)


previous = -1
previous_dp = -1
Indel_positions = []
Inside_indel = False
truncated_indel = False

for variant in variants:
    if variant.chrom in chrom:
        if truncated_indel == True:
            if variant.info['DP'] == 0:
                previous = variant.pos
                previous_dp = int(variant.info['DP']) 
                continue 
            else:
                truncated_indel = False
        #Get the number of lines
        if variant.pos != previous:
            Assembly_info['number_of_lines'] += 1
        #Get the number of matches
        if variant.alts == None and variant.info['DP'] != 0:
            Assembly_info['matches'] += 1
        if variant.alts != None and variant.info ['DP'] != 0:
            Assembly_info['mismatches'] += 1
        if int(variant.info['DP']) == 0 and variant.pos == previous + 1:
            if Inside_indel == False and truncated_indel == False:
                Assembly_info['indels'] += 1
                indel_chrom_start_end=[variant.chrom]
                indel_chrom_start_end.append(variant.pos)
                Inside_indel = True
        if int(variant.info['DP']) == 0 and variant.pos != previous + 1:
            if Inside_indel == True:
                Assembly_info['indels'] -= 1
                Inside_indel = False
            truncated_indel = True
            previous = variant.pos
            previous_dp = int(variant.info['DP']) 
            continue 
        if int(variant.info['DP']) != 0 and Inside_indel == True and truncated_indel == False: 
            if variant.pos == previous + 1:
                indel_chrom_start_end.append(previous)
                Indel_positions.append(indel_chrom_start_end)
                Inside_indel = False 
            else:
                indel_chrom_start_end=[]
                Inside_indel = False 
                Assembly_info['indels'] -= 1
    previous = variant.pos
    previous_dp = int(variant.info['DP']) 


#For getting indel_lengths 
Indel_length=[]
for chrom, start, end in Indel_positions:
    indel_length = (int(end) - int(start)) + 1
    Indel_length.append(indel_length)

#Feed indel length into here to get 100 and 100000 bp histogram outputs 
for chrom, start, end in Indel_length:
    if start in Indel_length_dict_100:
        Indel_length_dict_100[start] += 1
    ranges = get_text_range_in_hundreds(start, end)
    if ranges not in Indel_length_dict_10000.keys():
        Indel_length_dict_10000[ranges] = 1
    else:
        Indel_length_dict_10000[ranges]+=1
    if start in Every_indel_length_dict:
        Every_indel_length_dict[start] += 1
    else:
        Every_indel_length_dict[start] = 1

print(Assembly_info)
print(Indel_length_dict_100)
print(Indel_length_dict_10000)

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

#Indel_positions
with open(output2, 'w') as f:
    for item in Indel_positions:
        f.write(str(item) + "\n")
    f.close()

#100 dictionary output
f = open(output3, "w")
f.write("{\n")
for k in Indel_length_dict_100.keys():
    f.write("'{}':'{}'\n".format(k, Indel_length_dict_100[k]))
f.write("}")
f.close()

#10,000 dictionary output
f = open(output4, "w")
f.write("{\n")
for k in Indel_length_dict_10000.keys():
    f.write("'{}':'{}'\n".format(k, Indel_length_dict_10000[k]))
f.write("}")
f.close()

#every indel pos dictionary output
f = open(output5, "w")
f.write("{\n")
for k in Every_indel_length_dict.keys():
    f.write("'{}':'{}'\n".format(k, Every_indel_length_dict[k]))
f.write("}")
f.close()
