import argparse 
from argparse import ArgumentParser
import pandas as pd
import numpy as np
parser=ArgumentParser()
parser.add_argument("-i", "--input", help="Add in the VCF file under this flag", required=True)
parser.add_argument("-o", "--output", action='store', help="Add in a text file to get the dictionary output and the positions of the indels", required=False)
parser.add_argument("-i2", "--input2", help="Add in the VCF file under this flag", required=False)
args=parser.parse_args()

fasta_file=args.input
chrom_file=args.input2
output=args.output

file1=open(chrom_file, 'r')

Lines=file1.readlines()

chrom=[]
for line in Lines:
    chrom.append("{}".format(line.strip()))



def parseFasta(string, makeUppercase=False):
    splitString = string.split(">")[1:]
    names = [s.split()[0] for s in splitString]
    seqs = [s[s.index("\n"):].replace("\n","").replace(" ","") for s in splitString]
    if makeUppercase: seqs = [s.upper() for s in seqs]
    return (names,seqs)


with open(fasta_file) as f:
    allText = f.read()
 

names, seqs = parseFasta(allText)


names_starts_and_ends = []


for i in range(len(names)):
    if names[i] in chrom:
        print(names[i])
        inside = False
        for j in range(len(seqs[i])):
            if inside == False:
                if seqs[i][j].islower():
                    inside= True
                    current = [names[i], j]
            else:
                if not seqs[i][j].islower():
                    current.append(j)
                    names_starts_and_ends.append(current)
                    inside = False
         



with open(output, 'w') as f:
    for item in names_starts_and_ends:
        f.write(str(item) + "\n")
    f.close()
   