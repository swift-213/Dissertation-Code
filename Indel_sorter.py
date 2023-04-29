import argparse 
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", help="Add in the VCF file under this flag", required=True)
parser.add_argument("-o", "--output", action='store', help="Add in a text file to get the dictionary output and the positions of the indels", required=False)
parser.add_argument("-o2", "--output2", action='store', required=False)

args = parser.parse_args()


#Inputs 
indels_in_CDS_length=args.input

#Outputs 
output = args.output
output2 = args.output2


file1 = open(indels_in_CDS_length, 'r')
Lines = file1.readlines()
#print(Lines)
indel_length=[]
for line in Lines:
    indel_length.append("{}".format(line.strip()))

indel_length_num=[int(x) for x in indel_length]


Indel_length_dict = dict.fromkeys([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100], 0)
Indel_length_dict_longer_than_100 = dict.fromkeys([], 0)


for line in indel_length_num:
    if line in Indel_length_dict:
        Indel_length_dict[line] += 1
    if int(line) > 100:
        if int(line) not in Indel_length_dict_longer_than_100.keys():
            Indel_length_dict_longer_than_100[int(line)] = 1
        else:
            Indel_length_dict_longer_than_100[int(line)]+=1

#for indel length dict up tp 100
sorted_Indel_length_dict = sorted(Indel_length_dict.items(), key=lambda x:x[1])
converted_sort = dict(sorted_Indel_length_dict)
print(converted_sort)

#for indel length dict greater than 100
sorted_Indel_length_dict_longer_than_100 = sorted(Indel_length_dict_longer_than_100.items(), key=lambda x:x[1])
converted_sort2 = dict(sorted_Indel_length_dict_longer_than_100)
print(converted_sort2)


#Output files
f = open(output, "w")
f.write("{\n")
for k in converted_sort.keys():
    f.write("'{}':'{}'\n".format(k, converted_sort[k]))
f.write("}")
f.close()

f = open(output2, "w")
f.write("{\n")
for k in converted_sort2.keys():
    f.write("'{}':'{}'\n".format(k, converted_sort2[k]))
f.write("}")
f.close()

