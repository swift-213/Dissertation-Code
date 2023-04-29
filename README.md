# Dissertation

Scripts I have written in the process my dissertation on Indel Diveristy

## Cloning and Requirements

To download the repository + install reqirements use:

```
github clone https://github.com/swift-213/Dissertation.git

cd Dissertation

python3 -m pip install -r dissertation_requirements.txt

```

## Assembly_info.py

### Usage:
```
python3 Assembly_info.py -i vcf -i2 autosomes -o ${Assembly_info} -o2 indel_positions -o3 indel_binner_100 -o4 indel_binner_10000 -o5 every_indel_length
```
The only required input is the VCF file. If you wish to specify specific contigs or chromosomes to be run use -i2 and use a textfile with each contig name on a new row.

### Inputs

| Flag | Input | Required | 
|-|-|-|
|-i| VCF file| Yes |
|-i2|Textfile specifying which contigs to use| Yes |

### Outputs
| Flag | Output | Printed to terminal | Output type | Required | 
|-|-|-|-|-|
|-o | Number of lines, number of homozygous sites, number of SNPs, number of Indels| Yes | Text file | No |
|-o2 | Chromosome, start, end position of all Indels | No | Textfile | No |
|-o3 | Indels binned by length from 1-100 | Yes | Textfile | No |
|-o4 | Number of indels in ranges 1-100, 101-200 etc | Yes | Textfile | No |  
|-o5 | Every specific Indel length binned | No | Textfile | No |

If no output file is specified with a flag the outputs with Yes in printed to terminal will print outputs to terminal but no outputs will be saved
