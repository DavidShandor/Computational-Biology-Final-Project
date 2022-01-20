import os
import re
from Bio import SeqIO

## read file and import details
gene_bank_file = 'BS168.gb'
 # making sure that the path is valid
assert(os.path.exists(gene_bank_file))

with open(gene_bank_file, "r") as input_handle:
    gen = SeqIO.parse(input_handle, "genbank")
    print(type(gen))
    record_gb = next(gen) # content of 1st record

print(type(record_gb))


# Are there any letters that are not nucleotides?
pattern = '[^ACGTacgt]' # pattern of non-nucleotide
test_string = str(record_gb.seq)
result = re.match(pattern, test_string)
print(result)


features = record_gb.features
print(len(features))

Genome_dict = {
    'CDS':0,
    'Genes':0,
    'Other':0
}

gene_dict = []
protein_dict = []

def asd(_feature):
    lenght = _feature.end.position - _feature.start.position
    return {
        'feature': _feature,
        'len': lenght
    }

for f in features:
    if f.type == 'gene':
        Genome_dict['Genes'] +=1
        gene_dict.append(asd(f))
    elif f.type =='CDS': 
        Genome_dict['CDS'] +=1
        protein_dict.append(asd(f))
    else:
        Genome_dict['Other'] += 1
        gene_dict.append(asd(f))

