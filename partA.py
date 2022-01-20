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

CDS = genes = regulators = 0

for f in features:
    if f.type == 'gene': genes +=1 
    elif f.type =='CDS': CDS +=1 
    else: regulators += 1