import os
import re
from Bio import SeqIO

## read file and import details
january = 'CoronaJanuary2022.gb'
july = 'CoronaJuly2020.gb'
# making sure that the path is valid
assert (os.path.exists(january))
assert (os.path.exists(july))


def parser(files):
    with open(files, "r") as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        print(type(gen))
        record_gb = next(gen)  # content of 1st record

    print(type(record_gb))

    # Are there any letters that are not nucleotides?
    pattern = '[^ACGTacgt]'  # pattern of non-nucleotide
    test_string = str(record_gb.seq)
    result = re.match(pattern, test_string)
    print(result)
    features = record_gb.features
    print(len(features))




def main():
    parser(january)
    parser(july)





if __name__ == '__main__':
    main()
