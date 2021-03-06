import os
from Bio import SeqIO
import pandas as pd
from genetic_data_file import table, Hidrophobic, Hidrophilic


class GeneticDataGenerator:
    """Contains all the genetic data form Genbank files and Uniprot """

    def read_uni_dataframe(self):
        if self.unifile is not None:
            assert (os.path.exists(self.unifile)), 'UniProtKB file does not exist  # sanity check'
            self.uni_df = pd.read_excel(self.unifile)

    def parse_genebank_file(self):
        """Parse the GeneBank file and return its features as DataFrame"""
        if self.genebank_file is not None:
            assert (os.path.exists(self.genebank_file))

            with open(self.genebank_file, "r") as input_handle:
                gen = SeqIO.parse(input_handle, "genbank")
                record_gb = next(gen)

            self.gb_features = record_gb.features
            self.sequence = str(record_gb.seq)

    def create_gb_dataframe(self):

        if self.gb_features is not None:

            gene_type = []
            location_start = []
            location_end = []
            location_strand = []
            qualifiers = {}

            for i, feature in enumerate(self.gb_features):
                gene_type.append(feature.type)
                location_start.append(feature.location.start.position)
                location_end.append(feature.location.end.position)
                location_strand.append(feature.strand)
                qualifiers[i] = feature.qualifiers

            main_df = pd.DataFrame()
            main_df['type'] = pd.Series(gene_type)
            main_df['strand'] = pd.Series(location_strand)
            main_df['start'] = pd.Series(location_start)
            main_df['end'] = pd.Series(location_end)

            main_df = main_df.assign(sequence=main_df.apply(lambda x: self.sequence[x.start:x.end], axis=1))

            df_temp = pd.DataFrame(qualifiers)
            df_temp = df_temp.transpose()

            # convert locus tag from list to str
            if self.verbose:
                df_temp['locus_tag'] = df_temp['locus_tag'].str[0]
            else:
                df_temp['gene'] = df_temp['gene'].str[0]

            # concat main and temp.
            _df = pd.concat([main_df, df_temp], axis=1)

            # Drop irrelevant columns
            if self.cols_to_drop is not None:
                _df = _df.drop(columns=self.cols_to_drop)

            self.gb_df = _df

    def __init__(self,
                 genebank_file: str = None,
                 unifile: str = None,
                 _table: dict = table,
                 hidrophobic_amino: list = Hidrophobic,
                 hidrophilic_amino: list = Hidrophilic,
                 cols_to_drop: list = None,
                 verbose: bool = True):

        self.gb_features = None
        self.uni_df = None
        self.gb_df = None
        self.sequence = None
        self.verbose = verbose
        self.genebank_file = genebank_file
        self.unifile = unifile
        self.table = table
        self.hidrophobic = hidrophobic_amino
        self.hidrophilic = hidrophilic_amino
        self.cols_to_drop = cols_to_drop
        self.read_uni_dataframe()
        self.parse_genebank_file()
        self.create_gb_dataframe()
