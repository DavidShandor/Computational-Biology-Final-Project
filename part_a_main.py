from pprint import pprint

from part_a import GatherInfoAboutAGenome


gene_info = GatherInfoAboutAGenome('Computational-Biology-Final-Project\BS168.gb',
                                   'part_a_automated.txt')
all_genes = gene_info.get_all_genes_type_and_amount()
print(f'All genes: {all_genes}')

protein_info, non_protein_info, protein_lengths, non_protein_lengths, all_genes_lengths = gene_info.characterization_of_gene_lengths()
print(f'protein info: {protein_info}\n'
      f'non protein info: {non_protein_info}\n'
      f'protein: {protein_lengths}\n'
      f'non protein: {non_protein_lengths}\n'
      f'all: {all_genes_lengths}')

gene_info.build_histograms('Non Protein',
                           non_protein_lengths)

gene_info.build_histograms('Protein',
                           protein_lengths)

gene_info.build_histograms('All',
                           all_genes_lengths)

gc_percentage_full_gene, gc_average_percentage_proteins, proteins_gc_percentage, five_highest_gc_genes, five_lowest_gc_genes = gene_info.calculate_gc_percentage_in_genes(protein_info['length'])
print(f'gc_full: {gc_percentage_full_gene}')
print(f'gc_average: {gc_average_percentage_proteins}')
print(f'five highest: {five_highest_gc_genes}')
pprint(f'five lowest: {five_lowest_gc_genes}')

gene_info.build_histograms(histogram_title='GC percentage proteins',
                           histogram_value=proteins_gc_percentage,
                           x_label='GC percentage',
                           y_label='Number of proteins',
                           x_low_lim=0,
                           x_high_lim=1,
                           y_low_lim=0,
                           y_high_lim=100)

data_file_errors = gene_info.consistent_checks_data_file()
print(f'Data file errors: {data_file_errors}')

gene_info.create_dataframe_for_errors_in_data(columns=['Gene name', 'Start position', 'End position', 'Strand', 'Error'],
                                              data=data_file_errors)