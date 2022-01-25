from pprint import pprint

from parts_functions import GatherInfoAboutAGenome

gene_info_about_a_genome = GatherInfoAboutAGenome('BS168.gb',
                                                  'part_a_automated.txt')
print('*--------------- Part A ---------------*')
print('Question 1')
all_genes = gene_info_about_a_genome.get_all_genes_type_and_amount()
print(f'All genes and their count in file: {all_genes}')

print('\nQuestion 2')
protein_info, non_protein_info, protein_lengths, non_protein_lengths = gene_info_about_a_genome.characterization_of_gene_lengths()


print(f'Statistics of each group:')
print(f'proteins: {protein_info}\n'
      f'non proteins: {non_protein_info}\n')

print(f'The histograms are shown on screen')

gene_info_about_a_genome.build_histograms(hist_title='Non Protein length',
                                          hist_value=non_protein_lengths,
                                          x_label='Length',
                                          y_label='Number of genes',
                                          x_low_lim=0,
                                          x_high_lim=4000,
                                          y_low_lim=0,
                                          y_high_lim=200,
                                          bins_num=300,
                                          show_grid=True,
                                          graph_color='blue')

gene_info_about_a_genome.build_histograms(hist_title='Non Protein length',
                                          hist_value=protein_lengths,
                                          x_label='Length',
                                          y_label='Number of genes',
                                          x_low_lim=0,
                                          x_high_lim=4000,
                                          y_low_lim=0,
                                          y_high_lim=200,
                                          bins_num=300,
                                          show_grid=True,
                                          graph_color='green')

gene_info_about_a_genome.build_histograms(hist_title='All genes length',
                                          hist_value=protein_lengths,
                                          x_label='Length',
                                          y_label='Number of genes',
                                          x_low_lim=0,
                                          x_high_lim=4000,
                                          y_low_lim=0,
                                          y_high_lim=200,
                                          bins_num=300,
                                          show_grid=True,
                                          graph_color='yellow')

print('\nQuestion 3')
gc_percentage_full_gene, gc_average_percentage_proteins, proteins_gc_percentage, five_highest_gc_genes, five_lowest_gc_genes, all_genes_data = gene_info_about_a_genome.calculate_gc_percentage_in_genes(protein_info['total_length'])
print(f'GC percentage of the whole genome: {gc_percentage_full_gene}')
print(f'GC average percentage in the proteins: {gc_average_percentage_proteins}')

print(f'The histogram is shown on screen')
gene_info_about_a_genome.build_histograms(hist_title='GC percentage proteins',
                                          hist_value=proteins_gc_percentage,
                                          x_label='GC percentage',
                                          y_label='Number of proteins',
                                          x_low_lim=0,
                                          x_high_lim=1,
                                          y_low_lim=0,
                                          y_high_lim=100,
                                          bins_num=300,
                                          show_grid=True,
                                          graph_color='blue')

print(f'Five highest GC percentage genes: {five_highest_gc_genes}')
print(f'Five lowest GC percentage genes: {five_lowest_gc_genes}')

print('\nQuestion 4')
data_file_errors = gene_info_about_a_genome.consistent_checks_data_file()
print(f'Data file errors: {data_file_errors}')

gene_info_about_a_genome.create_csv_file(columns=['Gene name', 'Start position', 'End position', 'Strand', 'Error'],
                                         data=data_file_errors,
                                         csv_file_name='gene_exceptions.csv')

gene_info_about_a_genome.create_csv_file(columns=['GC percentage', 'Gene name', 'Start position', 'End position', 'Strand'],
                                         data=all_genes_data,
                                         csv_file_name='part_a.csv')