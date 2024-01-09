import cellphonedb
cpdb_file_path = './cellphonedb-data-master/cellphonedb.zip'
test_meta_file_path = './chc_meta.txt'
test_counts_file_path = './chc_count.txt'

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_file_path,
  meta_file_path = test_meta_file_path,
  counts_file_path = test_counts_file_path,
  counts_data = 'hgnc_symbol',#注意[ensembl | gene_name | hgnc_symbol] Type of gene
  iterations = 100,
  threshold = 0.1,
  threads = 1,
  output_suffix = "cHC",
  output_path = "./chc_out_path")