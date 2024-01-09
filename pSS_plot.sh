conda activate cellphonedb

cellphonedb plot dot_plot --pvalues-path /home/hlw/hlw/jia_cellphonedb/pSS/pSS_output/pvalues.txt --means-path /home/hlw/hlw/jia_cellphonedb/pSS/pSS_output/means.txt --output-path /home/hlw/hlw/jia_cellphonedb/pSS/pSS_output --output-name pSS.dotplot.pdf

cellphonedb plot heatmap_plot --pvalues-path /home/hlw/hlw/jia_cellphonedb/pSS/pSS_output/pvalues.txt --output-path /home/hlw/hlw/jia_cellphonedb/pSS/pSS_output --pvalue 0.05 --count-name pSS.heatmap_count.pdf --log-name pSS.heatmap_log_count.pdf --count-network-name pSS.count_network.txt --interaction-count-name pSS.interaction_count.txt pSS_meta.txt
