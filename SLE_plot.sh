conda activate cellphonedb

cellphonedb plot dot_plot --pvalues-path /home/hlw/hlw/jia_cellphonedb/SLE/SLE_output/pvalues.txt --means-path /home/hlw/hlw/jia_cellphonedb/SLE/SLE_output/means.txt --output-path /home/hlw/hlw/jia_cellphonedb/SLE/SLE_output --output-name SLE.dotplot.pdf

cellphonedb plot heatmap_plot --pvalues-path /home/hlw/hlw/jia_cellphonedb/SLE/SLE_output/pvalues.txt --output-path /home/hlw/hlw/jia_cellphonedb/SLE/SLE_output --pvalue 0.05 --count-name SLE.heatmap_count.pdf --log-name SLE.heatmap_log_count.pdf --count-network-name SLE.count_network.txt --interaction-count-name SLE.interaction_count.txt SLE_meta.txt
