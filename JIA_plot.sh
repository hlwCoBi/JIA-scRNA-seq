conda activate cellphonedb

cellphonedb plot heatmap_plot --pvalues-path /home/hlw/hlw/jia_cellphonedb/JIA/JIA_output/pvalues.txt --output-path /home/hlw/hlw/jia_cellphonedb/JIA/JIA_output/ --pvalue 0.05 --count-name JIA.heatmap_count.pdf --log-name JIA.heatmap_log_count.pdf --count-network-name JIA.count_network.txt --interaction-count-name JIA.interaction_count.txt JIA_meta.txt
