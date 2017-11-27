################################
#03 Drug-Drug Networks
# Combined CMap - LINCS branch
################################

#read Ranked Expression Matrix
#for each drug
#Spearman correlation to each other
#return spearman correlation matrix
source("libraries/libraries.R")

#read ranked matrices 

combined_krubor_matrix= fread(input = "results/krubor_combined_lincs_1_and_2_ai.txt", data.table = FALSE)
rownames(combined_krubor_matrix) = combined_krubor_matrix$V1
combined_krubor_matrix = combined_krubor_matrix[,-1]
combined_spearman_correlation = cor(as.matrix(combined_krubor_matrix), method = "spearman")
combined_spearman_footrule_matrix   = spearman_footrule_matrix(combined_krubor_matrix)


#write correlation matrices
fwrite(x = as.data.frame(combined_spearman_correlation), 
       file = "results/lincs_1_and_2_spearman_correlation_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(combined_spearman_footrule_matrix), 
       file = "results/lincs_1_and_2_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)
#plot(density(combined_spearman_footrule_matrix))

#networks using optimal cut function
combined_correlation_nw     = funcion_corte_optimo(matrix = combined_spearman_correlation, 
                                               starting_quantile = 1, 
                                               resolucion = 0.01) #0.77

combined_spearman_footrule_nw     = funcion_corte_optimo(matrix = combined_spearman_footrule_matrix, 
                                                   starting_quantile = 1, 
                                                   resolucion = 0.01) #0.84



#network analysis
analysis_LINCS_combined_correlation_nw = NetworkAnalyzer(combined_correlation_nw$graph)
analysis_LINCS_combined_spearman_footrule_nw = NetworkAnalyzer(combined_spearman_footrule_nw$graph)

save(analysis_LINCS_combined_correlation_nw, analysis_LINCS_combined_spearman_footrule_nw,
  file = "results/drug_drug_networks_correlation_and_footrule_analysis_combn_LINCS.RData")

plot_degree(analyisis_combined_spearman_footrule_nw$g)
plot_nicely(analysis_combined_spearman_footrule_nw)
