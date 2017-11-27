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

combined_krubor_matrix= fread(input = "results/krubor_cmap.txt", data.table = FALSE)
combined_spearman_correlation = cor(as.matrix(combined_krubor_matrix), method = "spearman")
combined_spearman_footrule_matrix   = spearman_footrule_matrix(combined_krubor_matrix)
#plot(density(combined_spearman_footrule_matrix))

#networks using optimal cut function
combined_correlation_nw     = funcion_corte_optimo(matrix = combined_spearman_correlation, 
                                               starting_quantile = 1, 
                                               resolucion = 0.01)

combined_spearman_footrule_nw     = funcion_corte_optimo(matrix = combined_spearman_footrule_matrix, 
                                                   starting_quantile = 1, 
                                                   resolucion = 0.01)

#network analysis
analysis_combined_correlation_nw = NetworkAnalyzer(combined_correlation_nw$graph)
analysis_combined_spearman_footrule_nw = NetworkAnalyzer(combined_spearman_footrule_nw$graph)

save(analysis_combined_correlation_nw, analysis_combined_spearman_footrule_nw,
  file = "results/drug_drug_networks_correlation_and_footrule_analysis.RData")

plot_degree(analyisis_combined_spearman_footrule_nw$g)
plot_nicely(analysis_combined_spearman_footrule_nw)
