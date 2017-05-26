################################
#03 Drug-Drug Networks
################################

#read Ranked Expression Matrix
#for each drug
#Spearman correlation to each other
#return spearman correlation matrix
source("libraries/libraries.R")

#read ranked matrices 

cmap_matrix= fread(input = "results/krubor_cmap.txt", data.table = FALSE)
rownames(cmap_matrix) <- cmap_matrix[,1]
cmap_matrix <- cmap_matrix[,-1]

lincs1_matrix= fread(input = "results/krubor_lincs1.txt", data.table = FALSE)
rownames(lincs1_matrix) <- lincs1_matrix[,1]
lincs1_matrix <- lincs1_matrix[,-1]

lincs2_matrix= fread(input = "results/krubor_lincs2.txt", data.table = FALSE)
rownames(lincs2_matrix) <- lincs2_matrix[,1]
lincs2_matrix <- lincs2_matrix[,-1]

#correlation matrices

cmap_spearman_correlation = cor(as.matrix(cmap_matrix), method = "spearman")
lincs1_spearman_correlation = cor(as.matrix(lincs1_matrix), method = "spearman")
lincs2_spearman_correlation = cor(as.matrix(lincs2_matrix), method = "spearman")

fwrite(x = as.data.frame(cmap_spearman_correlation), 
       file = "results/cmap_spearman_correlation_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs1_spearman_correlation), 
       file = "results/lincs1_spearman_correlation_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs2_spearman_correlation), 
       file = "results/lincs2_spearman_correlation_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

#Spearman footrule

#read Ranked Expression Matrix
#for each drug
#Spearman footrule to each other
#return spearman matrix

cmap_spearman_footrule_matrix   = spearman_footrule_matrix(cmap_matrix)
lincs1_spearman_footrule_matrix = spearman_footrule_matrix(lincs1_matrix)
lincs2_spearman_footrule_matrix = spearman_footrule_matrix(lincs2_matrix)

fwrite(x = as.data.frame(cmap_spearman_footrule_matrix), 
       file = "results/cmap_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs1_spearman_footrule_matrix), 
       file = "results/lincs1_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs2_spearman_footrule_matrix), 
       file = "results/lincs2_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

## matrices to networks using the "optimal cut" functions
## objective: keep a single component of maximum modularity 

#correlation networks
cmap_correlation_nw     = funcion_corte_optimo(matrix = cmap_spearman_correlation, 
                       starting_quantile = 1, 
                       resolucion = 0.01)

lincs1_correlation_nw     = funcion_corte_optimo(matrix = lincs1_spearman_correlation, 
                                            starting_quantile = 1, 
                                            resolucion = 0.01)

lincs2_correlation_nw     = funcion_corte_optimo(matrix = lincs2_spearman_correlation, 
                                              starting_quantile = 1, 
                                              resolucion = 0.01)

cmap_correlation_nw_analysis = NetworkAnalyzer(cmap_correlation_nw$graph)
lincs1_correlation_nw_analysis = NetworkAnalyzer(lincs1_correlation_nw$graph)
lincs2_correlation_nw_analysis = NetworkAnalyzer(lincs2_correlation_nw$graph)

#footrule networks
cmap_footrule_nw     = funcion_corte_optimo_dist(matrix = cmap_spearman_footrule_matrix, 
                                            starting_quantile = 1, 
                                            resolucion = 0.01)

lincs1_footrule_nw     = funcion_corte_optimo_dist(matrix = lincs1_spearman_footrule_matrix, 
                                              starting_quantile = 1, 
                                              resolucion = 0.01)

lincs2_footrule_nw     = funcion_corte_optimo_dist(matrix = lincs2_spearman_footrule_matrix, 
                                              starting_quantile = 1, 
                                              resolucion = 0.01)

cmap_footrule_nw_analysis = NetworkAnalyzer(cmap_footrule_nw$graph)
lincs1_footrule_nw_analysis = NetworkAnalyzer(lincs1_footrule_nw$graph)
lincs2_footrule_nw_analysis = NetworkAnalyzer(lincs2_footrule_nw$graph)

save(cmap_correlation_nw_analysis, cmap_footrule_nw_analysis, 
     lincs1_correlation_nw_analysis, lincs1_correlation_nw_analysis, 
     lincs2_correlation_nw_analysis, lincs2_correlation_nw_analysis,
     file = "results/drug_drug_networks_correlation_and_footrule_analysis.RData")
