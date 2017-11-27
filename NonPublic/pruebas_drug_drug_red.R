cmap_spearman = cor(as.matrix(cmap_matrix_a), method = "spearman")
lincs1_spearman = cor(as.matrix(lincs1_matrix_a), method = "spearman")
lincs2_spearman = cor(as.matrix(lincs2_matrix_a), method = "spearman")

cmap_opt_nw = funcion_corte_optimo(cmap_spearman, resolucion = 0.01)
lincs1_opt_nw = funcion_corte_optimo(lincs1_spearman, resolucion = 0.01)
lincs2_opt_nw = funcion_corte_optimo(lincs2_spearman, resolucion = 0.01)

cmap_opt_nw_analysis = NetworkAnalyzer(cmap_opt_nw$graph)
lincs1_opt_nw_analysis = NetworkAnalyzer(lincs1_opt_nw$graph)
lincs2_opt_nw_analysis = NetworkAnalyzer(lincs2_opt_nw$graph)

cmap_opt_nw_analysis$infomap$membership
lincs1_opt_nw_analysis$infomap$membership
lincs2_opt_nw_analysis$infomap$membership

plot_nicely(cmap_opt_nw_analysis)
plot_nicely(lincs1_opt_nw_analysis)
plot_nicely(lincs2_opt_nw_analysis)

cmap_opt_nw_analysis$infomap_star <- infomap.community(cmap_opt_nw_analysis$g, 
                                                        e.weights = exp(E(cmap_opt_nw_analysis$g)$weight)
                                                        )
lincs2_opt_nw_analysis$infomap_star$membership
plot_nicely_star(cmap_opt_nw_analysis)
#############
cmap_spearman = cor(as.matrix(cmap_matrix_a), method = "spearman")
plot(density((cmap_spearman)))
quantile(x = cmap_spearman, 0.95)
cmap_spearman_matrix_quantiled = ifelse(test = cmap_spearman < quantile(cmap_spearman, 0.60), 
                                        0, 
                                        1)
g = graph_from_adjacency_matrix(adjmatrix = cmap_spearman_matrix_quantiled, 
                                diag = FALSE, 
                                mode = "undirected", 
                                weighted = TRUE)

plot(g)
E(g)
V(g)
components(g)
info_g = infomap.community(g)
info_g$membership

str(exprs(ai_lincs1_krubor))

lincs1_spearman = cor(exprs(ai_lincs1_krubor), method = "spearman")
plot(density((lincs1_spearman)))
#quantile(x = cmap_spearman, 0.95)
lincs1_spearman_matrix_quantiled = ifelse(test = lincs1_spearman < quantile(lincs1_spearman, 0.87), 
                                        0, #1
                                        1/lincs1_spearman)
g = graph_from_adjacency_matrix(adjmatrix = lincs1_spearman_matrix_quantiled, 
                                diag = FALSE, 
                                mode = "undirected", 
                                weighted = TRUE)
plot(g)
E(g)$weight
V(g)
gaga = components(g)
info_g = infomap.community(g)
info_g$membership

lincs2_spearman = cor(exprs(ai_lincs2_krubor), method = "spearman")
##########################################################################
eset_cmap_matrix_a = ExpressionSet(assayData = as.matrix(cmap_matrix_a), )



ScoreGSEA(eset_cmap_matrix_a, SignatureLength = 100)

xxxx = ScoreGSEA(MergingSet = ai_cmap_krubor, SignatureLength = 1000, ScoringDistance = "avg")

max(xxxx)
plot(density(xxxx))
xxxx

library(igraph)

gx = graph_from_adjacency_matrix(adjmatrix = xxxx, mode = "undirected", weighted = TRUE)




cmap_spearman_matrix_quantiled = ifelse(test = cmap_spearman_matrix > quantile(cmap_spearman_matrix, 0.33), 
                                        0, 
                                        1)
g = graph_from_adjacency_matrix(adjmatrix = cmap_spearman_matrix_quantiled, 
                                diag = FALSE, 
                                mode = "undirected", 
                                weighted = TRUE)

plot(g)
components(g)
info_g = infomap.community(g)
info_g$membership
