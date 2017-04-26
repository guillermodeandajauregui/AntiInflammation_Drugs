#take a correlation matrix 
#make a network such as correlation value >= cutoff --> link
#identify percolation cutoff

funcion_corte_optimo = function(matrix, 
                                starting_quantile = 1, 
                                resolucion = 0.05){
  quantil = starting_quantile
  num_components = length(matrix) #force a non-1 number
  
  while(num_components > 1){
    
    matrix_q = ifelse(test = matrix < quantile(matrix, quantil), 
                                              0, 
                                              1/matrix)
    g = graph_from_adjacency_matrix(adjmatrix = matrix_q, 
                                    diag = FALSE, 
                                    mode = "undirected", 
                                    weighted = TRUE)
    
    num_components = components(g)$no
    running_quantil = quantil
    print(running_quantil)
    quantil = quantil - resolucion
    
  }
  return(list(graph = g, quantil = running_quantil))
}

prueba = funcion_corte_optimo(matrix = cmap_spearman, resolucion = 0.01)
E(prueba$graph)
prueba = funcion_corte_optimo(matrix = lincs1_spearman, resolucion = 0.01)
E(prueba$graph)
prueba = funcion_corte_optimo(matrix = lincs2_spearman, resolucion = 0.01)
E(prueba$graph)
infomap.community(prueba$graph)
plot(prueba$graph)

