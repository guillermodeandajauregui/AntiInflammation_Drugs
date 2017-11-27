#1)take an ranked expression profile matrix
##and a list of pathways
##return a list of enriched pathways through FGSEA

#2)take a list of enriched pathways through FGSEA
##return a sample(ej, drug) - enriched pathway bipartite network

#1)take an ranked expression profile matrix
##and a list of pathways
##return a list of enriched pathways through FGSEA
fgsea_expmatrix <-function(expmatrix, 
                           pathway_list, 
                           nperm = 50000, 
                           maxSize = 500, 
                           padj_cutoff = 0.1)
{
  expmatrix_fgsea = lapply(X = colnames(expmatrix), function(x)
  {
    trank = expmatrix[,x]
    names(trank) = rownames(expmatrix)
    trank = trank[order(trank, decreasing = TRUE)]
    figsy =  fgsea(pathways = pathway_list, 
                   stats = trank, 
                   nperm = nperm, 
                   maxSize = maxSize)
    top_figsy = figsy[padj < padj_cutoff]
    return(top_figsy)
  }
  )
  names(expmatrix_fgsea) <- colnames(expmatrix)
  return(expmatrix_fgsea)
}

#2)take a list of enriched pathways through FGSEA
##return a sample(ej, drug) - enriched pathway bipartite network

graph_from_fgsea<-function(fgsea_result){
  g_fgsea = igraph::make_empty_graph(n = 0, directed = FALSE)
  g_fgsea = add_vertices(graph = g_fgsea, 
                         nv = length(names(fgsea_result)), 
                         attr = list(type = "DRUG", 
                                     name = names(fgsea_result), 
                                     color = "red")
  )
  
  g_fgsea = add_vertices(graph = g_fgsea, 
                         nv = length(unique(unlist(lapply(fgsea_result, function(x) x$pathway)))),
                         attr = list(type = "PATHWAY",
                                     color = "blue",
                                     name = unique(unlist(lapply(fgsea_result, function(x) x$pathway)))
                         )
  )
  
  for(i in names(fgsea_result)){
    enriched_pws = fgsea_result[[i]]$pathway
    for(j in enriched_pws){
      g_fgsea<-add.edges(graph = g_fgsea, edges = c(i,j))
    }
  }
  
  return(g_fgsea)
}