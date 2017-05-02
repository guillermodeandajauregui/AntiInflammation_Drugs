#################################
#Anti-inflammatory Drug Network #
#################################

################################
#04 PATHWAY GSEA
################################

#For each drug 
#GSEA (drug profile, list pws)
#express as network

source("libraries/libraries.R")
load("~/GOs_and_pathways.RData")

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


################################################################
stop()

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


test_cmap_fgsea = lapply(X = colnames(cmap_matrix_a), function(x){
  trank = cmap_matrix_a[,x]
  names(trank) = rownames(cmap_matrix_a)
  trank = trank[order(trank, decreasing = TRUE)]
  figsy =  fgsea(pathways = PWs, 
                 stats = trank, 
                 nperm = 50000, 
                 maxSize = 500)
  top_figsy = figsy[padj < 0.1]
  return(top_figsy)
}
)

names(test_cmap_fgsea) <- colnames(cmap_matrix_a)
g_fgsea = igraph::make_empty_graph(n = 0, directed = FALSE)
g_fgsea = add_vertices(graph = g_fgsea, 
                       nv = length(colnames(cmap_matrix_a)), 
                       attr = list(type = "DRUG", name = colnames(cmap_matrix_a), color = "red"))

g_fgsea = add_vertices(graph = g_fgsea, 
                       nv = length(unique(unlist(lapply(test_cmap_fgsea, function(x) x$pathway)))),
                       attr = list(type = "PATHWAY",
                                   color = "blue",
                                   name = unique(unlist(lapply(test_cmap_fgsea, function(x) x$pathway))))
                      )

for(i in colnames(cmap_matrix_a)){
  farm = test_cmap_fgsea[[i]]$pathway
  for(j in farm){
    g_fgsea<-add.edges(graph = g_fgsea, edges = c(i,j))
  }
  rm(i)
  rm(farm)
}

plot(g_fgsea, vertex.label = "")
degree(g_fgsea, v = V(g_fgsea)[type == "PATHWAY"])[which(degree(g_fgsea, v = V(g_fgsea)[type == "PATHWAY"])==max(degree(g_fgsea, v = V(g_fgsea)[type == "PATHWAY"])))]
