library(igraph)

NetworkAnalyzer = function(g, directed = FALSE){
  if(directed == FALSE){
    V(g)$degree <- degree(g)
    V(g)$betweenness <- betweenness(g, directed = FALSE)
    V(g)$closeness <- closeness(g)
    
    average_path_length = average.path.length(g, directed = FALSE)
    ClusteringCoefficient = transitivity(g, type = "global")
    GraphDiameter = diameter(g, directed = FALSE)
    
    component = components(g)
    V(g)$component = component$membership
    
    infomap = infomap.community(g)
    V(g)$infomap = infomap$membership
    
    results = list(g = g, 
                   average_path_length = average_path_length, 
                   ClusteringCoefficient = ClusteringCoefficient, 
                   GraphDiameter = GraphDiameter,
                   component = component,
                   infomap = infomap)
  }
}

###########################
#Edge Jaccard 
#Takes a list of graphs
#returns a Jaccard matrix
#of similarity in edge sets
###########################
library(igraph)

jaccard_edge<-function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

jaccard_edge_matrix <- function(a, b=a){
  matriz = matrix(nrow = length(a), ncol = length(b), dimnames = list(names(a), names(b)))
  for(i in seq_along(a)){
    for(j in seq_along(b)){
      matriz[i,j]<-jaccard_edge(a[[i]], b[[j]])
    }
  }
  return(matriz)
}

jaccard<-function(a,b)
{
  x<-intersect(a,b)
  y<-union(a,b)
  nx<-length(x)
  ny<-length(y)
  J<-as.numeric(nx/ny)
  print(J)
  return(J)
}

#jaccard_matrix: calculate Jaccard index between two lists of sets
#returns a Jaccard matrix 
#with J between each set in list 
jaccard_matrix<-function(x,y){
  matriz<-matrix(nrow = length(x), ncol = length(y))  
  colnames(matriz)<-names(y)
  rownames(matriz)<-names(x)
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      alfa<-x[[i]]
      beta<-y[[j]]
      matriz[i,j]<-jaccard(alfa,beta)
    }
  }
  #print(matriz)
  return(matriz)  
}

#function enrichment 
function_enrichment_pathway_infomap = function(nw_an, lista_pws){
  list_enriched_pws = list()
  for(i in 1:max(nw_an$infomap$membership)){
    genes_isla = names(V(nw_an$g))[which(V(nw_an$g)$infomap==i)]
    h8 = multiHyperGeoTest(collectionOfGeneSets = lista_pws,
                           universe = names(V(nw_an$g)), #genes IN network
                           hits = genes_isla, # genes in island 
                           minGeneSetSize = 10, #all pathways larger than 10
                           pAdjustMethod = "BH", 
                           verbose = FALSE)
    h8 = as.data.frame(h8)
    h8 = h8[h8$Adjusted.Pvalue < 0.1,]
    list_enriched_pws = c(list_enriched_pws, list(h8))
  }
  names(list_enriched_pws)<-1:max(nw_an$infomap$membership)
  return(list_enriched_pws)
}

plot_nicely<-function(Annie){#result of Network Analyzer
  paleta = rainbow(n = max(Annie$infomap$membership))
  V(Annie$g)$color <- paleta[Annie$infomap$membership]
  plot(Annie$g, vertex.label = "", vertex.size = 4, edge.width = 0.5)
}

plot_nicely_star<-function(Annie){#result of Network Analyzer with the extra "infomap_star" (weight = exp(weight))
  paleta = rainbow(n = max(Annie$infomap_star$membership))
  V(Annie$g)$color <- paleta[Annie$infomap_star$membership]
  plot(Annie$g, vertex.label = "", vertex.size = 4, edge.width = 0.5)
}

plot_degree <- function(g){
  degz = table(degree(g))
  degz = degz[order(degz)]
  plot(x = names(degz), y = degz, main = "degree distribution", type = "p",
       xlab = "degree", ylab = "frequency")
}

