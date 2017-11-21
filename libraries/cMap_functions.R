########################################################
#
#Functions for handling connectivity Map (cMap) data. 
# Guillermo de Anda - Jauregui
# guillermo.deandajaur@med.und.edu
#
########################################################

###########
#Libraries#
###########
library("data.table")
library("stringr")
library("GeneExpressionSignature")
library("igraph")
###########

######
# fread.Ranked.Matrix 
#Read Ranked Matrices
######

fread.Ranked.Matrix <- function(filepath){
  mx<-fread(filepath, data.table = FALSE)
  rownames(mx) = mx[,1]
  colnames(mx) = mx[1,]
  mx = mx[-1,-1]
  return(mx)
}

######
# drug.eset
# make an ExpressionSet of samples 
# perturbed by drugs on a list, 
# based on a file with sample information
######

drug.eset<- function(RankedMatrix, SampleInfoFile, DrugList){
  #Identify samples treated with drugs of interest
  DrugInstances = SampleInfoFile$instance_id[which(SampleInfoFile$normalized_name%in%DrugList)]
  #make a subset of the RankedMatrix, containing only these Ids
  Drug.RankMatrix = RankedMatrix[,as.character(DrugInstances)]
  #dummy eset, to force sample order
  eset = ExpressionSet(assayData = as.matrix(Drug.RankMatrix))
  #make drug annotation object for the ExpressionSet
  DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 4)] #id, drugname
  rownames(DrugData) = DrugData[,1]
  DrugData = DrugData[,-1, drop = FALSE]
  DrugAFD = new("AnnotatedDataFrame", data = DrugData)
  #make the ExpressionSet
  eset = ExpressionSet(assayData = as.matrix(Drug.RankMatrix), phenoData = DrugAFD)
  return(eset)
}

######
# cMap.eset
# make an ExpressionSet of samples 
# based on a file with cMap instance information
######

cMap.eset<- function(RankedMatrix, SampleInfoFile){
  #dummy eset, to force sample order
  eset = ExpressionSet(assayData = as.matrix(RankedMatrix))
  #make drug annotation object for the ExpressionSet
  cMapData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), ] 
  rownames(cMapData) = cMapData[,1]
  cMapData = DrugData[,-1, drop = FALSE]
  cMapAFD = new("AnnotatedDataFrame", data = cMapData)
  #make the ExpressionSet
  eset = ExpressionSet(assayData = as.matrix(RankedMatrix), phenoData = cMapAFD)
  return(eset)
}

######
# read.Drugs
# helper function to read and format a drug list file
# 
######

read.Drugs<- function(DrugList){
  drugList = read.table(file = DrugList, 
                        sep = "$" #a fake separator, to avoid problems with spaces in drug names
                        )
  drugList = drugList$V1
  drugList = str_to_upper(drugList)
}

#######
# Drug.Gene.Graph
# Takes a Kru-Bor merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, unweighted network from DRUGS to GENES
#######

Drug.Gene.Graph <- function(MergedDrugEset, Threshold){
  DrugGraph = exprs(MergedDrugEset)
  minz = Threshold
  maxz = max(DrugGraph) - Threshold
  DrugGraph <- ifelse(minz>=DrugGraph | DrugGraph>maxz, 1, 0)
  DrugGraph = graph_from_incidence_matrix(incidence = DrugGraph, 
                                          directed = TRUE, 
                                          mode = "in")
  return(DrugGraph)
}

######################################
# TopBottomComparer
# Compares Top and Bottom Genes
# Which are the neighbors of drugs in
# The bipartite graph
######################################

TopBottomComparer <-function(G1, G2, NodeList, listed = FALSE){
#Make sure NodeList follows a determined order; 
  #not important here,
  #but very important if you are comparing non-Graph objects
    IntersectList<- list()
  for(i in seq_along(NodeList)){
    G1_neighbors = neighbors(G1, v = NodeList[i])$name
    G2_neighbors = neighbors(G2, v = NodeList[i])$name
    q<- intersect(G1_neighbors, G2_neighbors)
    if(listed==TRUE){
      IntersectList<-append(IntersectList, list(q))
      names(IntersectList)[i] <- NodeList[i]
    }
    else{
      IntersectList<-c(IntersectList, length(q))
      names(IntersectList)[i] <- NodeList[i]
    }
  }
    return(IntersectList)
}

######################################
# TopBottomSetDiff
# Compares Top and Bottom Genes
# Which are the neighbors of drugs in
# The bipartite graph
#returns the unique genes of the first graph
######################################

TopBottomSetDiff <-function(G1, G2, NodeList, listed = FALSE){
  #Make sure NodeList follows a determined order; 
  #not important here,
  #but very important if you are comparing non-Graph objects
  IntersectList<- list()
  for(i in seq_along(NodeList)){
    G1_neighbors = neighbors(G1, v = NodeList[i])$name
    G2_neighbors = neighbors(G2, v = NodeList[i])$name
    q<- setdiff(G1_neighbors, G2_neighbors)
    if(listed==TRUE){
      IntersectList<-append(IntersectList, list(q))
      names(IntersectList)[i] <- NodeList[i]
    }
    else{
      IntersectList<-c(IntersectList, length(q))
      names(IntersectList)[i] <- NodeList[i]
    }
  }
  return(IntersectList)
}

#######
# Drug.Gene.Sign.Graph
# Takes a Kru-Bor merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, weighted network from DRUGS to GENES
#with weight being either "plus" or "minus"
#interpretable as an activation or repression of said GENE expression
#by DRUG
#######

Drug.Gene.Sign.Graph <- function(MergedDrugEset, Threshold){
  DrugGraph = exprs(MergedDrugEset)
  minz = Threshold
  maxz = max(DrugGraph) - Threshold
  
  DrugGraph<-ifelse(minz>=DrugGraph, 1, 
                    ifelse(DrugGraph>maxz, -1,
                           0)
  )
  DrugGraph = graph_from_incidence_matrix(incidence = DrugGraph, 
                                          directed = TRUE, 
                                          mode = "in",
                                          weighted = TRUE)
  return(DrugGraph)
}

#######
# Drug.Gene.Sign.Graph.for.multi
# Takes Expression matrix from 
# Kru-Bor merged eset, merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, weighted network from DRUGS to GENES
#with weight being either "plus" or "minus"
#interpretable as an activation or repression of said GENE expression
#by DRUG
#Re write of function for efficiency in Threshold evaluation
#######

Drug.Gene.Sign.Graph.for.multi <- function(ExprsMergedDrugEset, Threshold, MaxDrugraph){
  #DrugGraph = exprs(MergedDrugEset)
  #minz = Threshold
  maxz = MaxDrugraph - Threshold
  
  ExprsMergedDrugEset<-ifelse(Threshold>=ExprsMergedDrugEset, 1, 
                    ifelse(ExprsMergedDrugEset>maxz, -1,
                           0)
  )
  
  return(graph_from_incidence_matrix(incidence = ExprsMergedDrugEset, 
                                          directed = TRUE, 
                                          mode = "in",
                                          weighted = TRUE)
  )
  
}


#######
# ConnectedComponentMembership
# Takes a graph (and optionally, a components object, default to components(graph))
# and a minimum Component Size (default 1)
#returns a list with the NAMES of NODES in each component.
#######

ConnectedComponentMembership<-function(graph, 
                                       ComponentObject= components(graph), 
                                       ComponentSize = 1){
  q = lapply(seq_along(ComponentObject$csize)[ComponentObject$csize>=ComponentSize],
             function(x) V(graph)$name[q$membership%in%x]
  )
  return(q)
}

#######
# ActivatorsInhibitors
# Takes a graph (and optionally, a set of GENE vertices)
#returns a DataFrame with the number of activators and inhibitors
#For said gene.
#######

ActivatorsInhibitors <- function(DrugGraph, 
                                 nodes=V(DrugGraph)[V(DrugGraph)$type!=TRUE]){
  
  df = data.frame(gene = as.character(),
                  activators = as.numeric(),
                  inhibitors = as.numeric(),
                  stringsAsFactors = FALSE
  )
  
  for(i in seq_along(nodes)){
    gen = nodes[i]$name
    Nhoods = ego(graph = DrugGraph, 
                 order = 1, 
                 nodes = gen, 
                 mode = "in")
    Nhood_sb = induced.subgraph(graph = DrugGraph, 
                                vids = unlist(Nhoods)
    )
    activator = 0
    inhibitor = 0
    for(j in E(Nhood_sb)$weight){
      if(j==1){
        activator=activator+1
      }else
        if(j==-1){
          inhibitor=inhibitor+1
        }
    }
    k=list(gen, activator, inhibitor)
    df[i,]<-k
  }
  return(df)
}

#######
# NodeDirStrength
# takes network and node
# returns sum of edges adyacent to node
#if given more than one node, will sum edges adyacent to ALL
#probably you don't want that. Use lapply (X = NODES instead)
#######

NodeDirStrength = function(graph, node, act.or.inhib = "act"){
  hood = ego(graph = graph, 
             order = 1, 
             nodes = node, 
             mode = "in"
  )  
  
  nhood = induced.subgraph(graph = graph, 
                           vids = unlist(hood))
  if(act.or.inhib=="act"){
    a = 1
  }else
    if(act.or.inhib=="inhib"){
      a=-1
    }
  s = which(E(nhood)$weight==a)
  activation= sum(E(nhood)$weight[s])
  return(abs(activation))
}

#######
#SetNodesDirStrength
#A wrapper for the last function over a set of nodes
#if unlist = TRUE, returns vector of strengths
#else, a list
#######
SetNodesDirStrength<-function(graph, 
                              NodeSet = V(graph)[V(graph)$type==FALSE],
                                unlist = FALSE, 
                              act.or.inhib = "act"){
  if(unlist == FALSE){
    return(lapply(X = NodeSet, 
                  FUN = NodeDirStrength, 
                  act.or.inhib=act.or.inhib, 
                  graph = graph))
  }else
  return(unlist(lapply(X = NodeSet, 
                       FUN = NodeDirStrength, 
                       act.or.inhib=act.or.inhib, 
                       graph = graph)))
}


#######
# ConnectedNotes
# Takes a graph 
# Removes nodes with degree < 1
#
#######

ConnectedNodes = function(graph){
  g = graph
  V(g)$Degree = degree(g)
  gg = induced_subgraph(graph = g, v = V(g)[Degree>0])
  return(gg)
}
#######
# Drug.Gene.Sign.Graph.Connected.for.multi
# Takes Expression matrix from 
# Kru-Bor merged eset, merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, weighted network from DRUGS to GENES
#with weight being either "plus" or "minus"
#interpretable as an activation or repression of said GENE expression
#by DRUG
#Removing all GENES that do not connect to DRUGS
#Re write of function for efficiency in Threshold evaluation
#######

Drug.Gene.Sign.Graph.Connected.for.multi <- function(ExprsMergedDrugEset, Threshold, MaxDrugraph){
  #DrugGraph = exprs(MergedDrugEset)
  #minz = Threshold
  maxz = MaxDrugraph - Threshold
  
  ExprsMergedDrugEset<-ifelse(Threshold>=ExprsMergedDrugEset,  1, 
                              ifelse(ExprsMergedDrugEset>maxz, -1,
                                     0)
  )
  
  return(ConnectedNodes(graph = graph_from_incidence_matrix(incidence = ExprsMergedDrugEset, 
                                      directed = TRUE, 
                                      mode = "in",
                                      weighted = TRUE)
                        )
  )
  
}

#######
# nw_analysis_function
# Handy function to get some general network parameters
#
#######

nw_analysis_function <- function(network){
  g = ConnectedNodes(network)
  nodes = length(V(g))
  Drugs = length(V(g)[V(g)$type==TRUE])
  Genes = length(V(g)[V(g)$type==FALSE])
  Edges = length(E(g))
  cc = components(g)$no
  max.degree = max(degree(graph = g, 
                          v = V(g)[V(g)$type==FALSE]))
  max.activation = max(SetNodesDirStrength(graph = g, 
                                           NodeSet = V(g)[V(g)$type==FALSE], 
                                           unlist = TRUE, 
                                           act.or.inhib = "act")
  )
  max.inhibition = max(SetNodesDirStrength(graph = g, 
                                           NodeSet = V(g)[V(g)$type==FALSE], 
                                           unlist = TRUE, 
                                           act.or.inhib = "inhib")
  )
  results = list(nodes, Drugs, Genes, Edges, cc, max.degree, max.activation, max.inhibition)
  names(results) = c("nodes", "drugs", "genes", "edges", "cc", "degree", "act", "inhib")
  return(results)
}

# nw_analysis_function_df
# Handy function to get some general network parameters
# Returns a dataframe.
#
#######
nw_analysis_function_df<-function(network_list){
  nwan = lapply(network_list, FUN = nw_analysis_function)
  nwan2 = rbindlist(nwan)
  nwan3 = as.data.frame(nwan2)
  rownames(nwan3) = names(network_list)
  return(nwan3)
}

#######
# mc_nw_analysis_function_df
# Handy function to get some general network parameters
# Returns a dataframe.
# multicore version
#######
mc_nw_analysis_function_df<-function(network_list, cores){
  nwan = mclapply(network_list, FUN = nw_analysis_function, mc.cores = cores)
  nwan2 = rbindlist(nwan)
  nwan3 = as.data.frame(nwan2)
  rownames(nwan3) = names(network_list)
  return(nwan3)
}

#######
# nw_analysis_function2
# Handy function to get some general network parameters
# modified for graphs with unconnected nodes already removed
#
#######

nw_analysis_function2 <- function(g){
  #g = ConnectedNodes(network)
  nodes = length(V(g))
  Drugs = length(V(g)[V(g)$type==TRUE])
  Genes = length(V(g)[V(g)$type==FALSE])
  GENES = V(g)[V(g)$type==FALSE]
  Edges = length(E(g))
  cc = components(g)$no
  max.degree = max(degree(graph = g, 
                          v = GENES))
  max.activation = max(SetNodesDirStrength(graph = g, 
                                           NodeSet = GENES, 
                                           unlist = TRUE, 
                                           act.or.inhib = "act")
  )
  max.inhibition = max(SetNodesDirStrength(graph = g, 
                                           NodeSet = GENES, 
                                           unlist = TRUE, 
                                           act.or.inhib = "inhib")
  )
  results = list(nodes, Drugs, Genes, Edges, cc, max.degree, max.activation, max.inhibition)
  names(results) = c("nodes", "drugs", "genes", "edges", "cc", "degree", "act", "inhib")
  return(results)
}
#######

#######
# nw_analysis_function2_df
# Handy function to get some general network parameters
# Returns a dataframe.
#
#######
nw_analysis_function2_df<-function(network_list){
  nwan = lapply(network_list, FUN = nw_analysis_function2)
  nwan2 = rbindlist(nwan)
  nwan3 = as.data.frame(nwan2)
  rownames(nwan3) = names(network_list)
  return(nwan3)
}

#######
# mc_nw_analysis_function2_df
# Handy function to get some general network parameters
# Returns a dataframe.
# multicore version
#######
mc_nw_analysis_function2_df<-function(network_list, cores){
  nwan = mclapply(network_list, FUN = nw_analysis_function2, mc.cores = cores)
  nwan2 = rbindlist(nwan)
  nwan3 = as.data.frame(nwan2)
  rownames(nwan3) = names(network_list)
  return(nwan3)
}

#######
# ShuffleRanks
# Takes a RankedMatrix
# Shuffles rank values for each sample.
# 
#######

ShuffleRanks = function(y){
  yy = apply(y, MARGIN = 2, FUN = function(x) x[sample(x = x, size = length(x), replace = FALSE)])
  rownames(yy) <- rownames(y)
  return(yy)
}

#######
# ListShuffleRanks
# Takes a RankedMatrix
# Returns a list of n matrices with shuffled rank values
# 
#######

ListShuffleRanks = function(matrix, n=100){
  list_Shuffled = do.call(what = list, args = replicate(n = n, 
                                                   expr = ShuffleRanks(matrix), 
                                                   simplify = FALSE))  
}

#######
#Degree function
#Wrapper for ease of multiplexing
#######
degree_function<-function(graph){
  return(degree(graph = graph, v = V(graph)[V(graph)$type==FALSE]))
}

######
#degree_table_function
#takes a list of graphs
#calculates the degree of all nodes
#returns a data.frame
#with 

degree_table_function <- function(GraphList){
  return(as.data.frame(t(as.data.frame(
    lapply(
      GraphList, 
      degree_function)))))
}


#######
# NodeDirStrength2
#
#EXPERIMENTAL
#Better implementation
#
# takes network and node
# returns sum of edges adyacent to node
#if given more than one node, will sum edges adyacent to ALL
#probably you don't want that. Use lapply (X = NODES instead)
#######

NodeDirStrength2 = function(graph, node, act.or.inhib = "act"){
  if(act.or.inhib=="act"){
    a = 1
  }else
    if(act.or.inhib=="inhib"){
      a=-1
    }
  return(sum(E(graph)[to(node)]$weight[which((E(graph)[to(node)]$weight == a))]))
}

#######
#Act_OR_Inhib_Degree
#Variation of the function, in which a subgraph of 
#######
Act_OR_Inhib_Degree<-function(graph, act.or.inhib = "act"){
  if(act.or.inhib == "act"){
    g = subgraph.edges(graph,
                       eids = E(graph)[E(graph)$weight==1],
                       delete.vertices = FALSE)
    return(degree(graph = g, v = V(g)[V(g)$type==FALSE])
           )
  }else if(act.or.inhib == "inhib"){
    g = subgraph.edges(graph,
                       eids = E(graph)[E(graph)$weight==-1],
                       delete.vertices = FALSE)
    return(degree(graph = g, v = V(g)[V(g)$type==FALSE])
                )
            }
}

######
#signed_degree_table_function
#takes a list of graphs
#calculates the Activation or Inhibition degree of all nodes
#returns a data.frame
#with 

signed_degree_table_function <- function(GraphList, act.or.inhib){
  return(as.data.frame(t(as.data.frame(
    lapply(
      GraphList, 
      Act_OR_Inhib_Degree, 
      act.or.inhib = act.or.inhib)))))
}  

######
#NodeTyoe
#writes a named vector of DRUG or GENES 
# names = node name 
#with 


NodeType = function(G){
  nodes = V(G)
  VGType = nodes$type
  q = ifelse(test =  VGType == TRUE, "DRUG", "GENE")
  names(q) = nodes$name
  return(q)
}

######
#Projection_Function
#Takes a bipartite DRUG-GENE graph
#projects to a unipartite, weighted DRUG graph
#where two DRUGS are linked if they both UP or DOWN regulate the same GENE
#Each gene shared is +1 in their link weight
######


Projection_Function<-function(graph){
  
  players<-V(graph)[V(graph)$type==TRUE]
  Dg = make_empty_graph(directed = FALSE)
  Dg = add_vertices(graph = Dg, nv = length(players), name = names(players))
  plot(Dg)
  
  for(i in seq_along(players)){
    for(j in seq_along(players)){
      if(i != j){
        if(j > i){
          genes = intersect(neighbors(graph = graph, v = players[i], mode = "all")$name,
                            neighbors(graph = graph, v = players[j], mode = "all")$name)
          peso = 0
          print(genes)
          for(k in genes){
            if(E(graph)[V(graph)[players[i]]%--%k]$weight == E(graph)[V(graph)[players[i]]%--%k]$weight){
              peso = peso + 1
            }
            
          }
          if(peso!=0){
            Dg<-add.edges(graph = Dg,
                          edges = c(V(Dg)[V(Dg)$name==players[i]$name],
                                    V(Dg)[V(Dg)$name==players[j]$name]),
                          weight = peso
            )
            
          }
        }
      }
    }
  }
  return(Dg)
}

######
#Gene_Projection_Function
#Takes a bipartite DRUG-GENE graph
#projects to a unipartite, weighted GENE graph
#where two GENE are linked if they are both UP or DOWN regulated by the same DRUG
#Each gene shared is +1 in their link weight
######


Gene_Projection_Function<-function(graph){
  
  players<-V(graph)[V(graph)$type==FALSE]
  Dg = make_empty_graph(directed = FALSE)
  Dg = add_vertices(graph = Dg, nv = length(players), name = names(players))
  plot(Dg)
  
  for(i in seq_along(players)){
    for(j in seq_along(players)){
      if(i != j){
        if(j > i){
          drugs = intersect(neighbors(graph = graph, v = players[i], mode = "all")$name,
                            neighbors(graph = graph, v = players[j], mode = "all")$name)
          peso = 0
          print(drugs)
          for(k in drugs){
            if(E(graph)[V(graph)[players[i]]%--%k]$weight == E(graph)[V(graph)[players[i]]%--%k]$weight){
              peso = peso + 1
            }
            
          }
          if(peso!=0){
            Dg<-add.edges(graph = Dg,
                          edges = c(V(Dg)[V(Dg)$name==players[i]$name],
                                    V(Dg)[V(Dg)$name==players[j]$name]),
                          weight = peso
            )
            
          }
        }
      }
    }
  }
  return(Dg)
}