library(igraph)
source("libraries/libraries.R")

list.files(path = "results/")

#fgsea networks 
load(file = "results/fgsea_nws.RData")
#drug drug networks
load(file = "results/drug_drug_networks_correlation_and_footrule_analysis_cmap_lincs_separated.RData")
#drug ADR networks
drug_adr = igraph::read.graph(file = "results/drug_adr.gml", format = "gml")

#merging
##id DRUGS 
DRUGS = names(V(cmap_correlation_nw_analysis$g))
V(cmap_correlation_nw_analysis$g)$NodeType = "DRUG"

V(nw_fgsea_cmap)$NodeType = ifelse(names(V(nw_fgsea_cmap))%in%DRUGS, "DRUG", "PATHWAY")

V(drug_adr)$NodeType = ifelse(names(V(drug_adr))%in%DRUGS, "DRUG", "ADR")

multigraph_cmap = igraph::union(cmap_correlation_nw_analysis$g, 
                                NetworkAnalyzer(nw_fgsea_cmap)$g, 
                                as.undirected(graph =  drug_adr, 
                                              mode  =  "collapse", 
                                              edge.attr.comb = "max")
                                )

write.graph(graph = multigraph_cmap, file = "results/multigraph_cmap_01.gml", format = "gml")
write.graph(graph = multigraph_cmap, file = "results/multigraph_cmap_01.graphml", format = "graphml")

