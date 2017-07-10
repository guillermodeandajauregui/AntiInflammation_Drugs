##################
#hlt_pt_drug graph
##################

#auxiliar graph
#links MedDRA higher level terms to drugs
#through ADRs identified as significantly (PRR >= 1) associated to drug

library(igraph)
library(data.table)
#read medra hlt_pt relationships
hlt_pt = fread("medDRA_hlt_pt.txt")
#transform to graph
g_hlt_pt = graph_from_incidence_matrix(incidence = as.matrix(table(hlt_pt)), directed = FALSE)
#find pts in the drug_adr network
a = V(g_hlt_pt)$name[V(g_hlt_pt)$name%in%V(drug_adr)$name]
#extract subgraph of first neighbors
subg2 = induced.subgraph(graph=g_hlt_pt,
                         vids=unlist(neighborhood(graph=g_hlt_pt,
                                                  order=1,
                                                  nodes=a)))

#merge with drug-adr graph
g_hlt_pt_drug = igraph::union(subg2, as.undirected(drug_adr))

#add types
V(g_hlt_pt_drug)$type = "HLT"
V(g_hlt_pt_drug)[V(g_hlt_pt_drug)$name%in%V(drug_adr)$name[which(V(drug_adr)$type=="ADR")]]$type = "ADR"
V(g_hlt_pt_drug)[V(g_hlt_pt_drug)$name%in%V(drug_adr)$name[which(V(drug_adr)$type=="DRUG")]]$type = "DRUG"

#projection HLT - Drug

##get incidence matrices 
incidence_hlt_pt = get.incidence(subg2)
#incidence_hlt_pt[1:5, 1:5]
#nrow(incidence_hlt_pt)
#ncol(incidence_hlt_pt)
incidence_hlt_pt = t(incidence_hlt_pt) #so pts are in rows
incidence_hlt_pt = incidence_hlt_pt[order(rownames(incidence_hlt_pt)),] #so they have alfabetical order
#incidence_hlt_pt[1:5, 1:5]

drug_adr_2 = drug_adr
V(drug_adr_2)$type = bipartite.mapping(drug_adr_2)$type
incidence_drug_adr = get.incidence(drug_adr_2, attr = "PRR")
incidence_drug_adr = incidence_drug_adr[,order(colnames(incidence_drug_adr))]
#nrow(incidence_drug_adr)
#ncol(incidence_drug_adr)
#incidence_drug_adr[1:5, 1:5]
#table(colnames(incidence_drug_adr)==rownames(incidence_hlt_pt))

##project to hlt x rows 
projection_matrix = incidence_drug_adr%*%incidence_hlt_pt
projection_matrix[1:5, 1:5]

#make network
g_hlt_drug = graph_from_incidence_matrix(projection_matrix,  directed = FALSE, weighted = TRUE)
g_hlt_drug
igraph::edge_density(g_hlt_drug)
plot(g_hlt_drug)
components(g_hlt_drug)
degree(g_hlt_drug, v = V(g_hlt_drug)[type == FALSE])
tail(degree(g_hlt_drug, v = V(g_hlt_drug)[type == TRUE])[order(degree(g_hlt_drug, v = V(g_hlt_drug)[type == TRUE]))])

write.graph(g_hlt_drug, file = "results/drug_hlt.gml", format = "gml")

