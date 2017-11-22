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

#read matrices 
cmap_matrix= fread(input = "results/krubor_cmap_ai.txt", data.table = FALSE)
rownames(cmap_matrix) <- cmap_matrix[,1]
cmap_matrix <- cmap_matrix[,-1]

lincs1_matrix= fread(input = "results/krubor_lincs1_ai.txt", data.table = FALSE)
rownames(lincs1_matrix) <- lincs1_matrix[,1]
lincs1_matrix <- lincs1_matrix[,-1]

lincs2_matrix= fread(input = "results/krubor_lincs2_ai.txt", data.table = FALSE)
rownames(lincs2_matrix) <- lincs2_matrix[,1]
lincs2_matrix <- lincs2_matrix[,-1]

#fgsea enrichment 

PWs_all = PWs
PWs <- PWs[which(sapply(PWs, length)>10)]

fgsea_cmap = fgsea_expmatrix(expmatrix = cmap_matrix, 
                             pathway_list = PWs, 
                             nperm = 50000, 
                             maxSize = 500, 
                             padj_cutoff = 0.1)

fgsea_lincs1 = fgsea_expmatrix(expmatrix = lincs1_matrix, 
                             pathway_list = PWs, 
                             nperm = 50000, 
                             maxSize = 500, 
                             padj_cutoff = 0.1)

fgsea_lincs2 = fgsea_expmatrix(expmatrix = lincs2_matrix, 
                               pathway_list = PWs, 
                               nperm = 50000, 
                               maxSize = 500, 
                               padj_cutoff = 0.1)

#generate networks

nw_fgsea_cmap = graph_from_fgsea(fgsea_cmap)
nw_fgsea_lincs1 = graph_from_fgsea(fgsea_lincs1)
nw_fgsea_lincs2 = graph_from_fgsea(fgsea_lincs2)

save(fgsea_cmap, fgsea_lincs1, fgsea_lincs2, nw_fgsea_cmap, nw_fgsea_lincs1, nw_fgsea_lincs2,
     file = "results/fgsea_nws.RData")

#
plot(nw_fgsea_cmap, vertex.label = "")
plot(nw_fgsea_lincs1, vertex.label = "")
plot(nw_fgsea_lincs2, vertex.label = "")

nw_intersection = igraph::intersection(nw_fgsea_cmap, nw_fgsea_lincs1, nw_fgsea_lincs2, keep.all.vertices = FALSE)
?igraph::intersection
plot(nw_intersection, vertex.label = "")
nw_intersection

E(nw_intersection)

shared_drugs = intersect(V(nw_fgsea_cmap)$name[V(nw_fgsea_cmap)$type == "DRUG"],
          intersect(V(nw_fgsea_lincs1)$name[V(nw_fgsea_lincs1)$type == "DRUG"], 
                    V(nw_fgsea_lincs2)$name[V(nw_fgsea_lincs2)$type == "DRUG"])
          )

neighbors(graph = nw_fgsea_cmap, v = "HYDROCORTISONE")
