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

lincs_matrix= fread(input = "results/krubor_combined_lincs_1_and_2_ai.txt", data.table = FALSE)
rownames(lincs_matrix) <- lincs_matrix[,1]
lincs_matrix <- lincs_matrix[,-1]

#fgsea enrichment 

PWs_all = PWs
PWs <- PWs[which(sapply(PWs, length)>10)]

fgsea_lincs_combined = fgsea_expmatrix(expmatrix = lincs_matrix, 
                               pathway_list = PWs, 
                               nperm = 50000, 
                               maxSize = 500, 
                               padj_cutoff = 0.1)

#generate networks

nw_fgsea_lincs_combined = graph_from_fgsea(fgsea_lincs_combined)

save(fgsea_lincs_combined, nw_fgsea_lincs_combined,
     file = "results/fgsea_nws_lincs_combined.RData")

#