#################################
#Anti-inflammatory Drug Network #
#################################

################################
#03 Drug - Drug Network
################################

#read Ranked Expression Matrix
#for each drug
#Spearman footrule to each other
#return spearman matrix
source("libraries/libraries.R")

#read ranked matrices 



cmap_matrix_a= fread(input = "results/krubor_cmap_ai.txt", data.table = FALSE)
rownames(cmap_matrix_a) <- cmap_matrix_a[,1]
cmap_matrix_a <- cmap_matrix_a[,-1]
cmap_spearman_matrix = spearman_footrule_matrix(cmap_matrix_a)

lincs1_matrix_a= fread(input = "results/krubor_lincs1_ai.txt", data.table = FALSE)
rownames(lincs1_matrix_a) <- lincs1_matrix_a[,1]
lincs1_matrix_a <- lincs1_matrix_a[,-1]
lincs1_spearman_matrix = spearman_footrule_matrix(lincs1_matrix_a)

lincs2_matrix_a= fread(input = "results/krubor_lincs2_ai.txt", data.table = FALSE)
rownames(lincs2_matrix_a) <- lincs2_matrix_a[,1]
lincs2_matrix_a <- lincs2_matrix_a[,-1]
lincs2_spearman_matrix = spearman_footrule_matrix(lincs2_matrix_a)

fwrite(x = as.data.frame(cmap_spearman_matrix), 
       file = "results/cmap_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs1_spearman_matrix), 
       file = "results/lincs1_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(x = as.data.frame(lincs2_spearman_matrix), 
       file = "results/lincs2_spearman_footrule_matrix.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)
