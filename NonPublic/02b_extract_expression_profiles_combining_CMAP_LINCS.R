#################################
#Anti-inflammatory Drug Network #
#################################

################################
#02 extracting ExpMatrices 
#combining CMap and LINCS
################################

#extract cmap, lincs
#extract annotation files
#cbind
#combine annotation files
#eset
#krubor

source("libraries/libraries.R")
#
inputs = scan(file = "input_files.txt", what = "character", comment.char = "#")
#read list of drugs
ai_cmap   =readLines("results/ai_drugs_in_cMap.txt")
ai_lincs1 =readLines("results/ai_drugs_in_lincs_phase_1.txt")
ai_lincs2 =readLines("results/ai_drugs_in_lincs_phase_2.txt")

#CMap
cmap_path       = inputs[6]
cmap_matrix     = ReRank(fread.Ranked.Matrix(cmap_path))
cmap_annotation = fread(input = inputs[4], data.table = FALSE)

#make a subset of the RankedMatrix, containing only these Ids
ai_cmap_drugmatrix = cmap_matrix[,as.character(cmap_annotation$instance_id[which(cmap_annotation$normalized_name%in%ai_cmap)])]

DrugData_cmap = cmap_annotation[which(cmap_annotation$instance_id%in%colnames(ai_cmap_drugmatrix)), c(1, 4)]
colnames(DrugData_cmap)[2] <- "DRUG"
colnames(DrugData_cmap)[1] <- "sig_id"
#LINCS_1
lincs_phase1_path_gctx     = inputs[7]
lincs_phase1_path_siginfo  = inputs[9]
lincs_phase1_path_geneinfo = inputs[11]

lincs1_gctx = parse.gctx(lincs_phase1_path_gctx)
print("loaded gctx")
lincs1_siginfo = fread(lincs_phase1_path_siginfo)
print("loaded signature info")
lincs1_geneinfo = fread(lincs_phase1_path_geneinfo)
print("loaded gene info")

instances_ai = lincs1_siginfo[which(str_to_upper(lincs1_siginfo$pert_iname)%in%ai_lincs1), 
                             c("sig_id", "pert_iname")]

#subset matrix in gctx, and change row ids to gene symbol
subset_matrix_lincs1 = lincs1_gctx@mat[, match(instances_ai$sig_id, lincs1_gctx@cdesc$id)]
rownames(subset_matrix_lincs1) <- lincs1_geneinfo$pr_gene_symbol[match(rownames(subset_matrix_lincs1), 
                                                               as.character(lincs1_geneinfo$pr_gene_id))]

DrugData_lincs1 = as.data.frame(instances_ai[match(colnames(subset_matrix_lincs1), 
                                                   instances_ai$sig_id),])

DrugData_lincs1$pert_iname = str_to_upper(DrugData_lincs1$pert_iname)
colnames(DrugData_lincs1)[2] = "DRUG"

#LINCS_2
lincs_phase2_path_gctx     = inputs[8]
lincs_phase2_path_siginfo  = inputs[10]
lincs_phase2_path_geneinfo = inputs[12]

lincs2_gctx = parse.gctx(lincs_phase2_path_gctx)
print("loaded gctx")
lincs2_siginfo = fread(lincs_phase2_path_siginfo)
print("loaded signature info")
lincs2_geneinfo = fread(lincs_phase2_path_geneinfo)
print("loaded gene info")

instances_ai = lincs2_siginfo[which(str_to_upper(lincs2_siginfo$pert_iname)%in%ai_lincs2), 
                              c("sig_id", "pert_iname")]

#subset matrix in gctx, and change row ids to gene symbol
subset_matrix_lincs2 = lincs2_gctx@mat[, match(instances_ai$sig_id, lincs2_gctx@cdesc$id)]
rownames(subset_matrix_lincs2) <- lincs2_geneinfo$pr_gene_symbol[match(rownames(subset_matrix_lincs2), 
                                                                       as.character(lincs2_geneinfo$pr_gene_id))]

DrugData_lincs2 = as.data.frame(instances_ai[match(colnames(subset_matrix_lincs2), 
                                                   instances_ai$sig_id),])

DrugData_lincs2$pert_iname = str_to_upper(DrugData_lincs2$pert_iname)
colnames(DrugData_lincs2)[2] = "DRUG"

#cbind the three matrices

#remove non genes that are not in CMap and LINCS, and reorder matrices
ai_cmap_drugmatrix_matching_genes = ai_cmap_drugmatrix[match(intersect(rownames(ai_cmap_drugmatrix), 
                                                                       rownames(subset_matrix_lincs1)), 
                                                             rownames(ai_cmap_drugmatrix)
),]

subset_matrix_lincs1_matching_genes = subset_matrix_lincs1[match(intersect(rownames(ai_cmap_drugmatrix), 
                                                                           rownames(subset_matrix_lincs1)), 
                                                             rownames(subset_matrix_lincs1)
),]

subset_matrix_lincs2_matching_genes = subset_matrix_lincs2[match(intersect(rownames(ai_cmap_drugmatrix), 
                                                                           rownames(subset_matrix_lincs1)), 
                                                               rownames(subset_matrix_lincs2)
),]


combined_matrices = cbind(ai_cmap_drugmatrix_matching_genes, subset_matrix_lincs1_matching_genes, subset_matrix_lincs2_matching_genes)

combined_matrices[1:5, 1:5]

#combine the annotation files

drugdata_total = rbind(DrugData_cmap, DrugData_lincs1, DrugData_lincs2)
head(drugdata_total)

#make the eset
eset = ExpressionSet(assayData = combined_matrices) #dummy out to force column order
  #Drug annotation for eset
DrugData = as.data.frame(drugdata_total[match(colnames(exprs(eset)), 
                                              drugdata_total$sig_id),])
rownames(DrugData) = DrugData$sig_id
DrugData = DrugData[,-1, drop = FALSE]
DrugData$DRUG <- str_to_upper(DrugData$DRUG)
DrugAFD = new("AnnotatedDataFrame", data = DrugData)
#annotated expression set
eset = ExpressionSet(assayData = combined_matrices, phenoData = DrugAFD)

#KruBor
combined_krubor = RankMerging(eset)
exprs_combined_krubor = as.data.frame(exprs(combined_krubor))

fwrite(x = exprs_combined_krubor, 
       file = "results/krubor_combined_cmap_lincs_ai.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)
