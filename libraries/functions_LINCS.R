########################################################
#
#Functions for handling connectivity Map (LINCS) data. 
# 
########################################################

########################################################
#drugEset_lincs 
## Takes PATHS to LINCS level 5 gctx, signature info
## and gene info, plus a list of drugs (ALL CAPS).
## Returns an annotated ExpressionSet 
## with samples treated with drugs 
########################################################
drugEset_lincs_level5 <- function(gctx_path, 
                                  siginfo_path,
                                  geneinfo_path,
                                  DrugList)
{
  #read files 
  lincs_gctx = parse.gctx(gctx_path)
  print("loaded gctx")
  lincs_siginfo = fread(siginfo_path)
  print("loaded signature info")
  lincs_geneinfo = fread(geneinfo_path)
  print("loaded gene info")
  
  #find singatures that match druglist
  instances_ai = lincs_siginfo[which(str_to_upper(lincs_siginfo$pert_iname)%in%DrugList), 
                                c("sig_id", "pert_iname")]
  
  #subset matrix in gctx, and change row ids to gene symbol
  subset_matrix = lincs_gctx@mat[, match(instances_ai$sig_id, lincs_gctx@cdesc$id)]
  rownames(subset_matrix) <- lincs_geneinfo$pr_gene_symbol[match(rownames(subset_matrix), 
                                                           as.character(lincs_geneinfo$pr_gene_id))]
  
  #make object of class ExpressionSet
  
  eset = ExpressionSet(assayData = subset_matrix) #dummy out to force column order
  #Drug annotation for eset
  DrugData = as.data.frame(instances_ai[match(colnames(exprs(eset)), 
                                              instances_ai$sig_id),])
  rownames(DrugData) = DrugData$sig_id
  DrugData = DrugData[,-1, drop = FALSE]
  DrugData$pert_iname <- str_to_upper(DrugData$pert_iname)
  DrugAFD = new("AnnotatedDataFrame", data = DrugData)
  #annotated expression set
  eset = ExpressionSet(assayData = subset_matrix, phenoData = DrugAFD)
  return(eset)
}

########################################################
#drugEset_lincs_celltype 
## retrieves eset for particular cell lines only
########################################################


drugEset_lincs_l5_celltype <- function(gctx_path, 
                                  siginfo_path,
                                  geneinfo_path,
                                  DrugList,
                                  CellLines)
{
  #read files 
  lincs_gctx = parse.gctx(gctx_path)
  print("loaded gctx")
  lincs_siginfo = fread(siginfo_path)
  print("loaded signature info")
  lincs_geneinfo = fread(geneinfo_path)
  print("loaded gene info")
  
  #find singatures that match druglist
  instances_ai = lincs_siginfo[which(str_to_upper(lincs_siginfo$pert_iname)%in%DrugList), 
                               c("sig_id", "pert_iname", "cell_id")]
  
  instances_ai = instances_ai[which(instances_ai$cell_id%in%CellLines), 
                               c("sig_id", "pert_iname", "cell_id")]
  
  #subset matrix in gctx, and change row ids to gene symbol
  subset_matrix = lincs_gctx@mat[, match(instances_ai$sig_id, lincs_gctx@cdesc$id)]
  rownames(subset_matrix) <- lincs_geneinfo$pr_gene_symbol[match(rownames(subset_matrix), 
                                                                 as.character(lincs_geneinfo$pr_gene_id))]
  
  #make object of class ExpressionSet
  
  eset = ExpressionSet(assayData = subset_matrix) #dummy out to force column order
  #Drug annotation for eset
  DrugData = as.data.frame(instances_ai[match(colnames(exprs(eset)), 
                                              instances_ai$sig_id),])
  rownames(DrugData) = DrugData$sig_id
  DrugData = DrugData[,-1, drop = FALSE]
  DrugData$pert_iname <- str_to_upper(DrugData$pert_iname)
  DrugAFD = new("AnnotatedDataFrame", data = DrugData)
  #annotated expression set
  eset = ExpressionSet(assayData = subset_matrix, phenoData = DrugAFD)
  return(eset)
}
