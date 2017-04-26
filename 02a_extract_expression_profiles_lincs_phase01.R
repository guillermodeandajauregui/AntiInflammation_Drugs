#######################################
#THIS FILE WAS USED TO TEST FUNCTION
#IT IS NOT TRACKED IN GIT
#AND THE KRUBORS GENERATED WITH IT
#CONSIDER DRUGLISTS BOTH IN 
# LINCS AND FAERS
#######################################
drugEset_lincs_phase01 <-function(lincs_phase01, 
                                  lincs_annotation, 
                                  geneinfo, 
                                  drugList){
  #lincs_phase01 = path to gctx, level 5
  #drugList = list of drug Names 
  #lincs_annotation = signature info, level 5 #PATH
  #geneinfo = lincs phase01 gene info #PATH

  #1) read gctx
  bigx = parse.gctx(fname = lincs_phase01)
  print("finished opening gctx")
  #2) read signature annotation
  instances = fread(input = lincs_annotation, 
                    data.table = FALSE, 
                    header = TRUE)
  
  #2a) Read geneinfo
  geneinfo = fread(input = geneinfo, 
                   data.table = FALSE, 
                   header = TRUE)
  print("finished opening annotations")
  #3) id signatures de interes  - match w/ drugs
  instances_ai = instances[which(str_to_upper(instances$pert_iname)%in%drugList), 
                           c("sig_id", "pert_iname")]
  
  #4) subset gctx@mat con ids de interes
  subset_matrix = bigx@mat[, match(instances_ai$sig_id, bigx@cdesc$id)]
  
  rownames(subset_matrix) <- geneinfo$pr_gene_symbol[match(rownames(subset_matrix), 
                                                                                  as.character(geneinfo$pr_gene_id))]
  
  #6) make exprset with gctx@mat[genes,ids]
  eset = ExpressionSet(assayData = subset_matrix)
  #DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 4)] #id, drugname
  DrugData = instances_ai[match(colnames(exprs(eset)), 
                                instances_ai$sig_id),]
  rownames(DrugData) = DrugData[,1]
  DrugData = DrugData[,-1, drop = FALSE]
  DrugData$pert_iname <- str_to_upper(DrugData$pert_iname)
  DrugAFD = new("AnnotatedDataFrame", data = DrugData)
  eset = ExpressionSet(assayData = subset_matrix, phenoData = DrugAFD)
  
  return(eset)
}
#########################################################################
drugEset_lincs_phase02 <- function(gctx_path, 
                                   siginfo_path,
                                   geneinfo_path,
                                   DrugList)
  {
  #read files 
  lincs2_gctx = parse.gctx(lincs2_gctx)
    print("loaded gctx")
  lincs2_siginfo = fread(lincs2_siginfo)
    print("loaded signature info")
  lincs2_geneinfo = fread(lincs2_geneinfo)
    print("loaded gene info")
    
  #find singatures that match druglist
  instances_ai = lincs2_siginfo[which(str_to_upper(lincs2_siginfo$pert_iname)%in%DrugList), 
                           c("sig_id", "pert_iname")]
  
  #subset matrix in gctx, and change row ids to gene symbol
  subset_matrix = lincs2_gctx@mat[, match(instances_ai$sig_id, lincs2_gctx@cdesc$id)]
  rownames(subset_matrix) <- geneinfo$pr_gene_symbol[match(rownames(subset_matrix), 
                                                           as.character(geneinfo$pr_gene_id))]
  
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

lincs2_gctx     = inputs[8]
lincs2_siginfo  = inputs[10]
lincs2_geneinfo = inputs[12]

lincs2_gctx = parse.gctx(lincs2_gctx)
lincs2_siginfo = fread(lincs2_siginfo)
lincs2_geneinfo = fread(lincs2_geneinfo)

length(colnames(lincs2_gctx@mat))

drugList = ai_lincs2_faers
instances = lincs2_siginfo
geneinfo = lincs2_geneinfo

instances_ai = instances[which(str_to_upper(instances$pert_iname)%in%drugList), 
                         c("sig_id", "pert_iname")]

subset_matrix = lincs2_gctx@mat[, match(instances_ai$sig_id, lincs2_gctx@cdesc$id)]
rownames(subset_matrix) <- geneinfo$pr_gene_symbol[match(rownames(subset_matrix), 
                                                         as.character(geneinfo$pr_gene_id))]


eset = ExpressionSet(assayData = subset_matrix)
#DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 4)] #id, drugname

DrugData = as.data.frame(instances_ai[match(colnames(exprs(eset)), 
                              instances_ai$sig_id),])
rownames(DrugData) = DrugData$sig_id
DrugData = DrugData[,-1, drop = FALSE]

DrugData$pert_iname <- str_to_upper(DrugData$pert_iname)
DrugAFD = new("AnnotatedDataFrame", data = DrugData)
eset = ExpressionSet(assayData = subset_matrix, phenoData = DrugAFD)

krubor_lincs2 = RankMerging(eset)
fwrite(x = as.data.frame(exprs(krubor_lincs2)),
                         file = "results/krubor_lincs2_ai.txt",
                         row.names = TRUE, 
                         col.names = TRUE, 
                         sep ="\t", 
                         quote = FALSE
       )
                        
# #######################################################################
# lincs_annotation = fread(input = inputs[9], 
#                          data.table = FALSE, 
#                          header = TRUE)
# 
# geneinfo = fread(input = inputs[11], 
#                  data.table = FALSE, 
#                  header = TRUE)
# 
# 
# 
# #######################################################################
# #1) read gctx
# bigx = parse.gctx(fname = inputs[7])
# #2) read anotacion de signatures
# instances = fread(input = inputs[9], data.table = FALSE, header = TRUE)
# #3) id signatures de interes  - match w/ drugs
# instances_ai = instances[which(str_to_upper(instances$pert_iname)%in%ai_lincs), 
#                 c("sig_id", "pert_iname")]
# 
# #4) subset gctx@mat con ids de interes
#   #match(instances$sig_id, bigx@cdesc$id)[1:100]
# lincs_phase01_ai_matrix = bigx@mat[, match(instances_ai$sig_id, bigx@cdesc$id)]
# #5) no need to aggregate genes; 1 row per gene
# rownames(lincs_phase01_ai_matrix) <- geneinfo_lincs_phase1$pr_gene_symbol[match(rownames(lincs_phase01_ai_matrix), 
#                                                                                 as.character(geneinfo_lincs_phase1$pr_gene_id))]
# #6) make exprset with gctx@mat[genes,ids]
# eset = ExpressionSet(assayData = lincs_phase01_ai_matrix)
# #DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 4)] #id, drugname
# DrugData = instances_ai[match(colnames(exprs(eset)), instances_ai$sig_id),]
# rownames(DrugData) = DrugData[,1]
# DrugData = DrugData[,-1, drop = FALSE]
# DrugData$pert_iname <- str_to_upper(DrugData$pert_iname)
# DrugAFD = new("AnnotatedDataFrame", data = DrugData)
# eset = ExpressionSet(assayData = lincs_phase01_ai_matrix, phenoData = DrugAFD)
# #7) agregar tx al exprset 
# 
