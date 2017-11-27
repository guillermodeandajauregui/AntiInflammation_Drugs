#################################
#Anti-inflammatory Drug Network #
#################################

################################
#02 extracting ExpMatrices 
################################

#read list of drugs
#subset CMap/LINCS matrix
#aggregate to gene level 
#aggregate drugs through Kru-Bor
#write out matrix

################################
#libraries
################################
source("libraries/libraries.R")

#
#inputs = read.table(file = "input_files.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
inputs = scan(file = "input_files.txt", what = "character", comment.char = "#")
#read list of drugs
ai_cmap   =readLines("results/ai_drugs_in_cMap.txt")
ai_lincs1 =readLines("results/ai_drugs_in_lincs_phase_1.txt")
ai_lincs2 =readLines("results/ai_drugs_in_lincs_phase_2.txt")
ai_faers  =readLines("results/ai_drugs_in_FAERS.txt")
ai_lincs = union(ai_lincs1, ai_lincs2)

##select only drugs with adverse event info
ai_cmap_faers  = intersect(ai_cmap, ai_faers)
ai_lincs_faers = intersect(ai_lincs, ai_faers)
ai_lincs1_faers = intersect(ai_lincs1, ai_faers)
ai_lincs2_faers = intersect(ai_lincs2, ai_faers)

#subset CMap/LINCS matrix

###########################
#       CMap              #
###########################
print("processing CMap")

cmap_path       = inputs[6]
cmap_matrix     = ReRank(fread.Ranked.Matrix(cmap_path))
cmap_annotation = fread(input = inputs[4], data.table = FALSE)

#already aggregated at gene level

#subset CMap
ai_cmap_DrugEset = drug.eset(RankedMatrix = cmap_matrix, 
                             SampleInfoFile = cmap_annotation, 
                             DrugList = ai_cmap)

#merge all drug profiles by Kru-Bor
ai_cmap_krubor = RankMerging(ai_cmap_DrugEset)
exprs_ai_cmap_krubor = as.data.frame(exprs(ai_cmap_krubor))

fwrite(x = exprs_ai_cmap_krubor, 
       file = "results/krubor_cmap_ai.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)


###########################
#       LINCS             #
###########################

#phase 1
print("processing LINCS-Phase1")

lincs_phase1_path_gctx     = inputs[7]
lincs_phase1_path_siginfo  = inputs[9]
lincs_phase1_path_geneinfo = inputs[11]

ai_lincs1_eset  =    drugEset_lincs_level5(gctx_path = lincs_phase1_path_gctx,
                                        siginfo_path = lincs_phase1_path_siginfo,
                                        geneinfo_path = lincs_phase1_path_geneinfo,
                                        DrugList = ai_lincs1)

ai_lincs1_krubor = RankMerging(ai_lincs1_eset)

exprs_ai_lincs1_krubor = as.data.frame(exprs(ai_lincs1_krubor))

fwrite(x = exprs_ai_lincs1_krubor, 
       file = "results/krubor_lincs1_ai.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

#Phase 2
print("processing LINCS-Phase2")

lincs_phase2_path_gctx     = inputs[8]
lincs_phase2_path_siginfo  = inputs[10]
lincs_phase2_path_geneinfo = inputs[12]

ai_lincs2_eset  =    drugEset_lincs_level5(gctx_path = lincs_phase2_path_gctx,
                                           siginfo_path = lincs_phase2_path_siginfo,
                                           geneinfo_path = lincs_phase2_path_geneinfo,
                                           DrugList = ai_lincs2)

ai_lincs2_krubor = RankMerging(ai_lincs2_eset)

exprs_ai_lincs2_krubor = as.data.frame(exprs(ai_lincs2_krubor))

fwrite(x = exprs_ai_lincs2_krubor, 
       file = "results/krubor_lincs2_ai.txt", 
       quote = FALSE, 
       sep = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

save(ai_cmap_krubor, ai_lincs1_krubor, ai_lincs2_krubor, file = "results/krubors_for_mantra.RData")
print("check your matrices in the results directory")