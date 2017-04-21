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
library(data.table)
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

#subset CMap/LINCS matrix

###########################
#       CMap              #
###########################

cmap_path       = inputs[6]
cmap_matrix     = ReRank(fread.Ranked.Matrix(cmap_path))
cmap_annotation = fread(input = inputs[4], data.table = FALSE)

#already aggregated at gene level

#subset CMap
ai_cmap_DrugEset = drug.eset(RankedMatrix = cmap_matrix, 
                             SampleInfoFile = cmap_annotation, 
                             DrugList = ai_cmap_faers)

#merge all drug profiles by Kru-Bor
ai_cmap_krubor = RankMerging(ai_cmap_DrugEset)
exprs_ai_cmap_krubor = exprs(ai_cmap_krubor)


###########################
#       LINCS             #
###########################

lincs_phase1_anot = inputs[9]
lincs_phase2_anot = inputs[10]

lincs_phase1_matrix = inputs[7]
lincs_phase2_matrix = inputs[8]

##match column annotation to drug list
##subset matrix 
##make annotated eset
##kru-bor