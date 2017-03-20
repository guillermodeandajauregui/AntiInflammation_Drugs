################################
#
#Antiinflammatory drugs
#in LINCS, CMap, and FAERS
#
################################

################################
#libraries
################################

library(stringr)
library(RSQLite)
library(data.table)
################################
#Input files
################################

inputs = read.table(file = "input_files.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

#manually curated anti inflammatory drug annotation from ATC
ai_drugs = inputs[1]
#annotation files for LINCS, phase 1
lincs_phase_1 = inputs[2]
#annotation files for LINCS, phase 2
lincs_phase_2 = inputs[3]
#annotation files for cMap
cMap = inputs[4]
#FAERS database
FAERS = inputs[5]

################################
#Read files
################################
drugs = fread(input = ai_drugs, data.table = FALSE, header = FALSE)
drugList = drugs$V2
glucocorticoids = drugs$V2[83:99]

lincs_phase_1 = fread(input = lincs_phase_1, data.table = FALSE, header = FALSE)
lincs_phase_2 = fread(input = lincs_phase_2, data.table = FALSE, header = TRUE)
lincs_phase_1[1,]
colnames(lincs_phase_1)<-lincs_phase_1[1,]
lincs_phase_1<-lincs_phase_1[-1,]
lincs_phase_1_drugs = unique(str_to_upper(unique(lincs_phase_1$pert_desc)))
lincs_phase_2_drugs = unique(str_to_upper(unique(lincs_phase_2$sm_name)))

cMap = read.table(file = cMap)
cMap_drugs = unique(cMap$normalized_name)

################################
#Connect to FAERS db
################################
con = dbConnect(drv=SQLite(), dbname=FAERS)
FAERS_drugs = unique(dbGetQuery(con, "SELECT * FROM drugmap")$replacement)

################################
#Match drug list to resources
################################

drugs_in_lincs_phase_1 = drugList[drugList%in%lincs_phase_1_drugs]
drugs_in_lincs_phase_2 = drugList[drugList%in%lincs_phase_2_drugs]
drugs_in_cMap = drugList[drugList%in%cMap]
drugs_in_FAERS = drugList[drugList%in%FAERS_drugs]

if("results"%in%dir()==FALSE){
  system("mkdir results")
}
setwd("results")
write.table(drugs_in_lincs_phase_1, file = "ai_drugs_in_lincs_phase_1.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(drugs_in_lincs_phase_2, file = "ai_drugs_in_lincs_phase_2.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(drugs_in_cMap, file = "ai_drugs_in_cMap.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(drugs_in_FAERS, file = "ai_drugs_in_FAERS.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

