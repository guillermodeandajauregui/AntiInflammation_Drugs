################################################################
#Anti-inflammatory Drug Network                                #
################################################################

################################################################
#05 Drug- Adverse Drug Reaction Network
################################################################

##################################
#For each anti-inflammatory drug
#identify Adverse Drug Reactions
#(ADRs)
#With the following conditions:
#1) FDR =< 0.05 for PRR vs other antiinflammatories
#2) FDR =< 0.05 for PRR All_antiinflammatories vs Non-antiinflammatories
#express as network
##################################

#read normalized FAERS data 
#read anti-inflammatory drug list
#calculate Proportional Reporting Ratios:
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)
#calculate PRR (anti-inflammatories vs anti-inflammatories)
#construct network

##################################
#Libraries
##################################
library(data.table)
library(PhViD)
library(igraph)
##################################
#Inputs
##################################

inputs = scan("input_files.txt", what = "character", quiet = TRUE)

#read normalized FAERS data 
faers = fread(inputs[13])
  #FAERS normalization contains multimapping compounds; remove those
faers_sans4 = faers[V10 != "Type4"]

#read anti-inflammatory drug list
antiinflammatories = fread(inputs[1], header = FALSE)
aif                = antiinflammatories$V2

################################################################
#PRR CALCULATIONS
################################################################

#################
#01 Anti-inflammatories v Anti - inflammatories
#################

#calculate PRR (all anti-inflammatories vs all anti-inflammatories)

dt_aif <- copy(faers_sans4)
dt_aif <-dt_aif[V14%in%aif]

  #make a frequency count of drug - adverse effect 
freq_aif = dt_aif[, .(.N), by=list(V14, V9)]
  #make an object of class PhViD for calculations
phvid_aif = as.PhViD(freq_aif)

  #calculate Proportional Reporting Ratio
res_aif_aif <- PRR(phvid_aif)

#################
#02 All Anti-inflammatories v non Anti-inflammatories
#################
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)

#convert all anti-inflammatories to a single "ALL_ANTIINFLAMMATORIES"
#convert all NON anti-inflammatories to a single "NON"

dt_dosvar <- copy(faers_sans4)
dt_dosvar[!(V14%in%aif), V14:="NON"]
dt_dosvar[V14%in%aif, V14:="ALL_ANTIINFLAMMATORIES"]

#make a frequency count of drug - adverse effect 
freq_dosvar = dt_dosvar[, .(.N), by=list(V14, V9)]
#make an object of class PhViD for calculations
phvid_dosvar = as.PhViD(freq_dosvar)
#calculate Proportional Reporting Ratio
res_dosvar <- PRR(phvid_dosvar)

################################################################
#Network integration
################################################################

#adverse effects identified in DOSVAR
ADR_dosvar = unique(res_dosvar$SIGNALS$`event effect`) #7820
ADR_dosvar_aif = unique(res_dosvar$SIGNALS$`event effect`[res_dosvar$SIGNALS$`drug code` == "ALL_ANTIINFLAMMATORIES"]) #4287

#adverse effects identified in aif v aif
ADR_aif_aif = unique(res_aif_aif$SIGNALS$`event effect`) #7841

#intersections of aif_aif and aif_v_non_aif
ADR_aifaif_2aif = intersect(ADR_aif_aif, ADR_dosvar_aif)
df_c =res_aif_aif$SIGNALS[res_aif_aif$SIGNALS$`event effect`%in%ADR_aifaif_2aif,] #nrow 10235

g_c = graph_from_data_frame(df_c, directed = TRUE)

write.graph(g_c, file = "results/drug_adr.gml", format = "gml")
